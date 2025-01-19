#include "mppiEstimation.h"

STATE::STATE() {
    p.setZero();
    R.setIdentity();
    v.setZero();
}

mppiEstimation::mppiEstimation(): T(10), dimU(6) {
    anchorPositions.setZero();
    N = 10;
    _g << 0, 0, 9.81;
    TOL = 1e-9;
    dt = 0;
    sigmaUvector.resize(6);
    sigmaUvector << 0, 0, 0, 0, 0, 0;
    sigmaU = sigmaUvector.asDiagonal();
    gammaU = 0;
    resultPuber = nh.advertise<geometry_msgs::PoseStamped>("mppi_pose", 1);
    u0 = Eigen::VectorXd::Zero(dimU);
    U0 = Eigen::MatrixXd::Zero(dimU, T);
    Xo.resize(T+1);
}

void mppiEstimation::setDt(const double deltaT) {
    dt = deltaT;
}

void mppiEstimation::setAnchorPositions(const Eigen::Matrix<double, 3, 8> &positions) {
    anchorPositions = positions;
    std::cout <<"Anchor positions: \n"<<anchorPositions<<std::endl;
}

STATE mppiEstimation::f(const STATE &state, const ImuData &imuData, const Eigen::VectorXd Ui) {
    Eigen::Matrix3d Rot = state.R;
    STATE next = state;
    Eigen::Vector3d accWorld = Rot*(imuData.acc - Ui.segment(0,3)) + _g;
    next.p += state.v*dt+0.5*accWorld*dt*dt;
    next.v += accWorld*dt;
    next.R = Rot*Exp((imuData.gyr - Ui.segment(3,3))*dt);
    return next;
}

Eigen::MatrixXd mppiEstimation::getNoise(const int &T) {
    return sigmaU * normGen.template generate<Eigen::MatrixXd>(dimU, T, urng);
}

void mppiEstimation::move(const ImuData &imuData) {
    xInit = f(xInit, imuData, u0);
    U0.leftCols(T-1) = Uo.rightCols(T-1); 
}

void mppiEstimation::solve(const std::vector<UwbData> &uwbData, const std::vector<ImuData> &imuData) {
    Eigen::MatrixXd Ui = U0.replicate(N, 1);
    Eigen::VectorXd costs(N);
    Eigen::VectorXd weights(N);

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        STATE Xi[T+1];
        Eigen::MatrixXd noise = getNoise(T);
        Ui.middleRows(i * dimU, dimU) += noise;

        Xi[0] = xInit;
        double cost = 0.0;

        for (int j = 0; j < T; ++j) {

            Xi[j+1] = f(Xi[j], imuData[j], Ui.block(i * dimU, j, dimU, 1));
            Eigen::VectorXd Hx(8);
            for (int k = 0; k < 8; ++k) {
                double dx = Xi[j].p(0) - anchorPositions(0,k);
                double dy = Xi[j].p(1) - anchorPositions(1,k);
                double dz = Xi[j].p(2) - anchorPositions(2,k);
                Hx(k) = std::sqrt(dx*dx + dy*dy + dz*dz);
            }

            double stepCost = (uwbData[j].ranges - Hx).squaredNorm();
            cost += stepCost;
        }
        costs(i) = cost;
    }

    double minCost = costs.minCoeff();
    weights = (-gammaU * (costs.array() - minCost)).exp();
    double totalWeight = weights.sum();
    weights /= totalWeight;

    Uo = Eigen::MatrixXd::Zero(dimU, T);
    for (int i = 0; i < N; ++i) {
        Uo += Ui.middleRows(i * dimU, dimU) * weights(i);
    }

    u0 = Uo.col(0);

    Xo[0] = xInit;
    for (int j = 0; j < T; ++j) {
        Xo[j+1] = f(Xo[j], imuData[j], Uo.col(j));
    }

    // visual_traj.push_back(xInit);

    publishPose(xInit);

    // publishPose(Xo);

}

void mppiEstimation::publishPose(const STATE &state) {
    geometry_msgs::PoseStamped pose;
    pose.header.frame_id = "map";
    pose.header.stamp = ros::Time::now();

    pose.pose.position.x = state.p(0);
    pose.pose.position.y = state.p(1);
    pose.pose.position.z = state.p(2);
    Eigen::Quaterniond q(state.R);
    pose.pose.orientation.w = q.w();
    pose.pose.orientation.x = q.x();
    pose.pose.orientation.y = q.y();
    pose.pose.orientation.z = q.z();

    resultPuber.publish(pose);
}

Eigen::MatrixXd mppiEstimation::Exp(const Eigen::Vector3d &omega) {
    double angle = omega.norm();
    Eigen::Matrix3d Rot;
    
    if (angle<TOL){
        Rot = Eigen::Matrix3d::Identity();
    }
    else{
        Eigen::Vector3d axis = omega/angle;
        double c = cos(angle);
        double s = sin(angle);

        Rot = c*Eigen::Matrix3d::Identity() + (1 - c)*axis*axis.transpose() + s*vectorToSkewSymmetric(axis);
    }
    return Rot;    
}

Eigen::MatrixXd mppiEstimation::vectorToSkewSymmetric(const Eigen::Vector3d &vector) {
    Eigen::Matrix3d Rot;
    Rot << 0, -vector.z(), vector.y(),
          vector.z(), 0, -vector.x(),
          -vector.y(), vector.x(), 0;
    
    return Rot;
}

mppiEstimation::~mppiEstimation() {}

Node::Node() {
    uwbSub = nh_.subscribe("/nlink_linktrack_tagframe0", 10, &Node::uwbCallback, this);
    imuSub = nh_.subscribe("/imu/data", 10, &Node::imuCallback, this);
    anchorPositions <<  0,    0,    8.86, 8.86, 0,    0,    8.86,  8.86,
            0,    8.00, 8.00, 0,    0,    8.00,  8.00,  0,
            0.0,  0.0,  0.0,  0.0,  2.20, 2.20,  2.20,  2.20;
    MppiEstimation.setAnchorPositions(anchorPositions);
    uwbData.ranges = Eigen::VectorXd(8);
}

void Node::uwbCallback(const nlink_parser::LinktrackTagframe0 &msg) {
    if(!uwbInit){
        uwbInitTime = msg.system_time/1000.00;
        uwbInit = true;
    }
    uwbData.timeStamp = msg.system_time/1000.00  - uwbInitTime;  
    for (int i = 0; i < 8; i++) {
        uwbData.ranges[i] = msg.dis_arr[i];
    }
    uwbDataQueue.push(uwbData);
    run();
}

void Node::imuCallback(const sensor_msgs::ImuConstPtr &msg) {   
    if (!imuInit){
        imuInitTime = msg->header.stamp.toSec();
        imuInit = true;
    }
    imuData.timeStamp = msg->header.stamp.toSec()-imuInitTime;  
    imuData.gyr(0) = msg->angular_velocity.x;
    imuData.gyr(1) = msg->angular_velocity.y;
    imuData.gyr(2) = msg->angular_velocity.z;
    imuData.acc(0) = msg->linear_acceleration.x;
    imuData.acc(1) = msg->linear_acceleration.y;
    imuData.acc(2) = msg->linear_acceleration.z;

    imuDataQueue.push(imuData);
    run();
}

std::pair<std::vector<ImuData>, std::vector<UwbData>> Node::interpolationAllT()
{
    int tCount = MppiEstimation.T;

    if (imuDataQueue.empty() || uwbDataQueue.empty()) {
        return {{}, {}};
    }

    std::vector<ImuData> imuVec(tCount);
    std::vector<UwbData> uwbVec(tCount);

    for (int i = 0; i < tCount; ++i) {
     
        if (imuDataQueue.empty() || uwbDataQueue.empty()) {
            return {{}, {}};
        }

        while (!uwbDataQueue.empty() &&
               (uwbDataQueue.front().timeStamp < imuDataQueue.front().timeStamp ||
                uwbDataQueue.front().timeStamp > imuDataQueue.back().timeStamp)) {
            uwbDataQueue.pop();

            if (uwbDataQueue.empty() || imuDataQueue.empty()) {
                return {{}, {}};
            }
        }

            if (imuDataQueue.empty() || uwbDataQueue.empty()) {
                return {{}, {}};
            }
            
        ImuData imuData1 = imuDataQueue.front();
        ImuData imuData2 = imuDataQueue.back();
        UwbData uwbData  = uwbDataQueue.front();

        if (imuData1.timeStamp <= uwbData.timeStamp &&
            uwbData.timeStamp <= imuData2.timeStamp)
        {
            double t1 = imuData1.timeStamp;
            double t2 = imuData2.timeStamp;
            double t  = uwbData.timeStamp;
            double delta_t = t2 - t1;
            double alpha = 0.0;

            if (delta_t == 0.0) {
                ROS_WARN("IMU timeStamps t1 and t2 are equal. alpha=0 forced.");
            } else {
                alpha = (t - t1) / delta_t;
            }

            imuData1.acc = (1.0 - alpha) * imuData1.acc + alpha * imuData2.acc;
            imuData1.gyr = (1.0 - alpha) * imuData1.gyr + alpha * imuData2.gyr;

            dt = t - beforeT;
            beforeT = t;

            uwbDataQueue.pop();
            imuDataQueue.pop();

            imuVec[i] = imuData1;
            uwbVec[i] = uwbData;
        }
    }

    return std::make_pair(imuVec, uwbVec);
}

void Node::run() {
    auto [tImu, tUwb] = interpolationAllT();

    if (tImu.size() == MppiEstimation.T
        && tUwb.size() == MppiEstimation.T)
    {
        MppiEstimation.setDt(dt);
        MppiEstimation.solve(tUwb, tImu);
        MppiEstimation.move(tImu[MppiEstimation.T - 1]);
    }
}