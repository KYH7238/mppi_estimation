#include "mppiEstimation_2tag.h"

std::vector<UwbData> uwbBufferLeft;
std::vector<UwbData> uwbBufferRight;

Eigen::VectorXd make16UwbRanges(double imuTime) {
    auto leftIt = std::min_element(uwbBufferLeft.begin(), uwbBufferLeft.end(),
        [&](const UwbData& a, const UwbData& b) {
            return std::abs(a.timeStamp - imuTime) < std::abs(b.timeStamp - imuTime);
        });
    auto rightIt = std::min_element(uwbBufferRight.begin(), uwbBufferRight.end(),
        [&](const UwbData& a, const UwbData& b) {
            return std::abs(a.timeStamp - imuTime) < std::abs(b.timeStamp - imuTime);
        });

    Eigen::VectorXd ranges(16);
    if (leftIt != uwbBufferLeft.end() && rightIt != uwbBufferRight.end()) {
        ranges.head(8) = leftIt->ranges;
        ranges.tail(8) = rightIt->ranges;
    } else {
        ranges.setZero();
    }
    return ranges;
}

STATE::STATE() {
    p << 0.0, 0.0, 0.0;
    R.setIdentity();
    v.setZero();
}

mppiEstimation::mppiEstimation(): T(5), dimU(6) {
    anchorPositions.setZero();
    N = 12000;
    _g << 0, 0, 9.81;
    TOL = 1e-9;
    dt = 0;
    sigmaUvector.resize(6);
    float linSigma = 10.0; // 10
    float angSigma = 1; // 1 
    sigmaUvector << linSigma, linSigma, linSigma, angSigma, angSigma, angSigma;
    sigmaU = sigmaUvector.asDiagonal();
    gammaU = 10.0;
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
    // start = std::chrono::high_resolution_clock::now();
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
            Eigen::Vector3d tagL = getTagPosition(Xi[j], 0.13);
            // Eigen::Vector3d tagR = getTagPosition(Xi[j], -0.13);
            Eigen::VectorXd HxL = (anchorPositions.colwise() - tagL).colwise().norm();
            // Eigen::VectorXd HxR = (anchorPositions.colwise() - tagR).colwise().norm();
            Eigen::VectorXd Hx(8);
            Hx << HxL;
            double stepCost = (uwbData[j].ranges.head(8) - Hx).norm();
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

    finish = std::chrono::high_resolution_clock::now();
    // elapsed_1 = finish - start;
    // elapsed = elapsed_1.count();
    u0 = Uo.col(0);

    Xo[0] = xInit;
    for (int j = 0; j < T; ++j) {
        Xo[j+1] = f(Xo[j], imuData[j], Uo.col(j));
    }
    publishPose(xInit);
    // publishPose(Xo[T]);
}

void mppiEstimation::publishPose(const STATE &state) {
    geometry_msgs::PoseStamped pose;
    pose.header.frame_id = "map";
    pose.header.stamp = ros::Time::now();

    pose.pose.position.x = state.p(0);
    pose.pose.position.y = state.p(1);
    pose.pose.position.z = state.p(2);
    Eigen::Quaterniond q(state.R);
    pose.pose.orientation.x = q.x();
    pose.pose.orientation.y = q.y();
    pose.pose.orientation.z = q.z();
    pose.pose.orientation.w = q.w();
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
    imuSub = nh_.subscribe("/vectornav_driver_node/imu/data", 10, &Node::imuCallback, this);
    uwbSub = nh_.subscribe("/nlink_linktrack_tagframe0", 10, &Node::uwbCallback, this);
    anchorPositions <<  -4.4, -4.4,  4.43,  4.43, -4.4, -4.4, 4.43,  4.43,    
                        -4.04, 4.76, 4.89, -4.08, -4.08, 4.8, 4.93, -4.08,
                         0.26, 0.26, 0.3,   0.34,  2.45, 2.5, 2.47,  2.45;
    MppiEstimation.setAnchorPositions(anchorPositions);
    imuInit = uwbInit = false;
    lastStamp = 0;
    std::thread process(&Node::processThread, this);
    ros::MultiThreadedSpinner spinner(2);
    spinner.spin();
    process.join();
}

void Node::imuCallback(const sensor_msgs::ImuConstPtr &msg) {
    std::lock_guard<std::mutex> lock(imuMutex);
    double t = msg->header.stamp.toSec();
    if (!imuInit) {
        imuInitTime = t;
        imuInit = true;
        return;
    }
    ImuData d;
    d.timeStamp = t - imuInitTime;
    d.gyr << msg->angular_velocity.x,
             msg->angular_velocity.y,
             msg->angular_velocity.z;
    d.acc << msg->linear_acceleration.x,
             msg->linear_acceleration.y,
             msg->linear_acceleration.z;
    imuBuffer.push_back(d);
    dataCond.notify_one();
}

void Node::uwbCallback(const nlink_parser::LinktrackTagframe0 &msg) {
    std::lock_guard<std::mutex> lock(uwbMutex);
    double t = msg.system_time/1000.0;
    if (!uwbInit) {
        uwbInitTime = t;
        uwbInit = true;
    }
    UwbData d;
    d.timeStamp = t - uwbInitTime;
    d.ranges = Eigen::VectorXd(8);
    for (int i = 0; i < 8; i++) d.ranges(i) = msg.dis_arr[i];
    if (msg.id == 0) uwbBufferLeft.push_back(d);
    else if (msg.id == 1) uwbBufferRight.push_back(d);
    dataCond.notify_one();
}

std::pair<std::vector<ImuData>, std::vector<UwbData>> Node::interpolationAllT_blocking()
{
    ros::Rate loopRate(uwbFreq);
    int tCount = MppiEstimation.T;
    std::vector<ImuData> imuVec;
    std::vector<UwbData> uwbVec;
    imuVec.reserve(tCount);
    uwbVec.reserve(tCount);

    size_t idx = 0;
    while (uwbVec.size() < tCount) {
        {
            std::lock_guard<std::mutex> imuLock(imuMutex); // vs unique_lock : lock_guard는 자동으로 락을 해제하지만 unique_lock은 수동으로 락을 해제해야 함 
            std::lock_guard<std::mutex> uwbLock(uwbMutex);// 중간에 해제하고 다시 잠금해야하면 unique_lock 사용하고 그냥 한 thread동안 잠구기만 하면 lock_guard 사용(더 가벼움)
            if (imuBuffer.size() < 2 || uwbBufferLeft.size() <= idx || uwbBufferRight.size() <= idx) { // UWB Data가 안 들어왔으면 다시 대기 
                continue;
            }
        }
        loopRate.sleep();
        std::unique_lock<std::mutex> imuLock(imuMutex);
        std::unique_lock<std::mutex> uwbLock(uwbMutex);
        // 동기화: left/right UWB의 timestamp가 같다고 가정
        double curStamp = uwbBufferLeft[idx].timeStamp;
        if (uwbBufferRight[idx].timeStamp != curStamp) {
            // 만약 timestamp가 다르면, 더 작은 쪽 idx만 증가시키고 continue
            if (uwbBufferRight[idx].timeStamp < curStamp) { ++idx; continue; }
            else { ++idx; continue; }
        }
        bool found = false;
        ImuData interImuMsg;
        for (size_t i = 0; i + 1 < imuBuffer.size(); ++i) {
            if (imuBuffer[i].timeStamp <= curStamp && curStamp <= imuBuffer[i+1].timeStamp) {
                interpolateImuData(imuBuffer[i], imuBuffer[i+1], curStamp, interImuMsg);
                interImuMsg.timeStamp = curStamp;
                found = true;
                break;
            }
        }
        if (!found) { ++idx; continue; }
        UwbData merged;
        merged.timeStamp = curStamp;
        merged.ranges = Eigen::VectorXd(16);
        merged.ranges.head(8) = uwbBufferLeft[idx].ranges;
        merged.ranges.tail(8) = uwbBufferRight[idx].ranges;
        imuVec.push_back(interImuMsg);
        uwbVec.push_back(merged);
        ++idx;
    }
    if (idx > 0) {
        uwbBufferLeft.erase(uwbBufferLeft.begin(), uwbBufferLeft.begin() + idx);
        uwbBufferRight.erase(uwbBufferRight.begin(), uwbBufferRight.begin() + idx);
    }
    if (!imuVec.empty()) {
        double minUsedTime = imuVec.front().timeStamp;
        imuBuffer.erase(
            std::remove_if(imuBuffer.begin(), imuBuffer.end(),
                [&](const ImuData& d){ return d.timeStamp < minUsedTime; }),
            imuBuffer.end());
    }
    return std::make_pair(imuVec, uwbVec);
}

void Node::processThread()
{
    ros::Rate loopRate(uwbFreq);
    while (ros::ok()) {
        {
            std::unique_lock<std::mutex> imuLock(imuMutex); // unique_lock은 수동으로 락을 해제해야 함 condition_variable이랑 같이 쓰이고 락을 해제하는 것을 잊으면 데드락 걸림
            std::unique_lock<std::mutex> uwbLock(uwbMutex);

            dataCond.wait_for(uwbLock, std::chrono::milliseconds(20), [&]{ //0.02s, 50Hz 대기
                return (imuBuffer.size() >= 2 && uwbBufferRight.size() >= 1);
            });// IMU 보간하려면 UWB 1개당 IMU 2개 필요, 이 값이 들어올 때까지 또는 최대 5ms (200hZ) 대기 
        }

        std::pair<std::vector<ImuData>, std::vector<UwbData>> syncPair;
        while (true)//dataCond에서 UWB 1개, IMU 2개를 기다리거나 혹은 5ms 를 기다리거나 해서 깨어나긴 하지만 다른 스레드에서 이를 소비해 버리거나(erase 혹은 remove하면)했을 수도 있으니까 main thread로 들어가기전에 한번 더 확인
        {
            {
                std::lock_guard<std::mutex> imuLock(imuMutex);
                std::lock_guard<std::mutex> uwbLock(uwbMutex);
                if (imuBuffer.size() < 2 || uwbBufferRight.empty()) {
                    continue;
                }
            }
            
            syncPair = interpolationAllT_blocking();
            std::cout << "sncpair size: " << syncPair.first.size() << " " << syncPair.second.size() << std::endl;
            if (syncPair.first.size() == MppiEstimation.T && syncPair.second.size() == MppiEstimation.T) {
                break;
            }
            loopRate.sleep();
        }

        double curStamp = syncPair.second.front().timeStamp;
        double dt = curStamp - lastStamp;
        if (lastStamp > 0 && dt > 0) { 
            MppiEstimation.setDt(dt);
        }
        lastStamp = curStamp;

        MppiEstimation.solve(syncPair.second, syncPair.first);
        MppiEstimation.move(syncPair.first[MppiEstimation.T - 1]);

        {
            std::lock_guard<std::mutex> imuLock(imuMutex);
            std::lock_guard<std::mutex> uwbLock(uwbMutex);
        }
        
        loopRate.sleep();
    }
}

Eigen::Vector3d mppiEstimation::getTagPosition(const STATE& state, double y_offset) {
    return state.p + state.R * Eigen::Vector3d(0, y_offset, 0);
}

void Node::interpolateImuData(const ImuData &a, const ImuData &b, double t, ImuData &out) {
    double w = (t - a.timeStamp) / (b.timeStamp - a.timeStamp);
    out.timeStamp = t;
    out.gyr = a.gyr + w * (b.gyr - a.gyr);
    out.acc = a.acc + w * (b.acc - a.acc);
}
