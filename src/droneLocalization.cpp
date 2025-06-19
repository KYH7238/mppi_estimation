#include "droneLocalization.h"
STATE::STATE() {
    p << 0.0, 0.0, 0.25;
    R.setIdentity();
    v.setZero();
    a_b.setZero();
    w_b.setZero();
}

EKF::EKF() {
    covP = Eigen::Matrix<double,15,15>::Identity() * 0.1;
    covQ.setIdentity()*0.001;
    double stdV, stdW;

    stdV = 0.0015387262937311438;
    stdW = 1.0966208586777776e-06;
    covQ.block<3,3>(9,9) = stdV*Eigen::Matrix3d::Identity();
    covQ.block<3,3>(12,12) = stdW*Eigen::Matrix3d::Identity();

    covR.setIdentity()*0.08;
    vecZ.setZero();
    vecH.setZero();
    jacobianMatF.setIdentity();
    _g = Eigen::Vector3d(0, 0, 9.81);
    lastStamp = 0;
    imuInit = false;    
    uwbInit = false;
    beforeT = 0;
    initialized = false;    
    uwbInitTime = 0;
    imuInitTime = 0;
    TOL = 1e-9;
    dt = 0;
    
    resultPuber = nh_.advertise<geometry_msgs::PoseStamped>("ekf_node", 1);
    imuSub = nh_.subscribe("/vectornav_driver_node/imu/data", 10, &EKF::imuCallback, this);
    uwbSub = nh_.subscribe("/nlink_linktrack_tagframe0", 10, &EKF::uwbCallback, this);
    anchorPositions <<  -4.4, -4.4,  4.43,  4.43, -4.4, -4.4, 4.43,  4.43,    
                        -4.04, 4.76, 4.89, -4.08, -4.08, 4.8, 4.93, -4.08,
                         0.26, 0.26, 0.3,   0.34,  2.45, 2.5, 2.47,  2.45;    
    std::thread(&EKF::processThread, this).detach();
}

void EKF::setImuVar(const double stdV, const double stdW) {
    covQ.block<3,3>(9,9) = stdV*Eigen::Matrix3d::Identity();
    covQ.block<3,3>(12,12) = stdW*Eigen::Matrix3d::Identity();
}

void EKF::setDt(const double delta_t) {
    dt = delta_t;
}

Eigen::MatrixXd EKF::Exp(const Eigen::Vector3d &omega) {
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

Eigen::MatrixXd EKF::vectorToSkewSymmetric(const Eigen::Vector3d &vector) {
    Eigen::Matrix3d Rot;
    Rot << 0, -vector.z(), vector.y(),
          vector.z(), 0, -vector.x(),
          -vector.y(), vector.x(), 0;
    return Rot;
}

void EKF::motionModel(const ImuData &imuData) {
    Eigen::Matrix3d Rot = State.R;
    Eigen::Vector3d accWorld = Rot*(imuData.acc - State.a_b) + _g;
    State.p += State.v*dt+0.5*accWorld*dt*dt;
    State.v += accWorld*dt;
    State.R = Rot*Exp((imuData.gyr - State.w_b)*dt);    
}

void EKF::motionModelJacobian(const ImuData &imuData) {
    jacobianMatF.setIdentity();
    Eigen::Matrix3d Rot = State.R;
    jacobianMatF.block<3, 3>(0,3) = 0.5*dt*dt*Eigen::Matrix3d::Identity();
    jacobianMatF.block<3, 3>(0,6) = dt*Eigen::Matrix3d::Identity();
    jacobianMatF.block<3, 3>(0,9) = -0.5*dt*dt*Rot;
    jacobianMatF.block<3, 3>(3,3) = Exp((imuData.gyr - State.w_b)*dt);
    jacobianMatF.block<3, 3>(3,12) = -Eigen::Matrix3d::Identity()*dt;
    jacobianMatF.block<3, 3>(6,9) = -Rot*dt;
    jacobianMatF.block<3, 3>(6,3) = -Rot*vectorToSkewSymmetric(imuData.acc-State.a_b)*dt;    
}

void EKF::prediction(const ImuData &imuData) {
    motionModelJacobian(imuData);
    motionModel(imuData);
    covP = jacobianMatF*covP*jacobianMatF.transpose() + covQ;
}

void EKF::measurementModel() {
    Eigen::Vector3d tagL = getTagPosition(State, 0.13);  
    Eigen::Vector3d tagR = getTagPosition(State, -0.13);  

    for (int i = 0; i < 8; ++i) {
        vecH(i)     = (tagL - anchorPositions.col(i)).norm();
        vecH(8 + i) = (tagR - anchorPositions.col(i)).norm();
    }
}

void EKF::measurementModelJacobian() {
    jacobianMatH.setZero();

    Eigen::Vector3d tagL = getTagPosition(State, 0.13);
    Eigen::Vector3d tagR = getTagPosition(State, -0.13);

    for (int i = 0; i < 8; ++i) {
        Eigen::Vector3d diffL = tagL - anchorPositions.col(i);
        double distL = vecH(i);
        if (distL > 1e-6) {
            jacobianMatH(i, 0) = diffL.x() / distL;
            jacobianMatH(i, 1) = diffL.y() / distL;
            jacobianMatH(i, 2) = diffL.z() / distL;
        }

        Eigen::Vector3d diffR = tagR - anchorPositions.col(i);
        double distR = vecH(8 + i);
        if (distR > 1e-6) {
            jacobianMatH(8 + i, 0) = diffR.x() / distR;
            jacobianMatH(8 + i, 1) = diffR.y() / distR;
            jacobianMatH(8 + i, 2) = diffR.z() / distR;
        }
    }
}

STATE EKF::correction() {
    measurementModel();
    measurementModelJacobian();
    Eigen::Matrix<double, 16, 16> residualCov;
    Eigen::Matrix<double, 15, 16> K;
    Eigen::Matrix<double, 15, 1> updateState;
    residualCov = jacobianMatH*covP*jacobianMatH.transpose() + covR;
    if (residualCov.determinant() == 0 || !residualCov.allFinite()) {
        std::cerr << "residualCov is singular or contains NaN/Inf" << std::endl;
        return State;
    }

    K = covP*jacobianMatH.transpose()*residualCov.inverse();
    updateState = K*(vecZ-vecH);
    State.p += updateState.segment<3>(0);
    State.R = State.R*Exp(updateState.segment<3>(3));
    State.v += updateState.segment<3>(6);
    State.a_b += updateState.segment<3>(9);
    State.w_b += updateState.segment<3>(12);
    covP = (Eigen::Matrix<double, 15, 15>::Identity()-(K*jacobianMatH))*covP;
    return State;    
}

void EKF::processThread() {
    std::unique_lock<std::mutex> lockU(uwbMutex, std::defer_lock);
    std::unique_lock<std::mutex> lockI(imuMutex, std::defer_lock);
    
    while (ros::ok()) {
        if (!dataCond.wait_for(lockU, std::chrono::milliseconds(100), [&]() {
            return !uwbBufferLeft.empty() && !uwbBufferRight.empty() && !imuBuffer.empty();
        })) {
            continue;
        }
        
        lockI.lock();
        
        processImuData();
        
        if (!processUwbData()) {
            lockI.unlock();
            continue;
        }
        
        lockI.unlock();
    }
}

void EKF::processImuData() {
    double tL = uwbBufferLeft.front().timeStamp;
    double tR = uwbBufferRight.front().timeStamp;
    double tU = std::max(tL, tR);
    
    while (!imuBuffer.empty() && imuBuffer.front().timeStamp <= tU) {
        ImuData imu = imuBuffer.front();
        imuBuffer.erase(imuBuffer.begin());
        
        double dt = imu.timeStamp - beforeT;
        setDt(dt);

        prediction(imu);
        beforeT = imu.timeStamp;
    }
    
    if (!imuBuffer.empty() && beforeT < tU) {
        interpolateImuData(tU);
    }
}

bool EKF::processUwbData() {
    if (uwbBufferLeft.empty() || uwbBufferRight.empty()) {
        return false;
    }
    
    UwbData dL = uwbBufferLeft.front();
    UwbData dR = uwbBufferRight.front();
    
    uwbBufferLeft.erase(uwbBufferLeft.begin());
    uwbBufferRight.erase(uwbBufferRight.begin());
    
    vecZ.segment<8>(0) = dL.ranges;
    vecZ.segment<8>(8) = dR.ranges;
    
    State = correction();
    publishPose(State);
    
    return true;
}

void EKF::interpolateImuData(double targetTime) {
    ImuData next = imuBuffer.front();
    double alpha = (targetTime - beforeT) / (next.timeStamp - beforeT);
    
    ImuData interp;
    interp.timeStamp = targetTime;
    interp.acc = (1-alpha)*prevAcc + alpha*next.acc;
    interp.gyr = (1-alpha)*prevGyr + alpha*next.gyr;
    
    setDt(targetTime - beforeT);
    prediction(interp);
    
    prevAcc = interp.acc;
    prevGyr = interp.gyr;
    beforeT = targetTime;
}

void EKF::uwbCallback(const nlink_parser::LinktrackTagframe0 &msg) {
    std::lock_guard<std::mutex> lock(uwbMutex);
    
    double t = msg.system_time/1000.0;
    if (!uwbInit) {
        uwbInitTime = t;
        uwbInit = true;
    }
    
    UwbData d;
    d.timeStamp = t - uwbInitTime;
    d.ranges = Eigen::VectorXd(8);
    for (int i = 0; i < 8; i++) {
        d.ranges(i) = msg.dis_arr[i];
    }
    
    if (msg.id == 0) {
        uwbBufferLeft.push_back(d);
        if (uwbBufferLeft.size() > MAX_BUFFER_SIZE) {
            uwbBufferLeft.erase(uwbBufferLeft.begin());
        }
    } else if (msg.id == 1) {
        uwbBufferRight.push_back(d);
        if (uwbBufferRight.size() > MAX_BUFFER_SIZE) {
            uwbBufferRight.erase(uwbBufferRight.begin());
        }
    }
    
    dataCond.notify_one();
}

void EKF::imuCallback(const sensor_msgs::ImuConstPtr &msg) {
    std::lock_guard<std::mutex> lock(imuMutex);
    double t = msg->header.stamp.toSec();

    if (!imuInit) {
        imuInit     = true;
        imuInitTime = t;
        beforeT     = 0.0;  
        prevAcc     = Eigen::Vector3d::Zero();
        prevGyr     = Eigen::Vector3d::Zero();
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
    if (imuBuffer.size() > MAX_BUFFER_SIZE)
        imuBuffer.erase(imuBuffer.begin());

    dataCond.notify_one();
}

Eigen::Vector3d EKF::getTagPosition(const STATE &state, double y_offset) {
    return state.p + state.R * Eigen::Vector3d(0, y_offset, 0);
}

void EKF::publishPose(const STATE &state) {
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