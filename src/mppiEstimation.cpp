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
    std::thread process(processThread);
    ros::MultiThreadedSpinner spinner(2);
    spinner.spin();
    process.join();

    imuInit = false;    
    uwbInit = false;
    initialized = false;
    lastStamp = 0;
    uwbInitTime = 0;
    imuInitTime = 0;
    dt = 0;
    beforeT = 0;
    cnt = 0;
    uwbFreq = 50;
}

ImuData Node::fromImuMsg (const sensor_msgs::Imu & msg) {
    ImuData data;
    data.stamp = msg.header.stamp.toSec();
    data.gyr[0] = msg.angular_velocity.x;
    data.gyr[1] = msg.angular_velocity.y;
    data.gyr[2] = msg.angular_velocity.z;
    data.acc[0] = msg.linear_acceleration.x;
    data.acc[1] = msg.linear_acceleration.y;
    data.acc[2] = msg.linear_acceleration.z;

    return data;
}

UwbData Node::fromUwbMsg (const nlink_parser::LinktrackTagframe0 &msg) {
    UwbData data;
    data.timeStamp = msg.system_time/1000.00  - uwbInitTime;  
    for (int i = 0; i < 8; i++) {
        data.ranges[i] = msg.dis_arr[i];
    }    
    return data;
}

void Node::interpolateImuData(const sensor_msgs::ImuConstPtr &firstData, const sensor_msgs::ImuConstPtr &secondData, double curStamp, sensor_msgs::Imu &interData) {
    double firstStamp = firstData->header.stamp.toSec();
    double secondStamp = secondData->header.stamp.toSec();
    double scale = (curStamp - firstStamp) / (secondStamp - firstStamp);
    interData = *firstData;
    interData.angular_velocity.x = scale * (secondStamp->angular_velocity.x - firstStamp->angular_velocity.x) + firstStamp->angular_velocity.x;
    interData.angular_velocity.y = scale * (secondStamp->angular_velocity.y - firstStamp->angular_velocity.y) + firstStamp->angular_velocity.y;
    interData.angular_velocity.z = scale * (secondStamp->angular_velocity.z - firstStamp->angular_velocity.z) + firstStamp->angular_velocity.z;
    interData.linear_acceleration.x = scale * (secondStamp->linear_acceleration.x - firstStamp->linear_acceleration.x) + firstStamp->linear_acceleration.x;
    interData.linear_acceleration.y = scale * (secondStamp->linear_acceleration.y - firstStamp->linear_acceleration.y) + firstStamp->linear_acceleration.y;
    interData.linear_acceleration.z = scale * (secondStamp->linear_acceleration.z - firstStamp->linear_acceleration.z) + firstStamp->linear_acceleration.z;
}

// void Node::processThread()
// {
//     ros::Rate loopRate(uwbFreq);
//     while(ros::ok()) {
//         std::unique_lock<mutex> imuLock(imuMutex);
//         std::unique_lock<mutex> uwbLock(uwbMutex);

//         if (!imuBuffer.size() && !uwbBuffer.size()) {
//             ROS_INFO_THROTTLE(10, "wait for uwb or imu msg ......");
//             imuLock.unlock();
//             uwbLock.unlock();
//             loopRate.sleep();
//             continue;
//         }

//         if (!initialized) {

//             if (!imuBuffer.size() || !uwbBuffer.size()) {
//                 ROS_INFO_THROTTLE(10, "wait for uwb or imu msg ......");
//                 imuLock.unlock();
//                 uwbLock.unlock();
//                 loopRate.sleep();
//                 continue;
//             }
//             Eigen::Vector<UwbData> UwbDatas;
//             for (auto &uwbMsg : uwbBuffer) {
//                 UwbDatas.push_back(fromUwbMsg(*uwbMsg));
//             }

//             auto iter = uwbBuffer.begin();
//             for (; iter != uwbBuffer.end(); iter++) {
//                 if ((*iter)->timeStamp > imuBuffer[0]->header.stamp.toSec())
//                     break;
//             }

//             if (uwbBuffer.begin() == iter || uwbBuffer.end() == iter) {
//                 if (uwbBuffer.begin() == iter) {
//                     imuBuffer.erase(imuBuffer.begin());
//                 }
//                 imuLock.unlock();
//                 uwbLock.unlock();
//                 loopRate.sleep();
//                 continue;
//             }
            
//             sensor_msgs::Imu interImu;
//             double curStamp = uwbBuffer[0] -> timeStamp;
//             interpolateImuData(imuBuffer[0], imuBuffer[1], curStamp, interImu);
//             lastUwbImu = fromImuMsg(interImu); //ImuData
//             lastUwbImu.timeStamp = curStamp;
//             lastImu = lastUwbImu;

//             syncImu.push_back(lastImu);
//             syncUwb.push_back(uwbBuffer[0]);

//             imuBuffer.erase(imuBuffer.begin());
//             uwbBuffer.erase(uwbBuffer.begin(), iter);

//             imuLock.unlock();
//             uwbLock.unlock();
//             loopRate.sleep();

//             initialized = true;
//             continue;            
//         }

//         if (uwbBuffer.size() !=0) {

//             if (!imuBuffer.size() || !uwbBuffer.size()) {
//                 ROS_INFO_THROTTLE(10, "wait for uwb or imu msg ......");
//                 imuLock.unlock();
//                 uwbLock.unlock();
//                 loopRate.sleep();
//                 continue;
//             }
//             Eigen::Vector<UwbData> UwbDatas;
//             for (auto &uwbMsg : uwbBuffer) {
//                 UwbDatas.push_back(fromUwbMsg(*uwbMsg));
//             }
//             auto iter = uwbBuffer.begin();
//             for (; iter != uwbBuffer.end(); iter++) {
//                 if ((*iter)->timeStamp > imuBuffer[0]->header.stamp.toSec())
//                     break;
//             }
//             if (uwbBuffer.begin() == iter || uwbBuffer.end() == iter) {
//                 if (uwbBuffer.begin() == iter) {
//                     imuBuffer.erase(imuBuffer.begin());
//                 }
//                 imuLock.unlock();
//                 uwbLock.unlock();
//                 loopRate.sleep();
//                 continue;
//             }

//             sensor_msgs::Imu interImu;
//             double curStamp = uwbBuffer[0] -> timeStamp;
//             interpolateImuData(imuBuffer[0], imuBuffer[1], curStamp, interImu);
//             lastUwbImu = fromImuMsg(interImu); //ImuData
//             lastUwbImu.timeStamp = curStamp;
//             lastImu = lastUwbImu;

//             syncImu.push_back(lastImu);
//             syncUwb.push_back(uwbBuffer[0]);

//             imuBuffer.erase(imuBuffer.begin());
//             uwbBuffer.erase(uwbBuffer.begin(), iter);            

//             imuLock.unlock();
//             uwbLock.unlock();
//             loopRate.sleep();

//             if (syncImu.size()== T && syncUwb.size() == T) {
//                 dt = syncUwb(0)->timeStamp - syncUwb(1)->timeStamp 
//                 MppiEstimation.setDt(dt);
//                 MppiEstimation.solve(syncUwb,syncImu);
//                 MppiEstimation.move(syncImu[9]);
//                 syncUwb.erase(syncUwb.begin());
//                 syncImu.erase(syncImu.begin());
//             }         

//             lastUwbImu.timeStamp = cur_stamp;
//             lastImu = lastUwbImu;
//             uwbBuffer.erase(uwbBuffer.begin());
//             imuBuffer.erase(imuBuffQer.begin(), iter);
//         }

//         imu_lock.unlock();
//         uwb_lock.unlock();
//         loop_rate.sleep();
//         continue;        
//     }
// }

void Node::processThread()
{
    ros::Rate loopRate(uwbFreq);

    while (ros::ok()) {
        {
            std::unique_lock<std::mutex> imuLock(imuMutex);
            std::unique_lock<std::mutex> uwbLock(uwbMutex);

            dataCond.wait_for(uwbLock, std::chrono::milliseconds(20), [&]{ //0.02s, 50Hz
                return (imuBuffer.size() >= 2 && uwbBuffer.size() >= 1);
            });
        }

        std::pair<std::vector<ImuData>, std::vector<UwbData>> syncPair;
        while (true)
        {
            {
                std::lock_guard<std::mutex> imuLock(imuMutex);
                std::lock_guard<std::mutex> uwbLock(uwbMutex);
                if (imuBuffer.size() < 2 || uwbBuffer.empty()) {
                    continue;
                }
            }
            
            syncPair = interpolationAllT_blocking();
            if (syncPair.first.size() == MppiEstimation.T && syncPair.second.size() == MppiEstimation.T) break;
            loopRate.sleep();
        }

        {
            MppiEstimation.solve(syncPair.second, syncPair.first);
            MppiEstimation.move(syncPair.first[MppiEstimation.T - 1]);
        }

        {
            std::lock_guard<std::mutex> imuLock(imuMutex);
            std::lock_guard<std::mutex> uwbLock(uwbMutex);
        }
        
        loopRate.sleep();
    }
}

std::pair<std::vector<ImuData>, std::vector<UwbData>> Node::interpolationAllT_blocking()
{
    int tCount = MppiEstimation.T;
    std::vector<ImuData> imuVec;
    std::vector<UwbData> uwbVec;
    imuVec.reserve(tCount);
    uwbVec.reserve(tCount);

    while (imuVec.size() < tCount) {
        {
            std::lock_guard<std::mutex> imuLock(imuMutex);
            std::lock_guard<std::mutex> uwbLock(uwbMutex);
            if (imuBuffer.size() < 2 || uwbBuffer.empty()) {
                continue;
            }
        }
        loopRate.sleep();

        std::unique_lock<std::mutex> imuLock(imuMutex);
        std::unique_lock<std::mutex> uwbLock(uwbMutex);

        while (!uwbBuffer.empty() &&
               (uwbBuffer.front()->header.stamp.toSec() < imuBuffer.front()->header.stamp.toSec() ||
                uwbBuffer.front()->header.stamp.toSec() > imuBuffer.back()->header.stamp.toSec()))
        {
            uwbBuffer.erase(uwbBuffer.begin());
        }

        if (imuBuffer.size() < 2 || uwbBuffer.empty()) {
            continue;  
        }

        sensor_msgs::ImuConstPtr firstImu = imuBuffer[0];
        sensor_msgs::ImuConstPtr secondImu = imuBuffer[1];
        double curStamp = uwbBuffer[0]->header.stamp.toSec();
        dt = curStamp - lastStamp;
        MppiEstimation.setDt(dt)
        sensor_msgs::Imu interImuMsg;
        interpolateImuData(firstImu, secondImu, curStamp, interImuMsg);

        ImuData interImu = fromImuMsg(interImuMsg);
        interImu.timeStamp = curStamp;

        UwbData curUwb = fromUwbMsg(*uwbBuffer[0]);

        imuVec.push_back(interImu);
        uwbVec.push_back(curUwb);
        lastStamp = curStamp;
    
        imuBuffer.erase(imuBuffer.begin());
        uwbBuffer.erase(uwbBuffer.begin());
    }

    return std::make_pair(imuVec, uwbVec);
}

void Node::uwbCallback(const nlink_parser::LinktrackTagframe0 &msg) {
    std::unique_lock<mutex> lock(uwbMutex);
    if(!uwbInit){
        uwbInitTime = msg.system_time/1000.00;
        uwbInit = true;
    }
    UwbData data;
    data.timeStamp = msg.system_time/1000.00  - uwbInitTime;  
    for (int i = 0; i < 8; i++) {
        data.ranges[i] = msg.dis_arr[i];
    }
    uwbBuffer.push_back(data);
}

void Node::imuCallback(const sensor_msgs::ImuConstPtr &msg) {   
    std::unique_lock<mutex> lock(imuMutex);
    if (!imuInit){
        imuInitTime = msg->header.stamp.toSec();
        imuInit = true;
    }
    else {
    sensor_msgs::ImuConstPtr data;
    data.timeStamp = msg->header.stamp.toSec()-imuInitTime;  
    imuBuffer.push_back(data); 
    }
}

