#pragma once
#include <iostream>
#include <cstring>
#include <ros/ros.h>
#include <EigenRand/EigenRand>
#include <ctime>
#include <vector>
#include <omp.h>
#include <nlink_parser/LinktrackTagframe0.h>
#include <geometry_msgs/PoseStamped.h>
#include <sensor_msgs/Imu.h>
#include <Eigen/Dense>
#include <chrono>
#include <queue>
#include <utility>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <fstream>

class STATE {
public:
    Eigen::Vector3d p;
    Eigen::Matrix3d R;
    Eigen::Vector3d v;
    Eigen::Vector3d a_b;
    Eigen::Vector3d w_b;
    STATE();
};

class ImuData {
public:
    double timeStamp;
    Eigen::Vector3d acc;
    Eigen::Vector3d gyr;
};

class UwbData {
public:
    double timeStamp;
    Eigen::VectorXd ranges;
};

class EKF {
private:
    static constexpr size_t MAX_BUFFER_SIZE = 1000;
    std::vector<UwbData> uwbBufferLeft;
    std::vector<UwbData> uwbBufferRight;
    std::vector<ImuData> imuBuffer;
    
    Eigen::Matrix<double, 15, 15> covP;
    Eigen::Matrix<double, 15, 15> covQ;
    Eigen::Matrix<double, 16, 16> covR;
    Eigen::Matrix<double, 16, 1> vecZ;
    Eigen::Matrix<double, 16, 1> vecH;
    Eigen::Matrix<double, 15, 15> jacobianMatF;
    Eigen::Matrix<double, 16, 15> jacobianMatH;
    
    Eigen::Vector3d _g;
    double lastStamp;
    bool imuInit;
    bool uwbInit;
    double beforeT;
    bool initialized;
    double uwbInitTime;
    double imuInitTime;
    double TOL;
    double dt;
    STATE State;
    Eigen::Vector3d prevAcc;
    Eigen::Vector3d prevGyr;

    ros::NodeHandle nh_;
    ros::Publisher resultPuber;
    ros::Subscriber imuSub;
    ros::Subscriber uwbSub;
    std::mutex uwbMutex;
    std::mutex imuMutex;
    std::condition_variable dataCond;
    Eigen::Matrix<double, 3, 8> anchorPositions;

    void processThread();
    void processImuData();
    bool processUwbData();
    void interpolateImuData(double targetTime);
    void motionModel(const ImuData &imuData);
    void motionModelJacobian(const ImuData &imuData);
    void prediction(const ImuData &imuData);
    void measurementModel();
    void measurementModelJacobian();
    STATE correction();
    Eigen::Vector3d getTagPosition(const STATE &state, double y_offset);
    void publishPose(const STATE &state);
    Eigen::MatrixXd Exp(const Eigen::Vector3d &omega);
    Eigen::MatrixXd vectorToSkewSymmetric(const Eigen::Vector3d &vector);
    void uwbCallback(const nlink_parser::LinktrackTagframe0 &msg);
    void imuCallback(const sensor_msgs::ImuConstPtr &msg);

public:
    EKF();
    ~EKF() {}
    void setImuVar(const double stdV, const double stdW);
    void setDt(const double delta_t);
};