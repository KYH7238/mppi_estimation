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

class STATE {
public:
    Eigen::Vector3d p;
    Eigen::Matrix3d R;
    Eigen::Vector3d v;
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

class mppiEstimation {
public:
    mppiEstimation();
    ~mppiEstimation();
    void setAnchorPositions(const Eigen::Matrix<double, 3, 8> &positions);
    void setDt(const double deltaT);
    STATE f(const STATE &state, const ImuData &imuData, const Eigen::VectorXd Ui);
    Eigen::MatrixXd getNoise(const int &T);
    void move(const ImuData &imuData);
    void solve(const std::vector<UwbData> &uwbData, const std::vector<ImuData> &imuData);
    void publishPose(const STATE &state);
    Eigen::MatrixXd Exp(const Eigen::Vector3d &omega);
    Eigen::MatrixXd vectorToSkewSymmetric(const Eigen::Vector3d &vector);
    Eigen::Matrix<double, 3, 8> anchorPositions;
    Eigen::MatrixXd U0;
    STATE xInit;
    Eigen::MatrixXd Uo;
    std::vector<STATE> Xo; 
    Eigen::VectorXd u0;
    Eigen::Vector3d _g;
    int T;
    float dt;
    int N;
    int dimU;
    double gammaU;    
    double TOL;
    Eigen::MatrixXd sigmaU;
    Eigen::VectorXd sigmaUvector;
    ros::NodeHandle nh;
    ros::Publisher resultPuber;
    ros::Subscriber uwbSub;
    ros::Subscriber imuSub;
    std::mt19937_64 urng{static_cast<std::uint_fast64_t>(std::time(nullptr))};    
    Eigen::Rand::NormalGen<double> normGen{0.0, 1.0};    
};

class Node {
public:
    ros::NodeHandle nh_;
    ros::Subscriber uwbSub;
    ros::Subscriber imuSub;
    mppiEstimation MppiEstimation;
    std::queue<UwbData> uwbDataQueue;
    std::queue<ImuData> imuDataQueue;

    std::vector<sensor_msgs::ImuConstPtr> imuBuffer;
    std::vector<nlink_parser::LinktrackTagframe0::ConstPtr> uwbBuffer;

    std::vector<sensor_msgs::ImuConstPtr> syncImu;
    std::vector<nlink_parser::LinktrackTagframe0::ConstPtr> syncUwb;

    std::mutex imuMutex, uwbMutex; 
    std::condition_variable dataCond;
    ImuData lastUwbImu;
    ImuData lastImu;    
    UwbData uwbData;
    ImuData imuData;

    bool imuInit;    
    bool uwbInit;
    bool initialized;
    double uwbInitTime;
    double imuInitTime, lastStamp;
    float dt;
    double beforeT;
    int cnt;
    int uwbFreq;

    Eigen::Matrix<double, 3, 8> anchorPositions;
    Node();
    void uwbCallback(const nlink_parser::LinktrackTagframe0 &msg);
    void imuCallback(const sensor_msgs::ImuConstPtr &msg);

    ImuData fromImuMsg (const sensor_msgs::Imu & msg);
    UwbData fromUwbMsg (const nlink_parser::LinktrackTagframe0 &msg); 
    std::pair<std::vector<ImuData>, std::vector<UwbData>> interpolationAllT_blocking();
    void interpolateImuData(const sensor_msgs::ImuConstPtr &firstData, const sensor_msgs::ImuConstPtr &secondData, double curStamp, sensor_msgs::Imu & interData);
    void processThread();
};


