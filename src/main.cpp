#include <ros/ros.h>
#include <nlink_parser/LinktrackTagframe0.h>
#include <sensor_msgs/Imu.h>
#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include <deque> 
#include "mppiEstimation.cpp"

class MppiEstimation {
public:
    MppiEstimation() {

        anchorpositions << 0,    0,    8.86, 8.86, 0,    0,    8.86,  8.86,
                           0,    8.00, 8.00, 0,    0,    8.00,  8.00,  0,
                           0.0,  0.0,  0.0,  0.0,  2.20, 2.20,  2.20,  2.20;

        auto model = Quadrotor();

        param.dt = 0.02;
        param.T = 10;
        param.x_init.resize(model.dim_x);
        param.x_init << 4.36, 4.0, 0.0, 0.0, 0.0, 0.0;
        param.N = 6000;
        param.gamma_u = 10.0;
        Eigen::VectorXd sigma_u(model.dim_u);
        sigma_u << 1.0, 1.0, 10.8;
        // sigma_u << 0.05, 0.05, 9.85;
        param.sigma_u = sigma_u.asDiagonal();

        param.y_meas.resize(8, param.T); 
        param.y_meas.setZero(); 

        solver = new MPPI(model);
        solver->U_0 = Eigen::MatrixXd::Zero(model.dim_u, param.T);
        solver->U_0.row(2).array() += model.g; 
        solver->init(param);

        maxiter = 200;
    }

    ~MppiEstimation() {
        delete solver;
    }

    void pushMeasurement(const Eigen::VectorXd &uwb_data) {
        uwbDataQueue.push_back(uwb_data);

        if ((int)uwbDataQueue.size() > param.T) {
            uwbDataQueue.pop_front();
        }
    }

    void updateYMeasFromQueue() {

        int idx = 0;
        for (auto it = uwbDataQueue.begin(); it != uwbDataQueue.end(); ++it, ++idx) {

            param.y_meas.col(idx) = *it;
        }
    }

    void run() {
        if ((int)uwbDataQueue.size() == param.T) {
            updateYMeasFromQueue(); 
            for (int i = 0; i < maxiter; ++i) {

                solver->solveEstimation(anchorpositions, param.y_meas);           
                solver->move();
            }
        } else {
            ROS_WARN("Not enough UWB data: need %d steps, have %d steps", param.T, (int)uwbDataQueue.size());
        }
    }

    Quadrotor model;
    MPPIParam param;
    MPPI* solver;
    Eigen::Matrix<double,3,8> anchorpositions;
    int maxiter;

    std::deque<Eigen::VectorXd> uwbDataQueue; 
};


class Node {
private:
    ros::NodeHandle nh_;
    ros::Subscriber uwb_sub;
    // ros::Subscriber imu_sub;
    MppiEstimation mppiEstimation;
    Eigen::VectorXd uwbData;
    Eigen::VectorXd imuData;

public:
    Node(): uwbData(8) {
        uwb_sub = nh_.subscribe("/nlink_linktrack_tagframe0", 10, &Node::uwbCallback, this);
        // imu_sub = nh_.subscribe("/imu/data", 10, &Node::uwbCallback, this);
    }

    void uwbCallback(const nlink_parser::LinktrackTagframe0 &msg) {
        for (int i = 0; i < 8; i++) {
            uwbData(i) = msg.dis_arr[i];
        }

        mppiEstimation.pushMeasurement(uwbData);

        if ((int)mppiEstimation.uwbDataQueue.size() == mppiEstimation.param.T) {
            
            mppiEstimation.run();
        } else {
            // std::cout << (int)mppiEstimation.uwbDataQueue.size() <<std::endl;
        }
    }
};

int main(int argc, char** argv) {
    ros::init(argc, argv, "mppi_estimation");
    Node node;
    ros::spin();
    return 0;    
}