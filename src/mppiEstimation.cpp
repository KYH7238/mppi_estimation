#include <iostream>
#include <cstring>
#include <ros/ros.h>
#include <EigenRand/EigenRand>
#include <ctime>
#include <vector>
#include <chrono>
#include <iostream>
#include <omp.h>
#include "mppiEstimation.h"

mppiSmoothing::mppiSmoothing() {
    _g << 0, 0, 9.81;
    TOL = 1e-9;
    dt = 0;

    STATE.p << 0, 0, 0;
    STATE.R.setIdentity();
    STATE.v.setZero();

    resultPuber = nh.advertise<geometry_msgs::PoseStamped>("mppi_pose", 1);
}

void mppiSmoothing::interpolation() {

}

void mppiSmoothing::f(const ImuData<double> &imuData) {
    Eigen::Matrix3d Rot = STATE.R;
    // Eigen::Vector3d accWorld = Rot*(imuData.acc - STATE.a_b) + _g;
    Eigen::Vector3d accWorld = Rot*imuData.acc + _g;
    STATE.p += STATE.v*dt+0.5*accWorld*dt*dt;
    STATE.v += accWorld*dt;
    // STATE.R = Rot*Exp((imuData.gyr - STATE.w_b)*dt);
    STATE.R = Rot*Exp(imuData.gyr*dt);
}

Eigen::MatrixXd mppiSmoothing::getNoise(const int &T) {
    return sigmaU * normGen.template generate<Eigen::MatrixXd>(dimU, T, urng);
}

void mppiSmoothing::move() {
    xInit = xInit + (dt * f(xInit, u0));
    U0.leftcols(T-1) = Uo.rightCols(T-1); 
}

void mppiSmoothing::solve(const Eigen::Matrix<double,3,8> &anchors, const Eigen::MatrixXd &yMeas) {
    start = std::chrono::high_resolution_clock::now();

    Eigen::MatrixXd Ui = U_0.replicate(N, 1);
    Eigen::VectorXd costs(N);
    Eigen::VectorXd weights(N);

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        Eigen::MatrixXd Xi(dimX, T+1);
        Eigen::MatrixXd noise = getNoise(T);
        Ui.middleRows(i * dimU, dimU) += noise;

        Xi.col(0) = xInit;
        double cost = 0.0;

        for (int j = 0; j < T; ++j) {
            Xi.col(j+1) = Xi.col(j) + (dt * f(Xi.col(j), Ui.block(i * dimU, j, dimU, 1)));

            Eigen::VectorXd Hx(8);
            for (int i = 0; i < 8; ++i) {
                double dx = Xi(0,j) - anchors(0,i);
                double dy = Xi(1,j) - anchors(1,i);
                double dz = Xi(2,j) - anchors(2,i);
                Hx(i) = std::sqrt(dx*dx + dy*dy + dz*dz);
            }

            double stepCost = (yMeas.col(j) - Hx).squaredNorm();
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

    Xo.col(0) = xInit;
    for (int j = 0; j < T; ++j) {
        Xo.col(j+1) = Xo.col(j) + (dt * f(Xo.col(j), Uo.col(j)));
    }

    visual_traj.push_back(xInit);

    publishPose(xInit);
    publishPose(Xo);

}

void mppiSmoothing::publishPose(const Eigen::VectorXd &state) {
    geometry_msgs::PoseStamped pose;
    pose.header.frame_id = "map";
    pose.header.stamp = ros::Time::now();

    pose.pose.position.x = state(0);
    pose.pose.position.y = state(1);
    pose.pose.position.z = state(2);

    pose.pose.orientation.w = 1.0;
    pose.pose.orientation.x = 0.0;
    pose.pose.orientation.y = 0.0;
    pose.pose.orientation.z = 0.0;

    resultPuber.publish(pose);
}

Eigen::MatrixXd mppiSmoothing::exp(const Eigen::Vector3d &omega) {
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

Eigen::MatrixXd mppiSmoothing::vectorToSkewSymmetric(const Eigen::Vector3d &vector) {
    Eigen::Matrix3d Rot;
    Rot << 0, -vector.z(), vector.y(),
          vector.z(), 0, -vector.x(),
          -vector.y(), vector.x(), 0;
    
    return Rot;
}

mppiSmoothing::~mppiSmoothing() { 

}