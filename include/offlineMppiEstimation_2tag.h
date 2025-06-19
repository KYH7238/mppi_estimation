#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <EigenRand/EigenRand>
#include <random>
#include <vector>

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
    int id;                  // 0: left tag, 1: right tag (unused offline merged)
    double timeStamp;
    Eigen::VectorXd ranges;  // size = 8 for each tag (merged into 16 later)
};

class mppiEstimation {
public:
    mppiEstimation();
    ~mppiEstimation();

    void setAnchorPositions(const Eigen::Matrix<double,3,8>& positions);
    void setDt(double deltaT);
    STATE f(const STATE& state, const ImuData& imu, const Eigen::VectorXd& u);
    Eigen::MatrixXd getNoise(int T);
    void move(const ImuData& imu);
    void solve(const std::vector<UwbData>& uwbData,
               const std::vector<ImuData>& imuData);
    static Eigen::Vector3d getTagPosition(const STATE& state, double y_offset);

    Eigen::Matrix<double,3,8> anchorPositions;
    int T, dimU, N;
    double dt, gammaU, TOL;
    Eigen::MatrixXd sigmaU, U0, Uo;
    Eigen::VectorXd sigmaUvector, u0;
    std::vector<STATE> Xo;
    STATE xInit;
    Eigen::Vector3d _g;

private:
    Eigen::MatrixXd Exp(const Eigen::Vector3d& omega);
    Eigen::Matrix3d vectorToSkewSymmetric(const Eigen::Vector3d& v);
    Eigen::Rand::NormalGen<double> normGen;
    std::mt19937_64 urng;
};