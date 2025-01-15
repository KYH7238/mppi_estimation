#pragma once

class mppiSmoothing {
public:
    mppiSmoothing();
    ~mppiSmoothing();
    void setAnchorPositions(const Eigen::Matrix<double, 3, 8> &anchorpositions);
    void interpolation();
    void f(const STATE state, const ImuData<double> &imuData, const Eigen::VectorXd Ui);
    Eigen::MatrixXd getNoise(const int &T);
    void move(const ImuData<double> &imuData);
    void solve(const Eigen::Matrix<double,3,8> &anchors, const Eigen::MatrixXd &yMeas, const ImuData<double> &imuData);
    void publishPose(const Eigen::VectorXd &state);
    Eigen::MatrixXd exp(const Eigen::Vector3d &omega);
    Eigen::MatrixXd vectorToSkewSymmetric(const Eigen::Vector3d &vector);

    Eigen::Matrix<double, 3, 8> anchorPositions;
    Eigen::MatrixXd U0;
    Eigen::MatrixXd xInit;
    Eigen::MatrixXd xTarget;
    Eigen::MatrixXd Uo;
    Eigen::MatrixXd Xo;
    Eigen::VectorXd u0;
    Eigen::MatrixXd processNoise;
    Eigen::VectorXd measurementNoise;
    Eigen::VectorXd _g;

private:
    float dt;
    int T;
    int N;
    int dimX;
    int dimU;
    double gammaU;    
    double TOL;
    Eigen::MatrixXd sigmaU;
    Eigen::Vector3d sigmaUvector;
    ros::NodeHandle nh;
    ros::Publisher resultPuber;

    std::mt19937_64 urng{static_cast<std::uint_fast64_t>(std::time(nullptr))};    
    Eigen::Rand::NormalGen<double> normGen{0.0, 1.0};    

};

class ImuData {
    double timeStamp;
    Eigen::Matrix3d acc;
    Eigen::Matrix3d gyr;
};

class STATE {
public:
    Eigen::Vector3d p;
    Eigen::Matrix3d R;
    Eigen::Vector3d v;

    STATE getState();
    STATE setState(const STATE &setState);
}