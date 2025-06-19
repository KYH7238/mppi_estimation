#include "offlineMppiEstimation_2tag.h"
#include <Eigen/Dense>
#include <iostream>
#include <omp.h>

STATE::STATE() {
    p.setZero();
    R.setIdentity();
    v.setZero();
}

mppiEstimation::mppiEstimation()
: T(5), dimU(6), N(20000), dt(0), gammaU(10.0), TOL(1e-9), _g(0,0,9.81), urng(std::random_device{}()), normGen(0.0, 1.0)
{
    sigmaUvector = Eigen::VectorXd::Constant(dimU, 0.0);
    sigmaUvector << 10,10,10,1,1,1;
    sigmaU = sigmaUvector.asDiagonal();
    U0 = Eigen::MatrixXd::Zero(dimU, T);
    Uo = Eigen::MatrixXd::Zero(dimU, T);
    Xo.resize(T+1);
}

mppiEstimation::~mppiEstimation() {}

void mppiEstimation::setAnchorPositions(const Eigen::Matrix<double,3,8>& positions) {
    anchorPositions = positions;
}

void mppiEstimation::setDt(double deltaT) {
    dt = deltaT;
}

STATE mppiEstimation::f(const STATE& state, const ImuData& imu, const Eigen::VectorXd& u) {
    STATE next = state;
    Eigen::Matrix3d Rot = state.R;
    Eigen::Vector3d accWorld = Rot * (imu.acc - u.head(3)) + _g;
    next.p += state.v * dt + 0.5 * accWorld * dt * dt;
    next.v += accWorld * dt;
    next.R = Rot * Exp((imu.gyr - u.tail(3)) * dt);
    return next;
}

Eigen::MatrixXd mppiEstimation::getNoise(int T) {
    return sigmaU * normGen.template generate<Eigen::MatrixXd>(dimU, T, urng);
}

void mppiEstimation::move(const ImuData& imu) {
    xInit = f(xInit, imu, u0);
    U0.leftCols(T-1) = Uo.rightCols(T-1);
}

void mppiEstimation::solve(const std::vector<UwbData>& uwbData,
                           const std::vector<ImuData>& imuData)
{
    int M = imuData.size();
    Eigen::MatrixXd Ui = U0.replicate(N,1);
    Eigen::VectorXd costs(N);

    #pragma omp parallel for
    for(int i=0;i<N;++i) {
        std::vector<STATE> Xi(T+1);
        Eigen::MatrixXd noise = getNoise(T);
        Ui.middleRows(i*dimU, dimU) += noise;
        Xi[0]=xInit;
        double cost=0;
        for(int j=0;j<T;++j){
            Xi[j+1] = f(Xi[j], imuData[j], Ui.block(i*dimU, j, dimU,1));
            Eigen::Vector3d tagL = getTagPosition(Xi[j],+0.13);
            Eigen::Vector3d tagR = getTagPosition(Xi[j],-0.13);
            Eigen::VectorXd HxL = (anchorPositions.colwise()-tagL).colwise().norm();
            Eigen::VectorXd HxR = (anchorPositions.colwise()-tagR).colwise().norm();
            Eigen::VectorXd Hx(16); Hx<<HxL,HxR;
            cost += (uwbData[j].ranges - Hx).norm();
        }
        costs(i)=cost;
    }
    double minC=costs.minCoeff();
    Eigen::VectorXd w = (-gammaU*(costs.array()-minC)).exp();
    w/=w.sum();
    Uo.setZero();
    for(int i=0;i<N;++i)
        Uo += Ui.middleRows(i*dimU,dimU)*w(i);
    u0 = Uo.col(0);
    Xo[0]=xInit;
    for(int j=0;j<T;++j)
        Xo[j+1]=f(Xo[j], imuData[j], Uo.col(j));
}

Eigen::Vector3d mppiEstimation::getTagPosition(const STATE& state, double y_offset) {
    return state.p + state.R * Eigen::Vector3d(0,y_offset,0);
}

Eigen::MatrixXd mppiEstimation::Exp(const Eigen::Vector3d& omega) {
    double angle=omega.norm();
    if(angle<TOL) return Eigen::Matrix3d::Identity();
    Eigen::Vector3d a=omega/angle;
    double c=std::cos(angle), s=std::sin(angle);
    return c*Eigen::Matrix3d::Identity()+(1-c)*a*a.transpose()+s*vectorToSkewSymmetric(a);
}

Eigen::Matrix3d mppiEstimation::vectorToSkewSymmetric(const Eigen::Vector3d& v) {
    Eigen::Matrix3d m;
    m<<0,-v.z(),v.y(), v.z(),0,-v.x(), -v.y(),v.x(),0;
    return m;
}