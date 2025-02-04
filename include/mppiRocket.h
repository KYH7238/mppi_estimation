#include <iostream>
#include <Eigen/Dense>
#include <EigenRand/EigenRand>
#include <chrono>
#include <cmath>
#include <vector>
#include <functional>
#include <random>
#include <omp.h>
#include "matplotlibcpp.h"

class MPPI {
public:    
    MPPI();
    ~MPPI(); 
    Eigen::MatrixXd U_0, Uo, Xo;   
    Eigen::VectorXd u0, x_init, x_target;
    std::vector<double> iteration_costs;        
    Eigen::MatrixXd getNoise(const int &T);
    void solve();
    void move();
    Eigen::VectorXd f(const Eigen::VectorXd& x, const Eigen::VectorXd& u);
    double q(const Eigen::VectorXd& x, const Eigen::VectorXd& u);
    double p(const Eigen::VectorXd& x, const Eigen::VectorXd& x_target);
    void h(Eigen::MatrixXd& U_seq);
    void saveCost(double cost);
    void plotCost();

protected:
    Eigen::Rand::NormalGen<double> norm_gen{0.0, 1.0};  
    Eigen::Vector2d g_; 
    Eigen::MatrixXd sigma_u;
    std::mt19937_64 urng{static_cast<std::uint_fast64_t>(std::time(nullptr))};    
    int N, dim_x, dim_u, dim_g, dim_h, T;       
    double umax, gamma_u, l, dt, mass, I;       
    

};