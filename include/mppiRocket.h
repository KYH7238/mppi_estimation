#include <ros/ros.h>
#include <iostream>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
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
    // MPPI(ros::NodeHandle& nh);
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
    void h(Eigen::Ref<Eigen::MatrixXd> U_seq);
    void saveCost(double cost);
    void plotCost();
    void publishTrajectory(const Eigen::MatrixXd &trajectory);
    int T;
    double dt;
protected:
    Eigen::Rand::NormalGen<double> norm_gen{0, 1};  
    Eigen::Vector2d g_; 
    Eigen::MatrixXd sigma_u;
    std::mt19937_64 urng{static_cast<std::uint_fast64_t>(std::time(nullptr))};    
    int N, dim_x, dim_u, dim_g, dim_h;       
    double gamma_u, l, mass, I, F_min, F_max, C;       
    ros::Publisher traj_pub;
};