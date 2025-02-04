#include "mppiRocket.h"
namespace plt = matplotlibcpp;
MPPI::MPPI() {

    N = 1000;
    dim_x = 6;
    dim_u = 2;
    dim_g = 1;    
    dim_h = 2;  

    l = 0.7; 
    dt = 0.1;
    mass = 10.0;
    I = 3.33;
    g_ << -9.81, 0.0;
    umax = mass * 9.81 * 1.1;
    
    T = 100;
    gamma_u = 1.0;
    sigma_u = 64.745 * Eigen::MatrixXd::Identity(dim_u, dim_u); //F_max = 107.91, F_min = 21.58

    U_0 = 9.81 * 10.0 * Eigen::MatrixXd::Ones(dim_u, T);
    Xo = Eigen::MatrixXd::Zero(dim_x, T+1);

    x_init = Eigen::VectorXd::Zero(dim_x);
    x_init(0) = 10.0;
    x_init(1) = 8.0;

    x_target = Eigen::VectorXd::Zero(dim_x);
}

MPPI::~MPPI() { }

Eigen::MatrixXd MPPI::getNoise(const int &T) {
    return sigma_u * norm_gen.template generate<Eigen::MatrixXd>(dim_u, T, urng);
}

void MPPI::solve() {
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::MatrixXd Ui = Eigen::MatrixXd::Zero(N * dim_u, T);
    for (int i = 0; i < N; i++) {
        Ui.block(i * dim_u, 0, dim_u, T) = U_0;
    }

    Eigen::VectorXd costs = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        Eigen::MatrixXd noise = getNoise(T);
        Eigen::MatrixXd U_traj = Ui.block(i * dim_u, 0, dim_u, T);
        U_traj += noise;
        h(U_traj);

        Eigen::MatrixXd Xi = Eigen::MatrixXd::Zero(dim_x, T+1);
        Xi.col(0) = x_init;
        double cost = 0.0;
        for (int j = 0; j < T; j++) {
            cost += q(Xi.col(j), U_traj.col(j));
            Xi.col(j + 1) = f(Xi.col(j), U_traj.col(j));
        }
        cost += p(Xi.col(T), x_target);
        costs(i) = cost;

        Ui.block(i * dim_u, 0, dim_u, T) = U_traj;
    }

    double min_cost = costs.minCoeff();
    for (int i = 0; i < N; i++) {
        weights(i) = exp(-gamma_u * (costs(i) - min_cost));
    }
    double total_weight = weights.sum();
    weights /= total_weight;

    Uo = Eigen::MatrixXd::Zero(dim_u, T);
    for (int i = 0; i < N; i++) {
        Uo += Ui.block(i * dim_u, 0, dim_u, T) * weights(i);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "MPPI solve time: " << elapsed.count() << " seconds" << std::endl;

    u0 = Uo.col(0);
    Xo.col(0) = x_init;
    for (int j = 0; j < T; j++) {
        Xo.col(j + 1) = f(Xo.col(j), Uo.col(j));
    }

    saveCost(min_cost);
}

void MPPI::move() {
    x_init = f(x_init, u0);
    if (T > 1) {
        U_0.leftCols(T - 1) = Uo.rightCols(T - 1);
        U_0.col(T - 1) = Eigen::VectorXd::Zero(dim_u);
    }
}

Eigen::VectorXd MPPI::f(const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
        Eigen::VectorXd x_next(dim_x);
        Eigen::VectorXd f_dot(dim_x);
        double theta = x(4);
        Eigen::Matrix2d R;
        R << cos(theta), -sin(theta),
            sin(theta),  cos(theta);
        Eigen::Vector2d acc = g_ + (R * u) / mass;
        f_dot(0) = x(2);     
        f_dot(1) = x(3);      
        f_dot(2) = acc(0);    
        f_dot(3) = acc(1);    
        f_dot(4) = x(5);     
        f_dot(5) = -(l / I) * u(1);  
        x_next = x + dt * f_dot;
        return x_next;
    };

double MPPI::q(const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
    double control_cost = 2e-5 * u.squaredNorm();
    double state_cost = 5e-3 * (x.head(2).squaredNorm());
    return control_cost + state_cost;
};

double MPPI::p(const Eigen::VectorXd& x, const Eigen::VectorXd& x_target) {
    return 2000.0 * (x - x_target).norm();
};


void MPPI::h(Eigen::MatrixXd& U_seq) {
//------------------- conic constraint -------------------
// thrust u는 u(0) ≥ 0 이어야 하고, tan20 내에 있어야 함

    double input_angmax = tan(20.0 * M_PI / 180.0);
    for (int j = 0; j < U_seq.cols(); j++) {
        double u0 = U_seq(0, j);
        double u1 = U_seq(1, j);
        if (u0 < 0) {
            U_seq(0, j) = 0;
            u0 = 0;
        }
        double limit = input_angmax * u0;
        if (u1 > limit) {
            U_seq(1, j) = limit;
        } else if (u1 < -limit) {
            U_seq(1, j) = -limit;
        }
    }
};

void MPPI::saveCost(double cost) {
    iteration_costs.push_back(cost);
};

void MPPI::plotCost() {
    std::vector<double> iterations(iteration_costs.size());
    for (size_t i = 0; i < iteration_costs.size(); ++i) {
        iterations[i] = static_cast<double>(i + 1); 

    }

    // Plot
    plt::figure_size(800, 400);
    plt::plot(iterations, iteration_costs, {{"label", "Cost"}});
    // plt::yscale("log"); 
    plt::xlabel("Iteration");
    plt::ylabel("Cost");
    plt::legend();
    plt::show();
}

