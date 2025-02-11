// #include "mppiRocket.h"
// namespace plt = matplotlibcpp;
// MPPI::MPPI() {
//     N = 500;
//     dim_x = 6;
//     dim_u = 2;
//     dim_g = 1;    
//     dim_h = 2;  
//     l = 0.7; 
//     dt = 0.01;
//     mass = 10.0;
//     I = 3.33;
//     g_ << -9.81, 0.0;
//     F_min = 21.58;
//     F_max = 107.91;
//     T = 50;
//     gamma_u = 1.0;
//     sigma_u = 10 * Eigen::MatrixXd::Identity(dim_u, dim_u); 
//     U_0 = 9.81 * 10.0 * Eigen::MatrixXd::Ones(dim_u, T);
//     Xo = Eigen::MatrixXd::Zero(dim_x, T+1);
//     x_init = Eigen::VectorXd::Zero(dim_x);
//     x_init(0) = 8.0;
//     x_init(1) = 6.0;
//     C = 0.18;
//     x_target = Eigen::VectorXd::Zero(dim_x);
// }

// MPPI::~MPPI() { }

// Eigen::MatrixXd MPPI::getNoise(const int &T) {
//     return sigma_u * norm_gen.template generate<Eigen::MatrixXd>(dim_u, T, urng);
// }

// void MPPI::solve() {
//     auto start = std::chrono::high_resolution_clock::now();

//     Eigen::MatrixXd Ui = U_0.replicate(N, 1);
//     Eigen::VectorXd costs = Eigen::VectorXd::Zero(N);
//     Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);
//     #pragma omp parallel for
//     for (int i = 0; i < N; i++) {
//         Eigen::MatrixXd noise = getNoise(T);
//         // std::cout << noise << std::endl;
//         Ui.middleRows(i * dim_u, dim_u) += noise;
//         h(Ui.middleRows(i * dim_u, dim_u));
//         Eigen::MatrixXd Xi = Eigen::MatrixXd::Zero(dim_x, T+1);
//         Xi.col(0) = x_init;
//         double cost = 0.0;

//         for (int j = 0; j < T; j++) {
//             cost += q(Xi.col(j), Ui.block(i * dim_u, j, dim_u, 1));
//             Xi.col(j+1) = f(Xi.col(j), Ui.block(i * dim_u, j, dim_u, 1));
//         }

//         cost += p(Xi.col(T), x_target);
//         costs(i) = cost;
//     }

//     double min_cost = costs.minCoeff();
//     weights = (-gamma_u * (costs.array() - min_cost)).exp();    
//     double total_weight = weights.sum();
//     weights /= total_weight;

//     Uo = Eigen::MatrixXd::Zero(dim_u, T);
//     for (int i = 0; i < N; i++) {
//         Uo += Ui.middleRows(i * dim_u, dim_u) * weights(i);
//     }

//     auto finish = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed = finish - start;
//     std::cout << "time: " << elapsed.count() << " sec" << std::endl;

//     u0 = Uo.col(0);
//     Xo.col(0) = x_init;
//     for (int j = 0; j < T; j++) {
//         Xo.col(j+1) = f(Xo.col(j), Uo.col(j));
//     }
//     saveCost(min_cost);
// }

// void MPPI::move() {
//     x_init = f(x_init, u0);
//     if (T > 1) {
//         U_0.leftCols(T - 1) = Uo.rightCols(T - 1);
//     }
// }

// Eigen::VectorXd MPPI::f(const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
//         Eigen::VectorXd x_next(dim_x);
//         Eigen::VectorXd f_dot(dim_x);
//         double theta = x(4);
//         Eigen::Matrix2d R;
//         R << cos(theta), -sin(theta),
//              sin(theta),  cos(theta);
//         Eigen::Vector2d acc = g_ + (R * u) / mass;
//         f_dot(0) = x(2);     
//         f_dot(1) = x(3);      
//         f_dot(2) = acc(0);    
//         f_dot(3) = acc(1);    
//         f_dot(4) = x(5);     
//         f_dot(5) = -(l / I) * u(1);  
//         x_next = x + dt * f_dot;
//         return x_next;
//     };

// double MPPI::q(const Eigen::VectorXd& x, const Eigen::VectorXd& u) {

//     double angle = std::atan2(x(1), x(0));
//     double alpha_max = std::tan(45.0 * M_PI / 180.0);
//     double conic_cost = std::fabs(angle - alpha_max);
//     double la = 10;
//     double control_cost = 2e-5 * u.squaredNorm();
//     double state_cost = 5e-3 * (x(1)*x(1));
//     // double state_cost = 5e-3 * (x.head(2).squaredNorm());
//     double alpha = std::atan2(u(1), u(0))*180/M_PI;  
//     double alpha_cost = C * alpha* alpha;
//     return control_cost + state_cost + alpha_cost + la*conic_cost;
//     // return state_cost ;
// };

// double MPPI::p(const Eigen::VectorXd& x, const Eigen::VectorXd& x_target) {
//     Eigen::Vector2d p_diff = x.head(2) - x_target.head(2); 
//     Eigen::Vector2d v_diff = x.segment(2, 2) - x_target.segment(2, 2); 
//     double theta_diff = x(4) - x_target(4); 
//     double w_diff = x(5) - x_target(5);
//     return 60 * (p_diff.squaredNorm() + v_diff.squaredNorm() + theta_diff * theta_diff + w_diff * w_diff);
// }

// void MPPI::h(Eigen::Ref<Eigen::MatrixXd> U_seq) {
//     const double t = std::tan(20.0 * M_PI / 180.0);
//     for (int j = 0; j < U_seq.cols(); ++j) {
//         double u0 = U_seq(0, j);
//         double u1 = U_seq(1, j);
        
//         if(u0 < F_min) {
//             u0 = F_min;
//         } else if(u0 > F_max) {
//             u0 = F_max;
//         }
//         double a = u0;
//         double v = u1 / t;
        
//         if (std::fabs(v) > a) {  
//             double factor = 0.5 * (1.0 + a / std::fabs(v));
//             v = factor * v;
//             a = factor * std::fabs(v);  
//         }
//         U_seq(0, j) = a;
//         U_seq(1, j) = t * v;
//     }
// }


// void MPPI::saveCost(double cost) {
//     iteration_costs.push_back(cost);
// };

// void MPPI::plotCost() {
//     std::vector<double> iterations(iteration_costs.size());
//     for (size_t i = 0; i < iteration_costs.size(); ++i) {
//         iterations[i] = static_cast<double>(i + 1); 
//     }

//     plt::figure_size(800, 400);
//     plt::plot(iterations, iteration_costs, {{"label", "Cost"}});
//     plt::xlabel("Iteration");
//     plt::ylabel("Cost");
//     plt::legend();
//     plt::show();
// }

#include "mppiRocket.h"
namespace plt = matplotlibcpp;
MPPI::MPPI() {
    N = 500;
    dim_x = 6;
    dim_u = 2;
    dim_g = 1;    
    dim_h = 2;  
    l = 0.7; 
    dt = 0.01;
    mass = 10.0;
    I = 3.33;
    g_ << -9.81, 0.0;
    F_min = 21.58;
    F_max = 107.91;
    T = 50;
    gamma_u = 1.0;
    sigma_u = 10 * Eigen::MatrixXd::Identity(dim_u, dim_u); 
    U_0 = 9.81 * 10.0 * Eigen::MatrixXd::Ones(dim_u, T);
    Xo = Eigen::MatrixXd::Zero(dim_x, T+1);
    x_init = Eigen::VectorXd::Zero(dim_x);
    x_init(0) = 8.0;
    x_init(1) = 6.0;
    C = 0.18;
    x_target = Eigen::VectorXd::Zero(dim_x);
}

MPPI::~MPPI() { }

Eigen::MatrixXd MPPI::getNoise(const int &T) {
    return sigma_u * norm_gen.template generate<Eigen::MatrixXd>(dim_u, T, urng);
}

void MPPI::solve() {
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::MatrixXd Ui = U_0.replicate(N, 1);
    Eigen::VectorXd costs = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd weights = Eigen::VectorXd::Zero(N);
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        Eigen::MatrixXd noise = getNoise(T);
        // std::cout << noise << std::endl;
        Ui.middleRows(i * dim_u, dim_u) += noise;
        h(Ui.middleRows(i * dim_u, dim_u));
        Eigen::MatrixXd Xi = Eigen::MatrixXd::Zero(dim_x, T+1);
        Xi.col(0) = x_init;
        double cost = 0.0;

        for (int j = 0; j < T; j++) {
            cost += q(Xi.col(j), Ui.block(i * dim_u, j, dim_u, 1));
            Xi.col(j+1) = f(Xi.col(j), Ui.block(i * dim_u, j, dim_u, 1));
        }

        cost += p(Xi.col(T), x_target);
        costs(i) = cost;
    }

    double min_cost = costs.minCoeff();
    weights = (-gamma_u * (costs.array() - min_cost)).exp();    
    double total_weight = weights.sum();
    weights /= total_weight;

    Uo = Eigen::MatrixXd::Zero(dim_u, T);
    for (int i = 0; i < N; i++) {
        Uo += Ui.middleRows(i * dim_u, dim_u) * weights(i);
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "time: " << elapsed.count() << " sec" << std::endl;

    u0 = Uo.col(0);
    Xo.col(0) = x_init;
    for (int j = 0; j < T; j++) {
        Xo.col(j+1) = f(Xo.col(j), Uo.col(j));
    }
    saveCost(min_cost);
}

void MPPI::move() {
    x_init = f(x_init, u0);
    if (T > 1) {
        U_0.leftCols(T - 1) = Uo.rightCols(T - 1);
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

    double angle = std::atan2(x(1), x(0));
    double alpha_max = std::tan(45.0 * M_PI / 180.0);
    double conic_cost = std::fabs(angle - alpha_max);
    double la = 10;
    double control_cost = 2e-5 * u.squaredNorm();
    double state_cost = 5e-3 * (x(1)*x(1));
    // double state_cost = 5e-3 * (x.head(2).squaredNorm());
    double alpha = std::atan2(u(1), u(0))*180/M_PI;  
    double alpha_cost = C * alpha* alpha;
    return control_cost + state_cost + alpha_cost + la*conic_cost;
    // return state_cost ;
};

double MPPI::p(const Eigen::VectorXd& x, const Eigen::VectorXd& x_target) {
    Eigen::Vector2d p_diff = x.head(2) - x_target.head(2); 
    Eigen::Vector2d v_diff = x.segment(2, 2) - x_target.segment(2, 2); 
    double theta_diff = x(4) - x_target(4); 
    double w_diff = x(5) - x_target(5);
    return 53 * (p_diff.squaredNorm() + v_diff.squaredNorm() + theta_diff * theta_diff + w_diff * w_diff);
}

void MPPI::h(Eigen::Ref<Eigen::MatrixXd> U_seq) {
    const double t = std::tan(20.0 * M_PI / 180.0);
    for (int j = 0; j < U_seq.cols(); ++j) {
        double u0 = U_seq(0, j);
        double u1 = U_seq(1, j);
        
        if(u0 < F_min) {
            u0 = F_min;
        } else if(u0 > F_max) {
            u0 = F_max;
        }
        double a = u0;
        double v = u1 / t;
        
        if (std::fabs(v) > a) {  
            double factor = 0.5 * (1.0 + a / std::fabs(v));
            v = factor * v;
            a = factor * std::fabs(v);  
        }
        U_seq(0, j) = a;
        U_seq(1, j) = t * v;
    }
}


void MPPI::saveCost(double cost) {
    iteration_costs.push_back(cost);
};

void MPPI::plotCost() {
    std::vector<double> iterations(iteration_costs.size());
    for (size_t i = 0; i < iteration_costs.size(); ++i) {
        iterations[i] = static_cast<double>(i + 1); 
    }

    plt::figure_size(800, 400);
    plt::plot(iterations, iteration_costs, {{"label", "Cost"}});
    plt::xlabel("Iteration");
    plt::ylabel("Cost");
    plt::legend();
    plt::show();
}