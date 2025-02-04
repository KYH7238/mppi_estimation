#include "mppiRocket.h"

int main() {
    int max_iter = 30;
    MPPI mppi;
    double final_cost = 0;

    for (int i = 0; i < max_iter; i++) {
        mppi.solve();
        mppi.move();
        // final_cost = (mppi.x_init - mppi.x_target).norm();
        // std::cout << "norm: " << final_cost << std::endl;
        // std::cout << mppi.x_init <<std::endl;
        // if ((mppi.x_init - mppi.x_target).norm() < 70) break;
    }

    // std::cout << "Final state: " << mppi.x_init.transpose() << std::endl;

    mppi.plotCost();

    return 0;
}