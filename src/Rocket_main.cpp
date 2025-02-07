#include "mppiRocket.h"

int main() {
    int max_iter = 220;
    MPPI mppi;
    double final_cost = 0;

    for (int i = 0; i < max_iter; i++) {
        mppi.solve();
        mppi.move();
        final_cost = (mppi.x_init - mppi.x_target).head(2).norm();
        // std::cout << "err: " << final_cost << std::endl;
        std::cout << mppi.x_init.transpose() <<std::endl;
        // std::cout << mppi.x_init(1) <<std::endl;
        // std::cout << mppi.x_init <<std::endl;
        if (final_cost < 1) break;
    }
    mppi.plotCost();

    return 0;
}