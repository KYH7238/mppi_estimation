#include "mppiRocket.h"

int main() {
    int max_iter = 200;
    MPPI mppi;
    double final_cost = 0;

    for (int i = 0; i < max_iter; i++) {
        if (i >=120) {
            mppi.dt = 0.001;
        }
        mppi.solve();
        // mppi.move();
        final_cost = (mppi.x_init - mppi.x_target).head(2).norm();
        // std::cout << "err: " << final_cost << std::endl;
        std::cout << "iter: "<< i <<"\tstate: "<<mppi.x_init.head(2).transpose() <<std::endl;
        // std::cout << mppi.x_init(1) <<std::endl;
        // std::cout << mppi.x_init <<std::endl;
        if (mppi.x_init(0) < 1.0) break;
    }
    
    mppi.plotCost();
    return 0;
}