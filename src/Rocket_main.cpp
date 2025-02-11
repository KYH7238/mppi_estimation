#include "mppiRocket.h"
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <limits.h>

int main() {
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        perror("getcwd() error");
    }

    int max_iter = 1000;
    MPPI mppi;
    double final_cost = 0;

    std::ofstream outfile("/home/kim/drone_ws/src/mppi_estimation/data/output3.xls", std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "Error opening /tmp/output.xls for writing!" << std::endl;
        return -1;
    }
    
    outfile << "Y\tX\tVY\tVX\tW\tA\tUY\tUX\n";
    
    for (int i = 0; i < max_iter; i++) {
        mppi.solve();
        mppi.move();
        final_cost = (mppi.x_init - mppi.x_target).head(2).norm();
        std::cout << mppi.x_init.transpose() << std::endl;
        
        for (int j = 0; j < 6; j++) {
            outfile << mppi.x_init(j) << "\t";
        }
        outfile << mppi.u0(0) << "\t" << mppi.u0(1) << "\n";
        
        if (mppi.x_init(1) < 0.0)
            break;
    }
    
    outfile.close();
    mppi.plotCost();

    return 0;
}
