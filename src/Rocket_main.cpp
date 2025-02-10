#include "mppiRocket.h"
#include <ros/ros.h>

int main(int argc, char **argv)
{
    ros::init(argc, argv, "mppi_estimation_node");
    ros::NodeHandle nh;

    int max_iter = 560;
    MPPI mppi(nh);
    double final_cost = 0;

    for (int i = 0; i < max_iter && ros::ok(); i++)
    {
        // if (i >= 120) {
        //     mppi.dt = 0.001;
        // }
        mppi.solve();  
        // mppi.move();   
        // mppi.publishTrajectory(mppi.Xo);

        final_cost = (mppi.x_init - mppi.x_target).head(2).norm();
        std::cout <<"iter: " << i << "\tstate: " <<mppi.x_init.head(2).transpose() << std::endl;
        // ROS_INFO_STREAM("Iteration " << i << " - Current state: " << mppi.x_init.transpose());
        ros::spinOnce();
    }
    // mppi.plotCost();
    return 0;
}