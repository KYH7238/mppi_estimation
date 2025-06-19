#include "droneLocalization.h"
int main(int argc, char **argv) {
    ros::init(argc, argv, "ekf_node");
    EKF ekf;
    ros::spin();
    return 0;
}