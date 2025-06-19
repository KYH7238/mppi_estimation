# For offline MPPI code (without ROS)

    src/main_offline.cpp
    src/offlineMppiEstimation_2tag.cpp
    
    include/offlineMppiEstimation_2tag.h

    scripts/bag2txt.py
    scripts/plot.py
    
    data/data1_gt.txt
    data/imu.txt
    data/mppi_pose.txt
    data/uwb.txt
    data/imu.txt
    
    config/anchors.txt

---
> g++ -std=c++17 -O3 -march=native -fopenmp -I../include -I/usr/include/eigen3 -I/home/kim/EigenRand offlineMppiEstimation_2tag.cpp main_offline.cpp -o mppi

> ./mppi

# For online MPPI code (with ROS1)
    scripts/new_rpy.py
    scripts/2dplot.py
## 1 Tag equipped
    src/main.cpp
    src/mainEstimation.cpp
    
    include/Estimation.h

    scripts/plot_realtime.py    

    config/0612_hw1_gt.txt
    config/0612_hw2_gt.txt
    config/0612_hw3_gt.txt

    bag/20240612/0612_hw1.bag
    bag/20240612/0612_hw2.bag
    bag/20240612/0612_hw3.bag
> roslaunch mppi_esimation mppi_estimation.launch

## 2 Tag equipped
    src/main_2tag.cpp
    src/mainEstimation_2tag.cpp
    
    include/Estimation_2tag.h

    scripts/plot_realtime.py    

    bag/20250502/hw1.bag
    bag/20250502/hw2.bag
> roslaunch mppi_esimation mppi_estimation_2tag.launch

# For localization (EKF)
    src/droneLocalization.cpp
    src/main_ekf.cpp

    include/droneLocalization.h

> roslaunch mppi_estimation ekf_localization.launch