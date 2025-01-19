cmake_minimum_required(VERSION 3.5)
project(mppi_estimation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  nlink_parser
  geometry_msgs
  roscpp
  std_msgs
  rospy
  sensor_msgs
)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
find_package(PythonLibs REQUIRED)  


catkin_package(
  CATKIN_DEPENDS roscpp nlink_parser
)

include_directories(
  include
  /home/kim/EigenRand  
  ${EIGEN3_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
)

add_executable(mppi_estimationT
  src/main.cpp
  src/mppiEstimation.cpp
)

target_link_libraries(mppi_estimationT
  ${catkin_LIBRARIES}        
  OpenMP::OpenMP_CXX        
)
