cmake_minimum_required(VERSION 3.5)
project(mppi_estimation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  nlink_parser
  geometry_msgs
  std_msgs
  sensor_msgs
)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED) 

catkin_package(
  CATKIN_DEPENDS roscpp nlink_parser rospy
)

include_directories(
  include
  /home/kim/EigenRand                    
  /usr/include/matplotlib-cpp            
  ${EIGEN3_INCLUDE_DIRS}                  
  ${Python3_INCLUDE_DIRS}                 
  ${Python3_NumPy_INCLUDE_DIRS}           
  ${catkin_INCLUDE_DIRS}                   
)

# add_executable(mppi_rocket
#   src/Rocket_main.cpp
#   src/mppiRocket.cpp
# )

# add_executable(estimation
#   src/main.cpp
#   src/mppiEstimation.cpp
# )

add_executable(mppi_node src/main.cpp src/mppiEstimation.cpp)
target_link_libraries(mppi_node ${catkin_LIBRARIES})
add_dependencies(mppi_node ${catkin_EXPORTED_TARGETS})

add_executable(mppi_node_2tag src/main_2tag.cpp src/mppiEstimation_2tag.cpp)
target_link_libraries(mppi_node_2tag ${catkin_LIBRARIES})
add_dependencies(mppi_node_2tag ${catkin_EXPORTED_TARGETS})

add_executable(ekf_node src/main_ekf.cpp src/droneLocalization.cpp)
target_link_libraries(ekf_node ${catkin_LIBRARIES})
add_dependencies(ekf_node ${catkin_EXPORTED_TARGETS})

target_link_libraries(mppi_node
  ${catkin_LIBRARIES}
  OpenMP::OpenMP_CXX
  ${Python3_LIBRARIES}
  ${Python3_NumPy_LIBRARIES}
)

target_link_libraries(mppi_node_2tag
  ${catkin_LIBRARIES}
  OpenMP::OpenMP_CXX
  ${Python3_LIBRARIES}
  ${Python3_NumPy_LIBRARIES}
)
