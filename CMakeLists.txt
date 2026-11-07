cmake_minimum_required(VERSION 3.14)
project(planeperc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-Wall -Wextra)

# Header-only core library: voxel grid, KD-tree, SOR, RANSAC, LS fit, flatness (all from scratch)
add_library(planeperc INTERFACE)
target_include_directories(planeperc INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

add_executable(synthetic_eval apps/synthetic_eval.cpp)
target_link_libraries(synthetic_eval planeperc)

add_executable(unit_tests tests/unit_tests.cpp)
target_link_libraries(unit_tests planeperc)
enable_testing()
add_test(NAME unit_tests COMMAND unit_tests)

# ROS 2 node: only built when sourced inside a ROS 2 environment
find_package(ament_cmake QUIET)
if(ament_cmake_FOUND)
  find_package(rclcpp REQUIRED)
  find_package(sensor_msgs REQUIRED)
  add_executable(plane_node ros/plane_node.cpp)
  target_link_libraries(plane_node planeperc)
  ament_target_dependencies(plane_node rclcpp sensor_msgs)
  install(TARGETS plane_node DESTINATION lib/${PROJECT_NAME})
  ament_package()
endif()
