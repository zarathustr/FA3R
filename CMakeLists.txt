cmake_minimum_required(VERSION 2.6)
project(fa3r)

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

add_executable(test FA3R.cpp main.cpp)



