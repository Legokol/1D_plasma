cmake_minimum_required(VERSION 3.23)
project(1D_plasma)

set(CMAKE_CXX_STANDARD 17)

add_executable(1D_plasma main.cpp ES1.cpp ES1.hpp Constants.hpp FourierTransform.hpp)
