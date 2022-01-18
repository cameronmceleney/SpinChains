cmake_minimum_required(VERSION 3.21)
project(SpinChains)

set(CMAKE_CXX_STANDARD 23)

if(WIN32 OR WIN64)
    #This is the absolute path to my Eigen folder on Windows. You should change line7 to be the path to yours
    set(EIGEN_DIR "C:/Program Files/eigen-3.4.0")
endif()

if(APPLE)
    #This is the absolute path to my Eigen folder on Mac; the default set by homebrew
    set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
endif()

include_directories(${EIGEN_DIR})

add_executable(SpinChains
        main.cpp
        common.h
        linspace.cpp linspace.h
        matrix_operations.cpp matrix_operations.h
        SpinChainEigenSolverClass.cpp SpinChainEigenSolverClass.h
        Numerical_Methods_Class.cpp Numerical_Methods_Class.h)
