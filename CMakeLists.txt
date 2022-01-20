# Latest version of CMake is 3.21 (as of 20 Jan 21)
cmake_minimum_required(VERSION 3.21)
# Project name here (SpinChains) should not be confused with the old local project (Chainspins) that is not available on Git
project(SpinChains)

# Written using C23 and no real attempt has been made to ensure backwards compatibility
set(CMAKE_CXX_STANDARD 23)

# Two IF statements test to find the user's Operating System, and then dictate the path to the Eigen folder. It was too much
# work to have relative paths to a cloud storage location that would hold the Eigen folder (And thus have the same relative path
# for both OS. This was the quick yet effective solution
if(WIN32 OR WIN64)
    # This is the absolute path to my Eigen folder on Windows; change as required.
    set(EIGEN_DIR "C:/Program Files/eigen-3.4.0")
endif()

if(APPLE)
    # This is the absolute path to my Eigen folder on Mac (used Homebrew to install_; change as required
    set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
endif()

# Includes the Eigen library. This lib is full of header files, so they should be called as such (eg #include <Eigen/Core>)
include_directories(${EIGEN_DIR})

# Links all files to the project. There are numerous classes to encapsulate data so ensure that all .h files have a corresponding .cpp as needed
add_executable(SpinChains
        main.cpp
        CommonLibs.h
        linspace.cpp linspace.h
        SpinChainEigenSolverClass.cpp SpinChainEigenSolverClass.h
        Numerical_Methods_Class.cpp Numerical_Methods_Class.h
        GlobalVariables.cpp GlobalVariables.h)
