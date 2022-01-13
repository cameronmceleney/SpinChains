cmake_minimum_required(VERSION 3.21)
project(SpinChains)

set(CMAKE_CXX_STANDARD 23)

#This is the absolute path to my Eigen folder on Windows. You should change line7 to be the path to yours
#set(EIGEN_DIR "C:\\Users\\Cameron McEleney\\OneDrive - University of Glasgow\\University\\PhD\\1st Year\\C++\\eigen-3.4.0")
#This is the absolute path to my Eigen folder on Mac; the default set by homebrew
set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
include_directories(${EIGEN_DIR})

add_executable(SpinChains main.cpp saveMatrix.h common.h linspace.cpp linspace.h saveMatrix.cpp)
