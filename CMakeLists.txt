# Latest version of CMake is 3.21 (as of 20 Jan 21)
cmake_minimum_required(VERSION 3.21)

# Project name here (SpinChains) should not be confused with the old local project (Chainspins) that is not available on Git
project(SpinChains)

# Written using C23 and no real attempt has been made to ensure backwards compatibility
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(isStudent TRUE)

# Manually specifying Eigen and FFTW directories. May be different on your system
if(WIN32 OR WIN64)
    # Windows requires manual install or use of Unix-based Homebrew
    set(EIGEN_DIR "C:/Program Files/Homebrew/eigen/3.4.0")
    set(FFTW_DIR "C:/Program Files/Homebrew/fftw/3.3.5")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3-3.a")
elseif(APPLE)
    if(isStudent)
        # Default Homebrew install locations on MacOS. Jiahui Wang
        set(EIGEN_DIR "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
        set(FFTW_DIR "/opt/homebrew/Cellar/fftw/3.3.10_1")
        set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3.a")
    else()
        # Default Homebrew install locations on MacOS. Cameron McEleney
        set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
        set(FFTW_DIR "/usr/local/Cellar/fftw/3.3.10_1")
        set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3.a")
    endif()
endif()

# Create lists of all sources
set(LIBRARY_HEADER_ONLY
        CommonLibs.h
        progressbar.hpp
)

set(LIBRARY_CLASSES
        GlobalVariables.cpp GlobalVariables.h
        linspace.cpp linspace.h
        NewSolver.cpp NewSolver.h
)

set(EXECUTABLE_CLASSES
        Numerical_Methods_Class.cpp Numerical_Methods_Class.h
        SpinChainEigenSolverClass.cpp SpinChainEigenSolverClass.h
)

# Links all files to the project. There are numerous classes to encapsulate data so ensure that all .h files have a corresponding .cpp as needed
add_executable(SpinChains main.cpp
        ${LIBRARY_HEADER_ONLY}
        ${LIBRARY_CLASSES}
        ${EXECUTABLE_CLASSES}
)

# Include directories for FFTW
include_directories(${EIGEN_DIR} ${FFTW_DIR}/include)

# FFTW requires linking
target_link_libraries(SpinChains PRIVATE ${FFTW_LIB})

# Optional: Uncomment for -O3 optimization
# target_compile_options(SpinChains PRIVATE -O3)