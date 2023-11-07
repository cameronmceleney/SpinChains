# Start of CMakeLists.txt
cmake_minimum_required(VERSION 3.21)

# Project name here (SpinChains) should not be confused with the old local project (Chainspins) that is not available on Git
project(SpinChains)

# Written using C23. No real attempt has been made to ensure backwards compatibility
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Manually specifying Eigen and FFTW directories. May be different on your system. Numbers in PATHS are version numbers
if(WIN32 OR WIN64)
    # Windows requires manual install or use of Unix-based Homebrew
    set(EIGEN_DIR "C:/Program Files/Homebrew/eigen/3.4.0")
    set(FFTW_DIR "C:/Program Files/Homebrew/fftw/3.3.5")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3-3.a")
elseif(APPLE)
    # Default Homebrew install locations on MacOS
    set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
    set(FFTW_DIR "/usr/local/Cellar/fftw/3.3.10_1")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3.a")
endif()

# Create lists of all user sources
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

# Links all user source files to project.
add_executable(SpinChains main.cpp
        ${LIBRARY_HEADER_ONLY}
        ${LIBRARY_CLASSES}
        ${EXECUTABLE_CLASSES}
)

# Include directories for 3rd party libraries
include_directories(${EIGEN_DIR} ${FFTW_DIR}/include)

# FFTW requires linking to a static library
target_link_libraries(SpinChains PRIVATE ${FFTW_LIB})

# Optional: Uncomment for -O3 optimization
# target_compile_options(SpinChains PRIVATE -O3)

# End of CMakeLists.txt