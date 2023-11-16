# Start of CMakeLists.txt
cmake_minimum_required(VERSION 3.17)

# Project name here (SpinChains) should not be confused with the old local project (Chainspins) that is not available on Git
project(SpinChains)

# Written using C23. No real attempt has been made to ensure backwards compatibility
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#if(CMAKE_BUILD_TYPE STREQUAL "Release (MinGW64)")
#    set(USE_RELMINGW ON)
#else()
#    set(USE_RELMINGW OFF)
#endif()

# Detect if using Clang
#if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#  # Add flags for Clang
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdenormal-fp-math=positive-zero")
#endif()

# Flags for optimisation and to make static executables
#,set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -static-libstdc++ -static-libgcc")

# Manually specifying Eigen and FFTW directories. May be different on your system. Numbers in PATHS are version numbers
if(WIN32 OR WIN64)
    # Windows requires manual install or use of Unix-based Homebrew
    set(EIGEN_DIR "C:/Program Files/Homebrew/eigen/3.4.0")
    set(FFTW_DIR "C:/Program Files/Homebrew/fftw/3.3.5")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3-3.a")

    # Set the path to the TBB library (replace with the actual path where TBB is installed)
    set(TBB_BASE "C:/msys64/mingw64")
    set(TBB_DIR "${TBB_BASE}/include/oneapi")
    set(TBB_LIBRARIES "${TBB_BASE}/lib/libtbb12.dll.a" "${TBB_BASE}/lib/libtbbmalloc.dll.a")
elseif(APPLE)
    # Default Homebrew install locations on MacOS
    set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
    set(FFTW_DIR "/usr/local/Cellar/fftw/3.3.10_1")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3.a")

    # Set the path to the TBB library (replace with the actual path where TBB is installed)
    set(TBB_BASE "/usr/local/Cellar/tbb/2021.10.0")
    set(TBB_DIR "${TBB_BASE}/include/oneapi")
    set(TBB_LIBRARIES "${TBB_BASE}/lib/libtbb.a" "${TBB_BASE}/lib/libtbbmalloc.a")
endif()

# Create lists of all user sources
set(LIBRARY_HEADER_ONLY
        src/CommonLibs.h
        include/progressbar.hpp
)

set(LIBRARY_CLASSES
        src/GlobalVariables.cpp include/GlobalVariables.h
        src/linspace.cpp include/linspace.h
        src/NewSolver.cpp include/NewSolver.h
)

set(EXECUTABLE_CLASSES
        src/Numerical_Methods_Class.cpp include/Numerical_Methods_Class.h
        src/SpinChainEigenSolverClass.cpp include/SpinChainEigenSolverClass.h
)

set(TBB_CLASSES
        src/TestClass.cpp include/TestClass.h
)

#set(TBB_ENABLE_IPO OFF CACHE BOOL "Disable interprocedural optimization for TBB")

# Links all user source files to project.
add_executable(${PROJECT_NAME} src/main.cpp
        ${LIBRARY_HEADER_ONLY}
        ${LIBRARY_CLASSES}
        ${EXECUTABLE_CLASSES}
        ${TBB_CLASSES}
)

# Include directories for 3rd party libraries

# FFTW requires linking to a static library

include_directories(${EIGEN_DIR} ${FFTW_DIR}/include ${TBB_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LIB} ${TBB_LIBRARIES})

# Optional: Uncomment for -O3 optimization
# target_compile_options(SpinChains PRIVATE -O0)

# End of CMakeLists.txt