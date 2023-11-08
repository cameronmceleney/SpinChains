# Start of CMakeLists.txt
cmake_minimum_required(VERSION 3.17)

# Project name here (SpinChains) should not be confused with the old local project (Chainspins) that is not available on Git
project(SpinChains)

# Written using C23. No real attempt has been made to ensure backwards compatibility
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_BUILD_TYPE STREQUAL "Release (MinGW64)")
    set(USE_RELMINGW ON)
else()
    set(USE_RELMINGW OFF)
endif()



# Manually specifying Eigen and FFTW directories. May be different on your system. Numbers in PATHS are version numbers
if(WIN32 OR WIN64)
    # Windows requires manual install or use of Unix-based Homebrew
    set(EIGEN_DIR "C:/Program Files/Homebrew/eigen/3.4.0")
    set(FFTW_DIR "C:/Program Files/Homebrew/fftw/3.3.5")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3-3.a")

    # Set the path to the TBB library (replace with the actual path where TBB is installed)
    set(TBB_DIR "C:/msys64/mingw64")
    # Find the TBB package
    if (USE_RELMINGW)
        find_package(TBB REQUIRED PATHS ${TBB_DIR} NO_DEFAULT_PATH)
        if(NOT TBB_FOUND)
            set(TBB_LIBRARIES "${TBB_DIR}/lib/libtbb.dll.a" "${TBB_DIR}/lib/libtbbmalloc.dll.a")
            include_directories("${TBB_DIR}/include")
        endif()
    endif()
elseif(APPLE)
    # Default Homebrew install locations on MacOS
    set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
    set(FFTW_DIR "/usr/local/Cellar/fftw/3.3.10_1")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3.a")
    set(TBB_DIR "/usr/local/Cellar/tbb/3.3.10_1")
    set(TBB_LIB "${TBB_DIR}/lib/libtbb.a")
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

set(TBB_CLASSES
        TestClass.cpp TestClass.h
)

#set(TBB_ENABLE_IPO OFF CACHE BOOL "Disable interprocedural optimization for TBB")

# Links all user source files to project.
add_executable(${PROJECT_NAME} main.cpp
        ${LIBRARY_HEADER_ONLY}
        ${LIBRARY_CLASSES}
        ${EXECUTABLE_CLASSES}
        ${TBB_CLASSES}
)

# Include directories for 3rd party libraries

# FFTW requires linking to a static library
if (USE_RELMINGW)
    if(TBB_FOUND)
        include_directories(${EIGEN_DIR} ${FFTW_DIR}/include)
        target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LIB} TBB::tbb)
    else()
        include_directories(${EIGEN_DIR} ${FFTW_DIR}/include "${TBB_DIR}/include")
        target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LIB} ${TBB_LIBRARIES})
    endif()
else()
    include_directories(${EIGEN_DIR} ${FFTW_DIR}/include)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LIB})
endif()

# Optional: Uncomment for -O3 optimization
# target_compile_options(SpinChains PRIVATE -O3)

# End of CMakeLists.txt