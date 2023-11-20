# Latest version of CMake is 3.21 (as of 20 Jan 21)
cmake_minimum_required(VERSION 3.21)
# Project name here (SpinChains) should not be confused with the old local project (Chainspins) that is not available on Git
project(SpinChains)

# Written using C23 and no real attempt has been made to ensure backwards compatibility
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Manually specifying Eigen and FFTW directories
if(WIN32 OR WIN64)
    set(EIGEN_DIR "C:/Program Files/Homebrew/eigen/3.4.0")
    set(FFTW_DIR "C:/Program Files/Homebrew/fftw/3.3.5")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3-3.a")
elseif(APPLE)
    set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
    set(FFTW_DIR "/usr/local/Cellar/fftw/3.3.10_1")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3.a")
endif()

#include_directories(${EIGEN_DIR} ${FFTW_DIR})

# Links all files to the project. There are numerous classes to encapsulate data so ensure that all .h files have a corresponding .cpp as needed
add_executable(SpinChains src/main.cpp
        src/CommonLibs.h
        src/linspace.cpp include/linspace.h
        src/SpinChainEigenSolverClass.cpp include/SpinChainEigenSolverClass.h
        #NMSuperClass.cpp NMSuperClass.h
        src/NMSuperClassTest.cpp include/NMSuperClassTest.h
        NMSubClasses/NMInitialisation.cpp NMSubClasses/NMInitialisation.h
        NMSubClasses/NMMethods.cpp NMSubClasses/NMMethods.h
        NMSubClasses/NMDataHandling.cpp NMSubClasses/NMDataHandling.h
        NMSubClasses/NMConfiguration.cpp NMSubClasses/NMConfiguration.h
        NMSubClasses/SystemDataContainer.h
        src/GlobalVariables.cpp include/GlobalVariables.h
        include/progressbar.hpp
        src/NewSolver.cpp include/NewSolver.h
        NMSubClasses/DipolarField.cpp NMSubClasses/DipolarField.h
        NMSubClasses/EffectiveField.cpp NMSubClasses/EffectiveField.h
        NMSubClasses/DemagField.cpp NMSubClasses/DemagField.h
        NMSubClasses/LLG.cpp NMSubClasses/LLG.h
        NMSubClasses/IDemagField.h
        NMSubClasses/INMMethods.h
        NMSubClasses/IDipolarField.h
        NMSubClasses/IEffectiveField.h
        NMSubClasses/ILLG.h
)

# Include directories for FFTW
include_directories(${EIGEN_DIR} ${FFTW_DIR}/include)

# Link FFTW to the executable
target_link_libraries(SpinChains PRIVATE ${FFTW_LIB})

# Optional: Uncomment for -O3 optimization
# target_compile_options(SpinChains PRIVATE -O3)