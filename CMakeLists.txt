################################################################

# Start of CMakeLists.txt

################################################################
# Latest version of CMake is 3.21 (as of 20 Jan 21)SolversSuperClass.
cmake_minimum_required(VERSION 3.17)
project(SpinChains)  # This project's name is not to be confused with (ChainSpins); an old local redundant project

set(CMAKE_CXX_STANDARD 23)  # No attempt has been made to ensure backwards compatibility from C23
set(CMAKE_CXX_STANDARD_REQUIRED ON)

################################################################
# if(CMAKE_BUILD_TYPE STREQUAL "Release (MinGW64)")
#     set(USE_RELMINGW ON)
# else()
#     set(USE_RELMINGW OFF)
# endif()

# # Detect if using Clang when debugging underflow/overflow
# if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdenormal-fp-math=positive-zero")
# endif()

# Manually specifying Eigen and FFTW directories. May be different on your system.
# Numbers in PATHS are version numbers

################################################################
# Operating system specific settings
if(WIN32 OR WIN64)
    # Windows requires manual installation, or use of UNIX-based Homebrew
    set(EIGEN_DIR "C:/Program Files/Homebrew/eigen/3.4.0")
    set(FFTW_DIR "C:/Program Files/Homebrew/fftw/3.3.5")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3-3.a")

    # Use TBB (not OpenMP) for parallelisation. Use the newer oneAPI version
    set(TBB_BASE "C:/msys64/mingw64")
    set(TBB_DIR "${TBB_BASE}/include/oneapi")
    set(TBB_LIB "${TBB_BASE}/lib/libtbb12.dll.a" "${TBB_BASE}/lib/libtbbmalloc.dll.a")
elseif(APPLE)
    # Default Homebrew installation directories on MacOS
    set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
    set(FFTW_DIR "/usr/local/Cellar/fftw/3.3.10_1")
    set(FFTW_LIB "${FFTW_DIR}/lib/libfftw3.a")

    # Use TBB (not OpenMP) for parallelisation. Use the newer oneAPI version
    set(TBB_BASE "/usr/local/Cellar/tbb/2021.10.0")
    set(TBB_DIR "${TBB_BASE}/include/oneapi")
    set(TBB_LIB "${TBB_BASE}/lib/libtbb.a" "${TBB_BASE}/lib/libtbbmalloc.a")
endif()

################################################################
# Create lists of all user sources
set(USER_GENERAL_LIBRARIES
        include/CommonLibs.h
        libs/progressbar.hpp
        libs/linspace.cpp libs/linspace.h
)

set(USER_OTHER_FILES
        Other/working_on/SpinChainEigenSolverClass.cpp Other/working_on/SpinChainEigenSolverClass.h
        Other/Legacy_Files/NewSolver.cpp Other/Legacy_Files/NewSolver.h
        test/IntelTBBTestIfWorking.cpp test/IntelTBBTestIfWorking.h
)

set(CONTAINER_LIBRARIES
        include/SimulationFlags.h
        include/SimulationParameters.h
        include/SimulationStates.h
        src/GlobalVariables.cpp include/GlobalVariables.h
)

set(SOLVERS_CLASS
        src/SolversSuperClass.cpp include/SolversSuperClass.h
        src/SolversInitialisation.cpp include/SolversInitialisation.h
        src/SolversConfiguration.cpp include/SolversConfiguration.h
        src/SolversImplementation.cpp include/SolversImplementation.h include/InterfaceSolversImplementation.h
        src/SolversDataHandling.cpp include/SolversDataHandling.h
)

set(SOLVER_CLASS_COMPONENTS
        src/DemagnetisationFields.cpp include/DemagnetisationFields.h
        src/DipolarFields.cpp include/DipolarFields.h
        src/EffectiveFields.cpp include/EffectiveFields.h
        src/MagnetisationDynamics.cpp include/MagnetisationDynamics.h
        src/DzyaloshinskiiMoriyaInteraction.cpp
        include/DzyaloshinskiiMoriyaInteraction.h
        src/SpinTransferTorque.cpp
        include/SpinTransferTorque.h
)

# Must be before add_executable()
# set(TBB_ENABLE_IPO OFF CACHE BOOL "Disable interprocedural optimization for TBB")

# Links all user sources to the executable
add_executable(${PROJECT_NAME} src/main.cpp
        ${USER_GENERAL_LIBRARIES}
        ${USER_OTHER_FILES}
        ${CONTAINER_LIBRARIES}
        ${SOLVERS_CLASS}
        ${SOLVER_CLASS_COMPONENTS}
)

################################################################
# Include and link directories for third-party libraries
# Both FFTW and TBB require linking to static libraries; Eigen does not.

include_directories(${EIGEN_DIR} ${FFTW_DIR}/include ${TBB_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LIB} ${TBB_LIB})

################################################################
# Optional: Uncomment for optimisation (-O2 and -03 being most common for me)
# target_compile_options(SpinChains PRIVATE -O3)

################################################################

# End of CMakeLists.txt

################################################################
