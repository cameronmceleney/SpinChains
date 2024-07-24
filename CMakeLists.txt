################################################################

# Start of CMakeLists.txt

################################################################
# Latest version of CMake is 3.21 (as of 20 Jan 21)SolversSuperClass.
cmake_minimum_required(VERSION 3.17)
project(SpinChains)  # This project's name is not to be confused with (ChainSpins); an old local redundant project

set(CMAKE_CXX_STANDARD 23)  # No attempt has been made to ensure backwards compatibility from C23
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")
endif()

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
    set(MINGW64_BASE "C:/msys64/mingw64")
    # Windows requires manual installation, or use of UNIX-based Homebrew
    set(EIGEN_DIR "${MINGW64_BASE}/include/eigen3")
    set(FFTW_DIR "${MINGW64_BASE}/fftw3.h")
    set(FFTW_LIB "${MINGW64_BASE}/lib/libfftw3.a")

    # Use TBB (not OpenMP) for parallelisation. Use the newer oneAPI version
    set(TBB_DIR "${MINGW64_BASE}/include/oneapi")  # Removed appended "/oneapi" from the end of the path
    set(TBB_LIB "${MINGW64_BASE}/lib/libtbb12.dll.a" "${MINGW64_BASE}/lib/libtbbmalloc.dll.a")
elseif(APPLE)
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        # Link to Mac-ARM version
        # Default Homebrew installation directories on MacOS
        set(EIGEN_DIR "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
        set(FFTW_BASE "/opt/homebrew/Cellar/fftw/3.3.10_1")
        set(FFTW_DIR "${FFTW_BASE}/include")
        set(FFTW_LIB "${FFTW_BASE}/lib/libfftw3.a")

        # Use TBB (not OpenMP) for parallelisation. Use the newer oneAPI.
        set(TBB_BASE "/opt/homebrew/Cellar/tbb/2021.13.0")
        set(TBB_DIR "${TBB_BASE}/include")
        # Must use dynamic libraries on MacOS (Xcode 15.2), use G++-13, and compile TBB via Homebrew with GCC/G++
        # with "HOMEBREW_CC=gcc-13 HOMEBREW_CXX=g++-13 brew install tbb"
        set(TBB_LIB "${TBB_BASE}/lib/libtbb.dylib" "${TBB_BASE}/lib/libtbbmalloc.dylib")
    else()
        # Link to MacOSX version
        # Default Homebrew installation directories on MacOS
        set(EIGEN_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
        set(FFTW_BASE "/usr/local/Cellar/fftw/3.3.10_1")
        set(FFTW_DIR "${FFTW_BASE}/include")
        set(FFTW_LIB "${FFTW_BASE}/lib/libfftw3.a")

        # Use TBB (not OpenMP) for parallelisation. Use the newer oneAPI.
        set(TBB_BASE "/usr/local/Cellar/tbb/2021.11.0")
        set(TBB_DIR "${TBB_BASE}/include")
        # Must use dynamic libraries on MacOS (Xcode 15.2), use G++-13, and compile TBB via Homebrew with GCC/G++
        # with "HOMEBREW_CC=gcc-13 HOMEBREW_CXX=g++-13 brew install tbb"
        set(TBB_LIB "${TBB_BASE}/lib/libtbb.dylib" "${TBB_BASE}/lib/libtbbmalloc.dylib")
    endif()
endif()

################################################################
# Create lists of all user sources
set(USER_GENERAL_LIBRARIES
        libs/CommonDefinitions.h
        libs/progressbar.hpp
        libs/linspace.cpp libs/linspace.h
        libs/revised_progressbar.hpp
        libs/CommonStructures.h
        libs/HashSpecialisation.h
)

set(USER_OTHER_FILES
        Other/working_on/SpinChainEigenSolverClass.cpp Other/working_on/SpinChainEigenSolverClass.h
        Other/Legacy_Files/NewSolver.cpp Other/Legacy_Files/NewSolver.h
        test/IntelTBBTestIfWorking.cpp test/IntelTBBTestIfWorking.h
        test/mappingSpeedTests.hpp
)

set(CONTAINER_LIBRARIES
        include/SimulationFlags.h
        include/SimulationParameters.h
        include/SimulationStates.h
        src/GlobalVariables.cpp include/GlobalVariables.h
        include/SimulationManager.h
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
        src/ExchangeField.cpp
        include/ExchangeField.h
        src/BiasFields.cpp
        include/BiasFields.h
        src/SolversManager.cpp
        include/SolversManager.h
        libs/CommonStructures.h
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

if(WIN32 OR WIN64)
    include_directories(${EIGEN_DIR} ${FFTW_DIR} ${TBB_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LIB} ${TBB_LIB})
elseif(APPLE)
    include_directories(${EIGEN_DIR} ${FFTW_DIR} ${TBB_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE stdc++ ${FFTW_LIB} ${TBB_LIB})
endif()


################################################################
# Optional: Uncomment for optimisation (-O2 and -O3 being most common for me)
target_compile_options(SpinChains PRIVATE -O3)

################################################################

# End of CMakeLists.txt

################################################################
