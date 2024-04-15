################################################################

# Start of CMakeLists.txt

################################################################
# Latest version of CMake is 3.21 (as of 20 Jan 21)SolversSuperClass.
cmake_minimum_required(VERSION 3.17)
project(SpinChains)  # This project's name is not to be confused with (ChainSpins); an old local redundant project

set(CMAKE_CXX_STANDARD 17)  # No attempt has been made to ensure backwards compatibility from C23
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_language(OBJCXX)

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")
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

# Path to your metal-cpp directory
set(METAL_CPP_DIR ${CMAKE_SOURCE_DIR}/metal-cpp)
if(APPLE)
    find_library(METAL Metal)
    find_library(FOUNDATION Foundation)
    find_library(QUARTZ_CORE QuartzCore)
    find_library(APP_KIT AppKit)
endif()

################################################################

# Add the include directories for metal-cpp and its dependencies
# (Metal-CPP, QuartzCore-CPP, and Foundation-CPP headers)


# OPTIONAL: Find the system Metal and Foundation frameworks
# find_library(METAL_FRAMEWORK Metal)
# find_library(FOUNDATION_FRAMEWORK Foundation)
# find_library(QUARTZCORE_FRAMEWORK QuartzCore)

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
        set(TBB_BASE "/opt/homebrew/Cellar/tbb/2021.11.0")
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
set(METAL_GENERAL_LIBRARY
    ${CMAKE_SOURCE_DIR}/libs/metal_config.h
        libs/LocalMetalResources.h
)

set(METAL_CHECKS_EXAMPLES
        #checks_examples/check_for_metal_device.cpp
        #checks_examples/check_for_metal_device.h
        checks_examples/CommonTools.cpp
        checks_examples/CommonTools.h
)

set(USER_GENERAL_LIBRARIES
        libs/CommonDefinitions.h
        libs/progressbar.hpp
        libs/linspace.cpp
        libs/linspace.h
        libs/revised_progressbar.hpp
        libs/CommonStructures.h
        libs/HashSpecialisation.h
        libs/CustomTimer.hpp

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

################################################################
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")

# Assuming your .metal file is located at src/my_shader.metal
set(METAL_FILE "${CMAKE_SOURCE_DIR}/src/coffeewitharch.metal")

# Mark the .metal file as a header-only file to include it in the project for IDEs
set_source_files_properties(${METAL_FILE} PROPERTIES HEADER_FILE_ONLY TRUE)

################################################################
# Links all user sources to the executable
add_executable(${PROJECT_NAME} src/main.cpp
        ${METAL_GENERAL_LIBRARY}
        ${METAL_CHECKS_EXAMPLES}
        ${USER_GENERAL_LIBRARIES}
        ${USER_OTHER_FILES}
        ${CONTAINER_LIBRARIES}
        ${SOLVERS_CLASS}
        ${SOLVER_CLASS_COMPONENTS}
        ${METAL_FILE}
)



target_include_directories(${PROJECT_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/libs
        ${METAL_CPP_DIR}
        ${METAL_CPP_DIR}/AppKit
        ${METAL_CPP_DIR}/Foundation
        ${METAL_CPP_DIR}/Metal
        ${METAL_CPP_DIR}/MetalFX
        ${METAL_CPP_DIR}/MetalKit
        ${METAL_CPP_DIR}/QuartzCore)

if (APPLE)
    target_include_directories(${PROJECT_NAME} PRIVATE ${EIGEN_DIR})
    target_include_directories(${PROJECT_NAME} PRIVATE ${FFTW_DIR})
    target_include_directories(${PROJECT_NAME} PRIVATE ${TBB_DIR})
endif ()
################################################################
# NEW. for create lists of all sources
# Print the CMAKE_CURRENT_BINARY_DIR variable
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")

# Path to Metal shader source file
set(METAL_SHADER_SRC ${CMAKE_SOURCE_DIR}/src/coffeewitharch.metal)

# Output path for the compiled Metal shader (AIR file)
set(METAL_SHADER_AIR ${CMAKE_CURRENT_BINARY_DIR}/coffeewitharch.air)

# Output path for the Metal library (.metallib file)
set(METAL_SHADER_METALLIB ${CMAKE_CURRENT_BINARY_DIR}/coffeewitharch.metallib)
add_definitions(-DMETAL_SHADER_METALLIB_PATH="${METAL_SHADER_METALLIB}")
target_compile_definitions(${PROJECT_NAME} PRIVATE
        NS_PRIVATE_IMPLEMENTATION
        CA_PRIVATE_IMPLEMENTATION
        MTL_PRIVATE_IMPLEMENTATION
        MTK_PRIVATE_IMPLEMENTATION
)


################################################################
# NEW. Custom commands

# Compile .metal shader to .air
add_custom_command(
    OUTPUT ${METAL_SHADER_AIR}
    COMMAND xcrun -sdk macosx metal -c ${METAL_SHADER_SRC} -o ${METAL_SHADER_AIR}
    DEPENDS ${METAL_SHADER_SRC}
    COMMENT "Compiling ${METAL_SHADER_SRC} to AIR"
)

# Compile .air to .metallib
add_custom_command(
    OUTPUT ${METAL_SHADER_METALLIB}
    COMMAND xcrun -sdk macosx metallib ${METAL_SHADER_AIR} -o ${METAL_SHADER_METALLIB}
    DEPENDS ${METAL_SHADER_AIR}
    COMMENT "Compiling ${METAL_SHADER_AIR} to Metal Library"
)

################################################################
# NEW. Custom target to ensure the Metal shader is compiled
add_custom_target(
    CompileMetalShader ALL
    DEPENDS ${METAL_SHADER_METALLIB}
)

################################################################
# NEW. Specify that the executable depends on the CompileMetalShader target
add_dependencies(${PROJECT_NAME} CompileMetalShader)

################################################################
# Include and link directories for third-party libraries
# Both FFTW and TBB require linking to static libraries; Eigen does not.

if(WIN32 OR WIN64)
    include_directories(${EIGEN_DIR} ${FFTW_DIR} ${TBB_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LIB} ${TBB_LIB})
elseif(APPLE)
    set(USER_FLAGS c++ c++abi)
    target_link_libraries(${PROJECT_NAME} PRIVATE
        ${METAL}
        ${FOUNDATION}
        ${QUARTZ_CORE}
        ${APP_KIT}
        ${FFTW_LIB} ${TBB_LIB}
        ${USER_FLAGS})
endif()


################################################################
# Optional: Uncomment for optimisation (-O2 and -O3 being most common for me)
# target_compile_options(SpinChains PRIVATE -O2)

################################################################

# End of CMakeLists.txt

################################################################
