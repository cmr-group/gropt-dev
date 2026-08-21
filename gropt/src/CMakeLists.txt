# Standalone C++ build (no Python / nanobind / pixi needed):
#   cmake -S gropt/src -B build-cpp -G Ninja -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-cpp
#   ./build-cpp/gropt
# Plotting/HDF5 addons are optional -- turn them off for a minimal build (no matplot++, HighFive, HDF5,
# gnuplot):  -DGROPT_PLOTTING=OFF
cmake_minimum_required(VERSION 3.17)
project(gropt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- core library (shared source list; NO plotting deps) -----------------------------------------
include(${CMAKE_CURRENT_LIST_DIR}/sources.cmake)   # -> GROPT_CORE_SOURCES

add_library(gropt_core STATIC ${GROPT_CORE_SOURCES})
target_include_directories(gropt_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/external
)
target_compile_definitions(gropt_core PUBLIC FMT_UNICODE=0)
target_compile_definitions(gropt_core PRIVATE POCKETFFT_NO_MULTITHREADING POCKETFFT_CACHE_SIZE=8)
# /EHsc: the standalone throws on bad input (recipe JSON, unknown modes); MSVC needs it for
# proper stack unwinding (GCC/Clang enable exceptions by default).
target_compile_options(gropt_core PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/MP /EHsc>)

option(GROPT_EIGEN_ASSERTIONS "Enable Eigen runtime assertions (for testing)" OFF)
if(GROPT_EIGEN_ASSERTIONS)
    target_compile_options(gropt_core PRIVATE $<IF:$<CXX_COMPILER_ID:MSVC>,/UNDEBUG,-UNDEBUG>)
    target_compile_definitions(gropt_core PRIVATE EIGEN_INITIALIZE_MATRICES_BY_NAN)
else()
    target_compile_definitions(gropt_core PRIVATE EIGEN_NO_DEBUG)
endif()

# ---- standalone executable -----------------------------------------------------------------------
add_executable(gropt main.cpp demo_diffusion.cpp)
target_link_libraries(gropt PRIVATE gropt_core)
target_compile_options(gropt PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/MP /EHsc>)

# ---- optional plotting addon (matplot++) ---------------------------------------------------------
option(GROPT_PLOTTING "Build the standalone with matplot++ plotting" ON)
if(GROPT_PLOTTING)
    add_subdirectory(external/matplotplusplus)
    target_link_libraries(gropt PRIVATE matplot)
    target_compile_definitions(gropt PRIVATE GROPT_PLOTTING)  # gates #ifdef GROPT_PLOTTING in the code
endif()

# ---- optional HDF5 solver-history dumps (HighFive) -----------------------------------------------
option(GROPT_HDF5 "Build the standalone with HighFive/HDF5 solver-history dumps" OFF)
if(GROPT_HDF5)
    add_subdirectory(external/highfive)
    target_link_libraries(gropt PRIVATE HighFive::HighFive)
    target_compile_definitions(gropt PRIVATE GROPT_HDF5)  # gates #ifdef GROPT_HDF5 in the code
endif()
