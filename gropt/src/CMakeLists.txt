# Standalone C++ build (no Python / nanobind / pixi needed):
#   cmake -S gropt/src -B build-cpp -G Ninja -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-cpp
#   ./build-cpp/gropt
# Optional addons: -DGROPT_PLOTTING=ON (matplot++, needs gnuplot at runtime), -DGROPT_HDF5=ON (HighFive).
cmake_minimum_required(VERSION 3.17)
project(gropt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(${CMAKE_CURRENT_LIST_DIR}/gropt_core.cmake)   # -> gropt_core (shared with the Python build)

# ---- standalone executable -----------------------------------------------------------------------
add_executable(gropt main.cpp demo_diffusion.cpp)
target_link_libraries(gropt PRIVATE gropt_core)

# ---- optional plotting addon (matplot++) ---------------------------------------------------------
option(GROPT_PLOTTING "Build the standalone with matplot++ plotting" OFF)
if(GROPT_PLOTTING)
    add_subdirectory(external/matplotplusplus)
    target_link_libraries(gropt PRIVATE matplot)
    target_compile_definitions(gropt PRIVATE GROPT_PLOTTING)  # gates #ifdef GROPT_PLOTTING in the code
endif()

# ---- optional HDF5 solver-history dumps (HighFive) -----------------------------------------------
option(GROPT_HDF5 "Build the standalone with HighFive/HDF5 solver-history dumps" OFF)
if(GROPT_HDF5)
    add_subdirectory(external/HighFive)
    target_link_libraries(gropt PRIVATE HighFive::HighFive)
    target_compile_definitions(gropt PRIVATE GROPT_HDF5)  # gates #ifdef GROPT_HDF5 in the code
endif()
