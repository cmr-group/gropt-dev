cmake_minimum_required(VERSION 3.10.0)
project(gropt VERSION 0.0.0 LANGUAGES CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external)
add_definitions(-DFMT_UNICODE=0)
add_definitions(-DH5_USE_EIGEN)

add_subdirectory(external/matplotplusplus)
add_subdirectory(external/highfive)


set(SRC_FILES
    fft_helper.cpp
    gropt_params.cpp
    gropt_utils.cpp
    ils.cpp
    ils_cg.cpp
    ils_nlcg.cpp
    ils_bicgstabl.cpp
    op_bvalue.cpp
    op_gradient.cpp
    op_identity.cpp
    op_main.cpp
    op_moment.cpp
    op_slew.cpp
    op_safe.cpp
    op_tv.cpp
    solver.cpp
    solver_groptsdmm.cpp
)

add_executable(gropt main.cpp ${SRC_FILES})
target_link_libraries(gropt PUBLIC matplot HighFive::HighFive)  

