cmake_minimum_required(VERSION 3.17)
project(gropt_python LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Detect nanobind — prefer installed package, fall back to find_package
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

# C++ source files (same as gropt/src/CMakeLists.txt minus main.cpp and fft_helper.cpp)
set(GROPT_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gropt/src")
set(GROPT_SOURCES
    ${GROPT_SRC_DIR}/gropt_params.cpp
    ${GROPT_SRC_DIR}/gropt_utils.cpp
    ${GROPT_SRC_DIR}/ils.cpp
    ${GROPT_SRC_DIR}/ils_cg.cpp
    ${GROPT_SRC_DIR}/ils_nlcg.cpp
    ${GROPT_SRC_DIR}/ils_bicgstabl.cpp
    ${GROPT_SRC_DIR}/op_main.cpp
    ${GROPT_SRC_DIR}/op_bvalue.cpp
    ${GROPT_SRC_DIR}/op_gradient.cpp
    ${GROPT_SRC_DIR}/op_identity.cpp
    ${GROPT_SRC_DIR}/op_moment.cpp
    ${GROPT_SRC_DIR}/op_slew.cpp
    ${GROPT_SRC_DIR}/op_safe.cpp
    ${GROPT_SRC_DIR}/op_tv.cpp
    ${GROPT_SRC_DIR}/sdmm_workspace.cpp
    ${GROPT_SRC_DIR}/solver.cpp
    ${GROPT_SRC_DIR}/solver_groptsdmm.cpp
)

nanobind_add_module(gropt_wrapper
    gropt/gropt_wrapper/gropt_wrapper.cpp
    ${GROPT_SOURCES}
)

target_include_directories(gropt_wrapper PRIVATE
    ${GROPT_SRC_DIR}
    ${GROPT_SRC_DIR}/external
)

target_compile_definitions(gropt_wrapper PRIVATE FMT_UNICODE=0)

install(TARGETS gropt_wrapper LIBRARY DESTINATION gropt)

nanobind_add_stub(gropt_wrapper_stub
    MODULE gropt_wrapper
    OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/gropt/gropt_wrapper.pyi"
    PYTHON_PATH "$<TARGET_FILE_DIR:gropt_wrapper>"
    DEPENDS gropt_wrapper
)
