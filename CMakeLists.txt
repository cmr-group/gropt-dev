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
    ${GROPT_SRC_DIR}/equilibrate.cpp
    ${GROPT_SRC_DIR}/gropt_params.cpp
    ${GROPT_SRC_DIR}/gropt_utils.cpp
    ${GROPT_SRC_DIR}/ils.cpp
    ${GROPT_SRC_DIR}/ils_cg.cpp
    ${GROPT_SRC_DIR}/ils_nlcg.cpp
    ${GROPT_SRC_DIR}/ils_bicgstabl.cpp
    ${GROPT_SRC_DIR}/op_main.cpp
    ${GROPT_SRC_DIR}/op_bvalue.cpp
    ${GROPT_SRC_DIR}/op_concomitant.cpp
    ${GROPT_SRC_DIR}/op_gradient.cpp
    ${GROPT_SRC_DIR}/op_identity.cpp
    ${GROPT_SRC_DIR}/op_moment.cpp
    ${GROPT_SRC_DIR}/op_slew.cpp
    ${GROPT_SRC_DIR}/op_safe.cpp
    ${GROPT_SRC_DIR}/op_eddy.cpp
    ${GROPT_SRC_DIR}/op_tv.cpp
    ${GROPT_SRC_DIR}/solver.cpp
    ${GROPT_SRC_DIR}/solver_groptsdmm.cpp
    ${GROPT_SRC_DIR}/solver_osqp.cpp
    ${GROPT_SRC_DIR}/workspace_osqp.cpp
    ${GROPT_SRC_DIR}/workspace_sdmm.cpp
    ${GROPT_SRC_DIR}/workspace_solver.cpp
)

add_library(gropt_core STATIC ${GROPT_SOURCES})
set_target_properties(gropt_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(gropt_core PUBLIC
    ${GROPT_SRC_DIR}
    ${GROPT_SRC_DIR}/external
)
target_compile_definitions(gropt_core PUBLIC FMT_UNICODE=0)

nanobind_add_module(gropt_wrapper
    gropt/gropt_wrapper/gropt_wrapper.cpp
)
target_link_libraries(gropt_wrapper PRIVATE gropt_core)

# Eigen assertion control: ON = keep assertions + NaN-init for testing,
#                          OFF = disable assertions for distribution wheels
option(GROPT_EIGEN_ASSERTIONS "Enable Eigen runtime assertions (for testing)" OFF)
if(GROPT_EIGEN_ASSERTIONS)
    target_compile_options(gropt_core PRIVATE
        $<IF:$<CXX_COMPILER_ID:MSVC>,/UNDEBUG,-UNDEBUG>
    )
    target_compile_definitions(gropt_core PRIVATE EIGEN_INITIALIZE_MATRICES_BY_NAN)
else()
    target_compile_definitions(gropt_core PRIVATE EIGEN_NO_DEBUG)
endif()

install(TARGETS gropt_wrapper LIBRARY DESTINATION gropt)

# Stub generation requires importing the built module, so it only works for native builds.
# Skip when cross-compiling on macOS (e.g. building x86_64 on an arm64 runner).
set(_generate_stubs ON)
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
    execute_process(
        COMMAND uname -m
        OUTPUT_VARIABLE _host_arch
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT CMAKE_OSX_ARCHITECTURES STREQUAL _host_arch)
        set(_generate_stubs OFF)
    endif()
endif()

if(_generate_stubs)
    nanobind_add_stub(gropt_wrapper_stub
        MODULE gropt_wrapper
        OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/gropt/gropt_wrapper.pyi"
        PYTHON_PATH "$<TARGET_FILE_DIR:gropt_wrapper>"
        DEPENDS gropt_wrapper
    )
endif()
