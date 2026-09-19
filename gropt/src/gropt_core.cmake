# The gropt_core static library, shared by the Python build (CMakeLists.txt) and the standalone build
# (gropt/src/CMakeLists.txt). main.cpp and demo_diffusion.cpp belong to the standalone executable only.
set(GROPT_CORE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/equilibrate.cpp
    ${CMAKE_CURRENT_LIST_DIR}/fft_tools.cpp
    ${CMAKE_CURRENT_LIST_DIR}/gropt_params.cpp
    ${CMAKE_CURRENT_LIST_DIR}/gropt_utils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ils_cg.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ils_nlcg.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ils_bicgstabl.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_bvalue.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_concomitant.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_gradient.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_identity.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_moment.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_slew.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_safe.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_eddy.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_tv.cpp
    ${CMAKE_CURRENT_LIST_DIR}/op_diffbasin.cpp
    ${CMAKE_CURRENT_LIST_DIR}/solver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/solver_groptsdmm.cpp
    ${CMAKE_CURRENT_LIST_DIR}/step_monitor.cpp
    ${CMAKE_CURRENT_LIST_DIR}/solver_osqp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/workspace_osqp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/workspace_sdmm.cpp
    ${CMAKE_CURRENT_LIST_DIR}/workspace_solver.cpp
)

add_library(gropt_core STATIC ${GROPT_CORE_SOURCES})
target_compile_features(gropt_core PUBLIC cxx_std_17)
set_target_properties(gropt_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(gropt_core PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
    ${CMAKE_CURRENT_LIST_DIR}/external
)
target_compile_definitions(gropt_core PUBLIC FMT_UNICODE=0)
target_compile_options(gropt_core PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/MP>) # MSVC: compile sources in parallel
# MSVC: C++ exception unwinding. Normally in CMake's default flags, but absent if CMAKE_CXX_FLAGS is overridden.
target_compile_options(gropt_core PUBLIC $<$<CXX_COMPILER_ID:MSVC>:/EHsc>)

# ON keeps Eigen assertions and NaN-initializes matrices (testing); OFF for distribution builds.
option(GROPT_EIGEN_ASSERTIONS "Enable Eigen runtime assertions (for testing)" OFF)
if(GROPT_EIGEN_ASSERTIONS)
    target_compile_options(gropt_core PRIVATE $<IF:$<CXX_COMPILER_ID:MSVC>,/UNDEBUG,-UNDEBUG>)
    target_compile_definitions(gropt_core PRIVATE EIGEN_INITIALIZE_MATRICES_BY_NAN)
else()
    target_compile_definitions(gropt_core PRIVATE EIGEN_NO_DEBUG)
endif()
