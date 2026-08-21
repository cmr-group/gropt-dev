# Single source of truth for the gropt core .cpp list. Included by BOTH the top-level Python build
# (CMakeLists.txt) and the standalone C++ build (gropt/src/CMakeLists.txt) so the list never drifts.
# Add a new core source here once; both builds pick it up. (Explicit list on purpose -- no file(GLOB).)
# main.cpp and demo_diffusion.cpp are NOT here: they belong only to the standalone executable.
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
