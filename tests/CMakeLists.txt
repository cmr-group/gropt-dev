add_executable(gropt_tests
    test_main.cpp
    test_op_transpose.cpp
    test_add_obj.cpp
    test_safe_axes.cpp
    test_safe_signed.cpp
    test_aux_blocks.cpp
    test_slack_lift.cpp
    test_slew_check.cpp
    test_raster.cpp
    test_diffusion_recipe.cpp
)

target_link_libraries(gropt_tests PRIVATE gropt_core)
target_compile_definitions(gropt_tests PRIVATE
    GROPT_RECIPE_LIBRARY="${PROJECT_SOURCE_DIR}/gropt/diffusion_recipes.json")

add_test(NAME gropt_tests COMMAND gropt_tests)
