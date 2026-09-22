add_executable(gropt_tests
    test_main.cpp
    test_op_transpose.cpp
    test_add_obj.cpp
    test_safe_axes.cpp
)

target_link_libraries(gropt_tests PRIVATE gropt_core)

add_test(NAME gropt_tests COMMAND gropt_tests)
