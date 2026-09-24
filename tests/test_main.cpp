// Entry point for the gropt C++ test executable.
//
// Each topic-area file exposes a `run_*_tests()` function that returns the
// number of failures. Add new files to tests/CMakeLists.txt and call them here.

#include <cstdio>

int run_op_transpose_tests();
int run_add_obj_tests();
int run_safe_axes_tests();
int run_safe_signed_tests();
int run_aux_block_tests();
int run_slack_lift_tests();
int run_slew_check_tests();

int main() {
    int failures = 0;
    failures += run_op_transpose_tests();
    failures += run_add_obj_tests();
    failures += run_safe_axes_tests();
    failures += run_safe_signed_tests();
    failures += run_aux_block_tests();
    failures += run_slack_lift_tests();
    failures += run_slew_check_tests();

    if (failures > 0) {
        std::fprintf(stderr, "\nFAILED: %d test(s)\n", failures);
        return 1;
    }
    std::printf("\nAll tests passed.\n");
    return 0;
}
