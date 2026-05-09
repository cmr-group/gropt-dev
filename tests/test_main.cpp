// Entry point for the gropt C++ test executable.
//
// Each topic-area file exposes a `run_*_tests()` function that returns the
// number of failures. Add new files to tests/CMakeLists.txt and call them here.

#include <cstdio>

int run_op_transpose_tests();

int main() {
    int failures = 0;
    failures += run_op_transpose_tests();

    if (failures > 0) {
        std::fprintf(stderr, "\nFAILED: %d test(s)\n", failures);
        return 1;
    }
    std::printf("\nAll tests passed.\n");
    return 0;
}
