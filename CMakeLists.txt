cmake_minimum_required(VERSION 3.17)
project(gropt_python LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Detect nanobind — prefer installed package, fall back to find_package
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

include(${CMAKE_CURRENT_SOURCE_DIR}/gropt/src/gropt_core.cmake)   # -> gropt_core (shared with gropt/src)

nanobind_add_module(gropt_wrapper
    gropt/gropt_wrapper/gropt_wrapper.cpp
)
target_link_libraries(gropt_wrapper PRIVATE gropt_core)

option(GROPT_BUILD_TESTS "Build C++ tests" OFF)
if(GROPT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
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
