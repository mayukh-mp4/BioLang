cmake_minimum_required(VERSION 3.16)
project(biolang CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Benchmarks need an optimised build, so default to Release.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_executable(biolang
  src/main.cpp src/lexer.cpp src/parser.cpp
  src/typechecker.cpp src/interpreter.cpp src/motif.cpp src/codegen.cpp)

# Unit tests for lexer, parser and type checker (no external framework).
# Build everything with `ninja`; run with `.\build\unit_tests.exe` or `ctest`.
enable_testing()
add_executable(unit_tests
  tests/unit_tests.cpp src/lexer.cpp src/parser.cpp
  src/typechecker.cpp src/interpreter.cpp src/motif.cpp)
target_include_directories(unit_tests PRIVATE src)
add_test(NAME unit_tests COMMAND unit_tests)
