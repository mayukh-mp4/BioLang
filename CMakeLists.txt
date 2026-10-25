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
  src/typechecker.cpp src/interpreter.cpp src/motif.cpp)
