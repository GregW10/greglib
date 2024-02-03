cmake_minimum_required(VERSION 3.21)
project(gregml)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(main main.cpp)
add_executable(test test.cpp)
