cmake_minimum_required(VERSION 4.0)
project(PSSketch)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mbmi -Wall -Wextra")

add_executable(PSSketch PSSketch.cpp)
