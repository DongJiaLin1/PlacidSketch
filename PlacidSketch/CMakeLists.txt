cmake_minimum_required(VERSION 3.15)
project(PlacidSketch CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_executable(placidsketch ${SOURCES})

target_include_directories(placidsketch PRIVATE ${CMAKE_SOURCE_DIR})
