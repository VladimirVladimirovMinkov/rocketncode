cmake_minimum_required(VERSION 4.0)
project(rocketncode)

set(CMAKE_CXX_STANDARD 23)

include_directories(lib/gas)

add_executable(rocketncode
        lib/gas/gas.cpp
        lib/gas/gas.h
        src/main.cpp)
