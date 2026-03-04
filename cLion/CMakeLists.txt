cmake_minimum_required(VERSION 3.20)
project(cLion CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(cLion
        main.cpp
        gas.cpp
        ambient.cpp
        pipe.cpp
        pipe.h)
