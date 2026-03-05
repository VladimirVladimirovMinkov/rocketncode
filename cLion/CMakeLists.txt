cmake_minimum_required(VERSION 4.1)
project(cLion LANGUAGES C CXX)

set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(cLion_lib STATIC library.c)

add_executable(cLion
        main.cpp
        gas.cpp
        ambient.cpp
        pipe.cpp)

target_link_libraries(cLion PRIVATE cLion_lib)

if (WIN32)
    set(RELEASE_BINARY_NAME mainWIN32)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(RELEASE_BINARY_NAME mainLINUX)
else ()
    set(RELEASE_BINARY_NAME main)
endif ()

set_target_properties(cLion PROPERTIES
        OUTPUT_NAME_DEBUG cLion
        OUTPUT_NAME_RELEASE ${RELEASE_BINARY_NAME}
        OUTPUT_NAME_RELWITHDEBINFO ${RELEASE_BINARY_NAME}
        OUTPUT_NAME_MINSIZEREL ${RELEASE_BINARY_NAME})

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set_target_properties(cLion PROPERTIES SUFFIX ".exe")
endif ()

if (WIN32 AND CMAKE_BUILD_TYPE STREQUAL "Release")
    add_custom_target(release-dual-artifacts ALL
            COMMAND ${CMAKE_COMMAND} --preset release-linux
            COMMAND ${CMAKE_COMMAND} --build --preset build-release-linux
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            BYPRODUCTS ${CMAKE_SOURCE_DIR}/cmake-build-release-linux/mainLINUX.exe
            COMMENT "Building release dual artifacts: mainWIN32.exe and mainLINUX.exe"
            VERBATIM)

    add_dependencies(release-dual-artifacts cLion)
endif ()
