set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

set(ZIG_TARGET "x86_64-linux-gnu" CACHE STRING "Zig target triple for Linux cross-compilation")

if(DEFINED ENV{ZIG_EXE})
    set(ZIG_EXE "$ENV{ZIG_EXE}")
else()
    find_program(ZIG_EXE
            NAMES zig
            PATHS
            "$ENV{LOCALAPPDATA}/Programs/Zig"
            "$ENV{ProgramFiles}/zig"
            "$ENV{USERPROFILE}/scoop/shims"
            "$ENV{USERPROFILE}/scoop/apps/zig/current")

    if(NOT ZIG_EXE)
        file(GLOB_RECURSE ZIG_CANDIDATES
                "$ENV{USERPROFILE}/Downloads/zig*/zig.exe"
                "$ENV{USERPROFILE}/Downloads/zig*/*/zig.exe")
        list(LENGTH ZIG_CANDIDATES ZIG_CANDIDATE_COUNT)
        if(ZIG_CANDIDATE_COUNT GREATER 0)
            list(GET ZIG_CANDIDATES 0 ZIG_EXE)
        endif()
    endif()
endif()

if(NOT ZIG_EXE)
    message(FATAL_ERROR
            "Zig compiler not found. Install Zig and ensure 'zig' is in PATH, "
            "or set ZIG_EXE to the full path to zig.exe.")
endif()

set(CMAKE_C_COMPILER "${ZIG_EXE}")
set(CMAKE_C_COMPILER_ARG1 cc)
set(CMAKE_C_COMPILER_TARGET "${ZIG_TARGET}")

set(CMAKE_CXX_COMPILER "${ZIG_EXE}")
set(CMAKE_CXX_COMPILER_ARG1 c++)
set(CMAKE_CXX_COMPILER_TARGET "${ZIG_TARGET}")

set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
