cmake_minimum_required(VERSION 3.0)

project(partition)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# Binaries
add_executable(partition
    ../shared/hmr_args.cpp
    ../shared/hmr_cgd_parse.cpp
    ../shared/hmr_ui.cpp
    src/group_conv.cpp
    src/hmr_args_type.cpp
    src/main.cpp
)
target_link_libraries(partition pthread)
