cmake_minimum_required(VERSION 3.0)

project(dump)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# Binaries
add_executable(dump
    ../shared/hmr_args.cpp
    ../shared/hmr_cgd_parse.cpp
    ../shared/hmr_cluster.cpp
    ../shared/hmr_reads_parse.cpp
    ../shared/hmr_tour.cpp
    ../shared/hmr_ui.cpp
    src/hmr_args_type.cpp
    src/main.cpp
    src/optimize_direction.cpp
    src/optimize_order.cpp
)
target_link_libraries(dump pthread)
