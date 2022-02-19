cmake_minimum_required(VERSION 3.0)

project(filter)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# Binaries
add_executable(filter
    ../shared/hmr_args.cpp
    ../shared/hmr_bam.cpp
    ../shared/hmr_bgzf.cpp
    ../shared/hmr_fasta.cpp
    ../shared/hmr_ui.cpp
    src/filter_bam.cpp
    src/filter_enzyme.cpp
    src/filter_fasta.cpp
    src/hmr_args_type.cpp
    src/main.cpp
)
target_link_libraries(filter pthread deflate)
