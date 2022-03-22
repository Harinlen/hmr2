cmake_minimum_required(VERSION 3.0)

project(correct)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# zlib
find_package(ZLIB)

# Binaries
add_executable(correct
    ../shared/hmr_args.cpp
    ../shared/hmr_bam.cpp
    ../shared/hmr_bgzf.cpp
    ../shared/hmr_bin_file.cpp
    ../shared/hmr_bin_queue.cpp
    ../shared/hmr_fasta.cpp
    ../shared/hmr_gz.cpp
    ../shared/hmr_mapping.cpp
    ../shared/hmr_path.cpp
    ../shared/hmr_text_file.cpp
    ../shared/hmr_ui.cpp
    src/args_correct.cpp
    src/contig_correct.cpp
    src/mapping_correct.cpp
    src/main.cpp
    src/mismatch_correct.cpp
)
target_link_libraries(correct pthread ZLIB::ZLIB)
