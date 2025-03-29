cmake_minimum_required(VERSION 3.10)
project(RandomBox LANGUAGES CXX)

# Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directory for headers
include_directories(headers)

# Create shared library for random_placer
add_library(random_placer src/random_placer.cpp)

# Ensure output directories exist inside the build directory
add_custom_target(create_output_dirs ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/output/crystals
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/output/placements
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/output/statistics
)

# Define executables and link with random_placer
add_executable(placement_n_dim_den src/run_D2.cpp)
target_link_libraries(placement_n_dim_den random_placer)

add_executable(brute_density src/run_D2_2.cpp)
target_link_libraries(brute_density random_placer)

add_executable(placement_bool src/run_D2_3.cpp)
target_link_libraries(placement_bool random_placer)
