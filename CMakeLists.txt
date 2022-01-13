cmake_minimum_required(VERSION 3.20)
project(numberLib)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_subdirectory("ext/bulk")

add_library(numberLib src/fieldArithmetic.cpp)
target_include_directories(numberLib PUBLIC include)

add_executable(nntransform main.cpp)
target_link_libraries(nntransform PUBLIC bulk_thread numberLib)
