cmake_minimum_required(VERSION 3.15)
project(EXOSIMS_benchmark)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(benchmark)

add_executable(EXOSIMS_benchmark keplerSTM_benchmarks.cpp)
target_link_libraries(EXOSIMS_benchmark benchmark::benchmark)