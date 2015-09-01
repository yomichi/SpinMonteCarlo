cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(exchange)

find_package(Boost 1.55 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

find_package(MPI)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_executable(potts potts_main.cpp)
target_link_libraries(potts ${Boost_LIBRARIES})
