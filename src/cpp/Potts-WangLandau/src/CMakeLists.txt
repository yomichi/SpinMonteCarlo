cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(wanglandau)

find_package(Boost 1.55 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_executable(potts potts.cpp wanglandau.cpp)
target_link_libraries(potts ${Boost_LIBRARIES})
