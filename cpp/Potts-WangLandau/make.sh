#!/bin/sh 

# set C++ compiler
CXX=g++

rm -rf build
mkdir build && cd build

cmake -DCMAKE_CXX_COMPILER=$CXX ../src
make potts
