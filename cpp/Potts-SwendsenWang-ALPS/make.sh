#!/bin/sh

rm -rf build

mkdir build && cd build

cmake -DENABLE_LOGGING=0 ../src
make potts
