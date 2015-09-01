#!/bin/sh

rm -f params.*
python init_param.py
mpiexec ../build/potts
