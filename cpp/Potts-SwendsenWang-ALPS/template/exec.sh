#!/bin/sh

rm -f params.*
python init_param.py
../build/potts params.in.xml
