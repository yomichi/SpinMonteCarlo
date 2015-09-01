#!/bin/bash
#PBS -q middle
#PBS -l nodes=4:ppn=12
#PBS -l walltime=24:00:00
#PBS -o /dev/null
#PBS -e /dev/null

export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

source /opt/MateriApps/env.sh
source /opt/MateriApps/alps/alpsvars.sh

prog=../build/ising

rm -f params.*
python init_param.py
mpiexec $prog --mpi params.in.xml > std.out 2> std.err
