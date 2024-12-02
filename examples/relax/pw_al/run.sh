#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee relax.output

if [[ ! -f relax.output ]] || 
   [[ ! -f OUT.ABACUS/running_cell-relax.log ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_cell-relax.log)" == " Total  Time  :"* ) ]] 
then
	echo "job failed!"
	exit 1
else
	echo "job succeeded!"
	exit 0
fi