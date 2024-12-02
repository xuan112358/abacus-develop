#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee output

if [[ ! -f output ]] || 
   [[ ! -f OUT.autotest/running_scf.log ]] ||
   [[ ! -f OUT.autotest/data-0-H ]] ||
   [[ ! -f OUT.autotest/data-0-S ]] ||
   [[ ! -f OUT.autotest/data-1-H ]] ||
   [[ ! -f OUT.autotest/data-1-S ]] ||
   [[ ! -f OUT.autotest/data-2-H ]] ||
   [[ ! -f OUT.autotest/data-2-S ]] ||
   [[ ! -f OUT.autotest/data-3-H ]] ||
   [[ ! -f OUT.autotest/data-3-S ]] ||
   [[ ! -f OUT.autotest/data-4-H ]] ||
   [[ ! -f OUT.autotest/data-4-S ]] ||
   [[ ! -f OUT.autotest/data-5-H ]] ||
   [[ ! -f OUT.autotest/data-5-S ]] ||
   [[ ! -f OUT.autotest/data-6-H ]] ||
   [[ ! -f OUT.autotest/data-6-S ]] ||
   [[ ! -f OUT.autotest/data-7-H ]] ||
   [[ ! -f OUT.autotest/data-7-S ]] ||
   [[ ! ( "$(tail -1 OUT.autotest/running_scf.log)" == " Total  Time  :"* ) ]] 
then
	echo "job failed!"
	exit 1
else
	echo "job succeeded!"
	exit 0
fi