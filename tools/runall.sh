#!/bin/bash

export OMP_PLACES=cores OMP_PROC_BIND=close

for f in ~/mtx/suitesparse/*/*.mtx; do
	t=1
	OMP_NUM_THREADS=$t ./../spmvrd -f $f -v > $f-00${t}threads.csv
        t=48
        OMP_NUM_THREADS=$t ./../spmvrd -f $f -v > $f-0${t}threads.csv
done

