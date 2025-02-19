#!/bin/bash

#PJM -N spmvrd-01threads
#PJM -g jh180024o
#PJM -L rscgrp=regular-o
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=1
#PJM -L elapse=12:00:00
#PJM -o spmvrd-48threads-stdout.txt
#PJM -e spmvrd-48threads-stderr.txt

export LC_ALL=C
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_WAIT_POLICY=active
export OMP_DISPLAY_ENV=true

export XOS_MMM_L_HPAGE_TYPE=hugetlbfs
export XOS_MMM_L_PAGING_POLICY=demand:demand:demand

BASE_DIR=

run () {
    local CONF=$1
    local RESULT_DIR=$2
    mkdir -p $RESULT_DIR

    while read -r group name
    do
        ./spmvrd -f $BASE_DIR/mtx/suitesparse/$group/$name -v -c $CONF -o ${RESULT_DIR}/${name}-001threads.csv
    done < ss490matrices.txt

    mv overhead001t.csv ${RESULT_DIR}
}

run 0 results/001threads/sector0
run 1 results/001threads/sector1
