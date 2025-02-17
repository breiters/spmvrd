#!/bin/bash

#PJM -N spmvrd-48threads
#PJM -g jh180024o
#PJM -L rscgrp=regular-o
#PJM -L node=1
#PJM --mpi proc=48
#PJM --omp thread=1
#PJM -L elapse=12:00:00
#PJM -o spmvrd-stdout.txt
#PJM -e spmvrd-stderr.txt

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
        ./spmvrd -f $BASE_DIR/mtx/suitesparse/$group/$name -v -c $CONF -o ${RESULT_DIR}/${name}-048threads.csv
    done < ss490matrices.txt
    
    mv overhead048t.csv ${RESULT_DIR}
}

run 0 results/048threads/sector0
run 1 results/048threads/sector1
