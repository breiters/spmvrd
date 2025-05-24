#!/bin/bash

if [ -z "$1" ]
    then echo "Need to pass number of threads as first argument"; exit;
fi

export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_WAIT_POLICY=active
# export OMP_DISPLAY_ENV=true
export OMP_NUM_THREADS=$1

THREADS=$(printf %03d $1)
BASE_DIR=/home/ru37geh

run () {
    local CONF=$1
    local RESULT_DIR=$2
    mkdir -p $RESULT_DIR

    while read -r group name
    do
        name=$(echo $name | tr -d [:space:])
        ./spmvrd -f $BASE_DIR/mtx/suitesparse/$group/$name.mtx -v -c $CONF -o ${RESULT_DIR}/${name}-${THREADS}threads.csv
    done < ss490matrices.txt
}

run 0 results/${THREADS}threads/
