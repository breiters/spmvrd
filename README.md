This repository includes code to compute private and shared reuse distances of SpMV using the CSR format with a row-wise thread work-sharing. Two different methods are provided with this code: method (A) computes reuse distances considering all memory references, and method (B) computes reuse distances only from references to the x-vector.

In order to compute the reuse distances using method (A), two passes are required. Using method (B), only a single pass is required.

- (1) set `USE_SCALED_REUSE=0` and set `USE_CALC_NOSC_REUSE=0` and compile the code (`$ make -B`). The resulting executable computes reuse distances with method (A) in case no cache partitioning is applied.
- (2) set `USE_SCALED_REUSE=0` and set `USE_CALC_NOSC_REUSE=1`and compile the code (`$ make -B`). The resulting executable computes reuse distances with method (A) in case cache partitioning is applied when the matrix data is assigned to a separate partition.
- (3) set `USE_SCALED_REUSE=1` and compile the code (`$ make -B`). The resulting executable computes reuse distances with method (B).

To run the executable, you need to specify the number of OpenMP threads (set `OMP_NUM_THREADS` accordingly), and provide the input sparse matrix in the Matrix Market format with the command line option `-f`. An example is provided below:

```bash
export OMP_NUM_THREADS=48
t=${OMP_NUM_THREADS}
./spmvrd -f example.mtx -v > reuse-distance-output.csv
```
