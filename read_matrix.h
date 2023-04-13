#pragma once

#include <assert.h>
#include <inttypes.h>
// #include <likwid-marker.h>
// #include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef __cplusplus
#    define restrict __restrict
#endif

template <typename VT, typename RT, typename CT>
struct matrix_csr {
    VT *restrict val;
    RT *restrict row_ptr;
    CT *restrict col_idx;
    uint64_t nrow;
    uint64_t ncol;
    uint64_t nnz;

    static matrix_csr read_matrix(const char *path)
    {
        char buf[BUFSIZ];
        bool first = true;
        // struct
        matrix_csr mat;
        uint64_t   i           = 0;
        uint64_t   current_row = 0;

        FILE *f = fopen(path, "r");
        while (fgets(buf, BUFSIZ - 1, f)) {
            if (buf[0] == '%')
                continue;
            if (first) {
                first = false;
                sscanf(buf, "%lu %lu %lu", &mat.nrow, &mat.ncol, &mat.nnz);
                mat.row_ptr = (RT *restrict)malloc(sizeof(RT) * (mat.nrow + 1));
                mat.col_idx = (CT *restrict)malloc(sizeof(CT) * mat.nnz);
                mat.val     = (VT *restrict)malloc(sizeof(VT) * mat.nnz);
                // printf("%s", buf);
            } else {
                uint64_t col;
                uint64_t row;
                double   val;
                if (sscanf(buf, "%lu %lu %lf", &row, &col, &val) != 3) {
                    fprintf(stderr, "%s", buf);
                    assert(false);
                }
                if (row != current_row) {
                    mat.row_ptr[current_row] = i;
                    current_row              = row;
                }
                mat.col_idx[i] = col - 1;
                mat.val[i]     = val;
                i++;
            }
            mat.row_ptr[mat.nrow] = mat.nnz;
        }
        // printf("last idx: %lu, nnz: %lu\n", i, mat.nnz);
        fclose(f);
        return mat;
    }
};
