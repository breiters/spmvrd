#pragma once

#include <cinttypes>

#ifdef __cplusplus
#    define restrict __restrict
#endif

template <typename VT, typename RT, typename CT>
struct matrix_csr {

    ~matrix_csr()
    {
        free(val);
        free(row_ptr);
        free(col_idx);
    }

    VT *restrict val{nullptr};
    RT *restrict row_ptr{nullptr};
    CT *restrict col_idx{nullptr};

    uint64_t nrow{0u};
    uint64_t ncol{0u};
    uint64_t nnz{0u};

    const char *name{nullptr};
    bool        symmetric{false};

    struct coo_entry {
        VT val;
        RT row;
        CT col;
    };

    static void read_matrix_coo(FILE *f, matrix_csr &mat, bool pattern, bool symmetric) {}

    static void matrix_csr_from_coo(matrix_csr &mat, std::vector<coo_entry> &v) {}

    static matrix_csr read_matrix(const char *path) { return {}; };
};
