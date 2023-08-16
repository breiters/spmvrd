#include "cache.h"
#include "matrix_csr.h"
#include "read_matrix.h"

#include <omp.h>

#include <array>
#include <string>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#define USE_SCALED_REUSE    0
#define USE_CALC_NOSC_REUSE 0

using rowptr_t = int64_t;
using colidx_t = int;
using val_t    = double;
using x_t      = val_t;
using y_t      = val_t;

void set_buckets_a64fx(const auto &matrix)
{
    (void)matrix;
    // required buckets for a64fx:
    // 4-way L1d 64KiB => 4 Buckets with distance 64KiB / 4
    // 16-way L2 8MiB => 16 Buckets with distance 8MiB / 16
    //
    int KiB = 1024;
    int MiB = 1024 * KiB;

    int L1ways = 4;
    int L2ways = 16;

    int L1d_capacity_per_way = 64 * KiB / 4;
    int L2_capacity_per_way  = 8 * MiB / 16;

    Bucket::min_dists.push_back(0);
    
    for (int i = 0; i < L1ways; i++)
        Bucket::min_dists.push_back(L1d_capacity_per_way * (i + 1) / MEMBLOCKLEN);

    for (int i = 0; i < L2ways; i++)
        Bucket::min_dists.push_back(L2_capacity_per_way * (i + 1) / MEMBLOCKLEN);

    // bucket for cold misses (infinite reuse distance)
    Bucket::min_dists.push_back(Bucket::INF_DIST);

    // sort buckets in ascending order
    std::sort(Bucket::min_dists.begin(), Bucket::min_dists.end());
}

void set_buckets_a64fx_scaled(const auto &matrix)
{
    // required buckets for a64fx:
    // 4-way L1d 64KiB => 4 Buckets with distance 64KiB / 4
    // 16-way L2 8MiB => 16 Buckets with distance 8MiB / 16
    //
    int KiB = 1024;
    int MiB = 1024 * KiB;

    int L1ways = 4;
    int L2ways = 16;

    int L1d_capacity_per_way = 64 * KiB / 4;
    int L2_capacity_per_way  = 8 * MiB / 16;

    Bucket::min_dists.push_back(0);

    double min_scale_nosc =
        (double)sizeof(val_t) / ((double)(sizeof(rowptr_t) + sizeof(y_t)) * (double)matrix.nrow / matrix.nnz +
                                 sizeof(val_t) + sizeof(colidx_t) + sizeof(x_t));

    double min_scale_sc = (double)sizeof(val_t) /
                          ((double)(sizeof(rowptr_t) + sizeof(y_t)) * (double)matrix.nrow / matrix.nnz + sizeof(x_t));

    for (int i = 0; i < L1ways; i++) {
        Bucket::min_type min = min_scale_sc * L1d_capacity_per_way * (i + 1) / MEMBLOCKLEN;
        Bucket::min_dists.push_back(min);
    }

    for (int i = 0; i < L2ways; i++) {
        Bucket::min_type min = min_scale_sc * L2_capacity_per_way * (i + 1) / MEMBLOCKLEN;
        Bucket::min_dists.push_back(min);
    }

    Bucket::min_dists.push_back(min_scale_nosc * 64 * KiB / MEMBLOCKLEN);
    Bucket::min_dists.push_back(min_scale_nosc * 8 * MiB / MEMBLOCKLEN);

    // bucket for cold misses (infinite reuse distance)
    Bucket::min_dists.push_back(Bucket::INF_DIST);

    // remove duplicated buckets
    auto &vec = Bucket::min_dists;
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

    // sort buckets in ascending order
    std::sort(Bucket::min_dists.begin(), Bucket::min_dists.end());

    for (auto b : Bucket::min_dists) {
        eprintf("Bucket order: %lu\n", b);
    }
}

#if 0
void set_buckets()
{
    unsigned long KiB = 1024;

    Bucket::min_dists.push_back(0);

    for (int i = 0; i < 20; i++) {
        Bucket::min_dists.push_back(KiB / MEMBLOCKLEN);
        KiB *= 1.5;
    }

    // bucket for cold misses (infinite reuse distance)
    Bucket::min_dists.push_back(Bucket::INF_DIST);

    // sort buckets in ascending order
    std::sort(Bucket::min_dists.begin(), Bucket::min_dists.end());
}

void make_numa(matrix_csr<double, uint32_t, uint32_t> &matrix)
{
    // printf("rowptr 0: %d\n", matrix.row_ptr[0]);
    uint32_t *colidx_ = (uint32_t *)malloc(sizeof(uint32_t) * matrix.nnz);
    uint32_t *rowptr_ = (uint32_t *)malloc(sizeof(uint32_t) * (matrix.nrow + 1));
#    pragma omp parallel
#    pragma omp for
    for (uint r = 0; r < matrix.nrow + 1; ++r) {
        rowptr_[r] = matrix.row_ptr[r];
        for (uint i = matrix.row_ptr[r]; i < matrix.row_ptr[r + 1]; ++i) {
            colidx_[i] = matrix.col_idx[i];
        }
    }
    rowptr_[matrix.nrow] = matrix.row_ptr[matrix.nrow];

    free(matrix.row_ptr);
    free(matrix.col_idx);
    matrix.row_ptr = rowptr_;
    matrix.col_idx = colidx_;
}
#endif

void reuse_scaled(int tid, PrivateCache &pc, SharedCache &sc, const auto &matrix)
{
#pragma omp for
    for (uint r = 0; r < matrix.nrow; ++r) {
        // fprintf(stderr, "row: %d\n", r);
        for (rowptr_t i = matrix.row_ptr[r]; i < matrix.row_ptr[r + 1]; ++i) {
            auto cl = cline<val_t, MEMBLOCKLEN>(matrix.col_idx[i]);
            // printf("row: %u, coldix: %u, cline: %lu val: %f i: %u\n", r,
            // matrix.col_idx[i], cl, matrix.val[i], i);
            pc.handle_cline(cl);
            sc.handle_cline_shared(tid, cl);
        }
    }
}

void reuse_calc(int tid, PrivateCache &pc, SharedCache &sc, const auto &matrix)
{
    //           x[0]...x[ncol] <padding> rowptr[0] ... rowptr[nrow] ...
    //
    // cacheline(x[0]) = 0 ... cacheline(ncol) = ncol * sizeof(val_t) / MEMBLOCKLEN ...
    auto cl_x_end     = cline<val_t, MEMBLOCKLEN>(matrix.ncol);
    auto cl_row_start = cl_x_end + 1;
    auto cl_y_start   = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(matrix.nrow + 1) + 1;
    auto cl_a_start   = cl_y_start + cline<val_t, MEMBLOCKLEN>(matrix.nrow) + 1;
    auto cl_col_start = cl_a_start + cline<val_t, MEMBLOCKLEN>(matrix.nnz) + 1;

    // row ptr[0]
    auto cl_row = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(0);
    pc.handle_cline(cl_row);
    sc.handle_cline_shared(tid, cl_row);

#pragma omp for
    for (uint r = 0; r < matrix.nrow; ++r) {
        // fprintf(stderr, "row: %d\n", r);
        // rowptr[r + 1]
        auto cl_row_plus1 = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(r + 1);
        pc.handle_cline(cl_row_plus1);

        // y[r]
        auto cl_y = cl_y_start + cline<val_t, MEMBLOCKLEN>(r);
        pc.handle_cline(cl_y);

        sc.handle_clines_shared(tid, cl_row_plus1, cl_y);

        for (rowptr_t i = matrix.row_ptr[r]; i < matrix.row_ptr[r + 1]; ++i) {
#if USE_CALC_NOSC_REUSE
            // a[i]
            auto cl_a = cl_a_start + cline<val_t, MEMBLOCKLEN>(i);
            pc.handle_cline(cl_a);
            // col_idx[i]
            auto cl_col = cl_col_start + cline<colidx_t, MEMBLOCKLEN>(i);
            pc.handle_cline(cl_col);
#endif /* USE_CALC_NOSC_REUSE */
            // x[col_idx[i]]
            auto cl_x = cline<val_t, MEMBLOCKLEN>(matrix.col_idx[i]);
            pc.handle_cline(cl_x);
#if USE_CALC_NOSC_REUSE
            sc.handle_clines_shared(tid, cl_a, cl_col, cl_x);
#else
            sc.handle_cline_shared(tid, cl_x);
#endif /* USE_CALC_NOSC_REUSE */
        }
    }
}

int main(int argc, char *argv[])
{
    char *matrix_path = nullptr;
    FILE *csv_file    = stdout;
    bool  verbose     = false;

    int opt;
    while ((opt = getopt(argc, argv, "f:o:v")) != -1) {
        switch (opt) {
        case 'f':
            matrix_path = optarg;
            break;
        case 'o':
            csv_file = fopen(optarg, "w+");
            if (!csv_file) {
                perror("fopen (csv file)");
                exit(EXIT_FAILURE);
            }
            break;
        case 'v':
            verbose = true;
            break;
        default: /* '?' */
usage:
            fprintf(stderr, "Usage: %s -f <matrix file>  [-o csv file] [-v]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (!matrix_path)
        goto usage;

    char overhead_csv_path[1024];
    snprintf(overhead_csv_path, 1024, "overhead%03dt.csv", omp_get_max_threads());

    FILE *overhead_csv_file = fopen(overhead_csv_path, "a");
    if (!overhead_csv_file) {
        perror("fopen (overhead csv file)");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "matrix: %s\n", matrix_path);

    matrix_csr<val_t, rowptr_t, colidx_t> matrix;
    read_matrix(matrix, matrix_path);

#if CALCULATE_NNZ_PER_ROW_VARIANCE
    double variance = 0.0;

    double avg_nnzs_per_row = (double)matrix.nnz / matrix.nrow;
    for (uint r = 0; r < matrix.nrow; ++r) {
	int nnzs = matrix.row_ptr[r + 1] - matrix.row_ptr[r];
        variance += ((double)nnzs - avg_nnzs_per_row) * ((double)nnzs - avg_nnzs_per_row);
    }
    variance /= matrix.nrow;
    FILE *varfile = fopen("variance.csv", "a");
    fprintf(varfile, "%s,%f,%f\n", matrix_path, variance, variance / avg_nnzs_per_row);
    fclose(varfile);
    exit(0);
#endif

    constexpr int threads_per_shared_cache = 12;
    constexpr int num_shared_caches        = 4;

    assert(omp_get_max_threads() <= MAX_THREADS);
    assert((cline<int32_t, 256>(0u) == 0));
    assert((cline<int32_t, 256>(64u) == 1));
    assert((cline<int32_t, 256>(128u) == 2));
    assert((cline<double, 256>(0u) == 0));
    assert((cline<double, 256>(32u) == 1));
    assert((cline<double, 256>(64u) == 2));

    // set_buckets();
#if USE_SCALED_REUSE
    set_buckets_a64fx_scaled(matrix);
#else
    set_buckets_a64fx(matrix);
#endif /* USE_SCALED_REUSE */

    std::array<SharedCache, num_shared_caches> shared_caches{};

    fprintf(csv_file, "matrix,nnz,nrows,cache_id,shared,mindist,count\n");
    double time;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        PrivateCache pc{};
        SharedCache &sc = shared_caches[tid / threads_per_shared_cache];

        if (verbose) {
#pragma omp barrier
#pragma omp single
            time = omp_get_wtime();
        }

        for (int rep = 0; rep < 2; ++rep) {
#if USE_SCALED_REUSE
            reuse_scaled(tid, pc, sc, matrix);
#else
            reuse_calc(tid, pc, sc, matrix);
#endif /* USE_SCALED_REUSE */

            if (rep == 0) {
                pc.reset_buckets();
                sc.reset_buckets_shared(tid);
#pragma omp barrier
#pragma omp single
                time = omp_get_wtime();
            }
        } /* parallel */

        if (verbose) {
#pragma omp barrier
#pragma omp single
            {
                fprintf(stderr, "matrix: %s, time: %f sec\n", matrix_path, omp_get_wtime() - time);
                fprintf(overhead_csv_file, "%s, %f\n", matrix_path, omp_get_wtime() - time);
            }
        }

#pragma omp critical
        pc.print_csv(csv_file, matrix, tid);
    }

    size_t i = 0u;
    for (auto &sc : shared_caches) {
        sc.print_csv(csv_file, matrix, i);
        ++i;
    }
    fclose(csv_file);
    fclose(overhead_csv_file);
}
