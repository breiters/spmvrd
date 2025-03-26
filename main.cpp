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

#ifndef USE_ONLY_X_IN_TEMPORAL_SECTOR
#    define USE_ONLY_X_IN_TEMPORAL_SECTOR 1
#endif /* USE_ONLY_X_IN_TEMPORAL_SECTOR */

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

#define WANT_L1_CACHE_MISSES 1
#if WANT_L1_CACHE_MISSES
    for (int i = 0; i < L1ways; i++)
        Bucket::min_dists.push_back(L1d_capacity_per_way * (i + 1) / MEMBLOCKLEN);
#endif

    for (int i = 0; i < L2ways; i++)
        Bucket::min_dists.push_back(L2_capacity_per_way * (i + 1) / MEMBLOCKLEN);

    // bucket for cold misses (infinite reuse distance)
    Bucket::min_dists.push_back(Bucket::INF_DIST);

    // remove duplicated buckets (if any)
    auto &vec = Bucket::min_dists;
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

    // sort buckets in ascending order
    std::sort(Bucket::min_dists.begin(), Bucket::min_dists.end());
}

void reuse_sector0(int tid, PrivateCache &pc, SharedCache &sc, const auto &matrix)
{
#pragma omp for schedule(static)
    for (unsigned r = 0; r < matrix.nrow; ++r) {
        // fprintf(stderr, "row: %d\n", r);
        for (rowptr_t i = matrix.row_ptr[r]; i < matrix.row_ptr[r + 1]; ++i) {
            auto cl_x = cline<val_t, MEMBLOCKLEN>(matrix.col_idx[i]);
            // printf("row: %u, coldix: %u, cline: %lu val: %f i: %u\n", r,
            // matrix.col_idx[i], cl_x, matrix.val[i], i);
#if WANT_L1_CACHE_MISSES
            pc.handle_cline(cl_x);
#endif
            sc.handle_cline_shared(tid, cl_x);
        }
        /* TODO: this should be incremented in handle_cline */
#if WANT_L1_CACHE_MISSES
        pc.row_count_++;
#endif
        sc.row_count_++;
    }
}

void reuse_calc2(int tid, PrivateCache &pc, SharedCache &sc, const auto &matrix)
{
    //           x[0]...x[ncol] <padding> rowptr[0] ... rowptr[nrow] ...
    //
    // cacheline(x[0]) = 0 ... cacheline(ncol) = ncol * sizeof(val_t) / MEMBLOCKLEN ...
    auto cl_x_end     = cline<val_t, MEMBLOCKLEN>(matrix.ncol);
    auto cl_row_start = cl_x_end + 1;
    auto cl_y_start   = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(matrix.nrow + 1) + 1;
    auto cl_a_start   = cl_y_start + cline<val_t, MEMBLOCKLEN>(matrix.nrow) + 1;
    auto cl_col_start = cl_a_start + cline<val_t, MEMBLOCKLEN>(matrix.nnz) + 1;

#if !USE_ONLY_X_IN_TEMPORAL_SECTOR || USE_CALC_NOSC_REUSE
    // row ptr[0]
    auto cl_row = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(0);
    pc.handle_cline(cl_row);
    sc.handle_cline_shared(tid, cl_row);
#endif

#pragma omp for schedule(static)
    for (unsigned r = 0; r < matrix.nrow; ++r) {
        // fprintf(stderr, "row: %d\n", r);
#if !USE_ONLY_X_IN_TEMPORAL_SECTOR || USE_CALC_NOSC_REUSE
        // rowptr[r + 1]
        auto cl_row_plus1 = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(r + 1);
        pc.handle_cline(cl_row_plus1);

        // y[r]
        auto cl_y = cl_y_start + cline<val_t, MEMBLOCKLEN>(r);
        pc.handle_cline(cl_y);

        sc.handle_clines_shared(tid, cl_row_plus1, cl_y);
#endif

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

void reuse_sector1(int tid, PrivateCache &pc, SharedCache &sc, const auto &matrix)
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

#pragma omp for schedule(static)
    for (unsigned r = 0; r < matrix.nrow; ++r) {
        // rowptr[r + 1]
        auto cl_row_plus1 = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(r + 1);
        pc.handle_cline(cl_row_plus1);

        // y[r]
        auto cl_y = cl_y_start + cline<val_t, MEMBLOCKLEN>(r);
        pc.handle_cline(cl_y);
        sc.handle_clines_shared(tid, cl_row_plus1, cl_y);

        for (rowptr_t i = matrix.row_ptr[r]; i < matrix.row_ptr[r + 1]; ++i) {
            // a[i]
            auto cl_a = cl_a_start + cline<val_t, MEMBLOCKLEN>(i);
            pc.handle_cline(cl_a);
            // col_idx[i]
            auto cl_col = cl_col_start + cline<colidx_t, MEMBLOCKLEN>(i);
            pc.handle_cline(cl_col);

            sc.handle_clines_shared(tid, cl_a, cl_col);
        }
    }
}

enum Config : unsigned { SECTOR0 = 0, SECTOR1, CONFIG_END };

int main(int argc, char *argv[])
{
    char       *matrix_path = nullptr;
    FILE       *csv_file    = stdout;
    bool        verbose     = false;
    enum Config conf        = SECTOR0;

    int opt;
    while ((opt = getopt(argc, argv, "f:o:c:v")) != -1) {
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

        case 'c':
            conf = static_cast<Config>(strtoul(optarg, nullptr, 10));
            if (conf >= CONFIG_END) {
                goto usage;
            }
            break;

        case 'v':
            verbose = true;
            break;

        default: /* '?' */
usage:
            fprintf(stderr, "Usage: %s -f <matrix file>  [-o csv file] [-c config] [-v]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (!matrix_path)
        goto usage;

    char overhead_csv_path[1024];
    snprintf(overhead_csv_path, 1024, "overhead-%03dthreads.csv", omp_get_max_threads());

    FILE *overhead_csv_file = fopen(overhead_csv_path, "a");
    if (!overhead_csv_file) {
        perror("fopen (overhead csv file)");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "matrix: %s\n", matrix_path);

    matrix_csr<val_t, rowptr_t, colidx_t> matrix;
    read_matrix(matrix, matrix_path);

    // matrix values are not required ==> free to make space for reuse distance algorithm
    free(matrix.val);
    matrix.val = nullptr;

    char *needle = strrchr(matrix_path, '/');
    matrix.name  = needle ? needle + 1 : matrix_path;

    constexpr int threads_per_shared_cache = 12;
    constexpr int num_shared_caches        = 4;

    assert(omp_get_max_threads() <= MAX_THREADS);

    assert((cline<int32_t, 256>(0u) == 0));
    assert((cline<int32_t, 256>(64u) == 1));
    assert((cline<int32_t, 256>(128u) == 2));
    assert((cline<double, 256>(0u) == 0));
    assert((cline<double, 256>(32u) == 1));
    assert((cline<double, 256>(64u) == 2));

#if USE_SCALED_REUSE
    set_buckets_a64fx_scaled(matrix);
#else
    set_buckets_a64fx(matrix);
#endif /* USE_SCALED_REUSE */

    std::array<SharedCache, num_shared_caches> shared_caches{};

    // fprintf(csv_file, "matrix,nnz,nrows,cache_id,shared,mindist,count\n");
    fprintf(csv_file, Cache::csv_header_);
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
            switch (conf) {
            case SECTOR0:
                reuse_sector0(tid, pc, sc, matrix);
                break;
            case SECTOR1:
                reuse_sector1(tid, pc, sc, matrix);
                break;
            default:
                /* unreachable */
                break;
            }

            if (rep == 0) {
                pc.reset_buckets();
                sc.reset_buckets_shared(tid);
#pragma omp barrier
#pragma omp single
                time = omp_get_wtime();
            }
        }

        if (verbose) {
#pragma omp barrier
#pragma omp single
            {
                double time_diff = omp_get_wtime() - time;
                fprintf(stderr, "matrix: %s, time: %f sec\n", matrix_path, time_diff);
                fprintf(overhead_csv_file, "%s, %f\n", matrix_path, time_diff);
            }
        }

#pragma omp critical
        pc.print_csv(csv_file, matrix, tid);
    } /* parallel */

    size_t i = 0u;
    for (auto &sc : shared_caches) {
        sc.print_csv(csv_file, matrix, i);
        ++i;
    }
    fclose(csv_file);
    fclose(overhead_csv_file);
}
