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

#ifndef RD_POLICY_X
#    define RD_POLICY_X 1
#endif /* RD_POLICY_X */

#ifndef RD_WANT_L1_CACHE_MISSES
#    define RD_WANT_L1_CACHE_MISSES 0
#endif /* RD_WANT_L1_CACHE_MISSES */

#ifndef RD_WANT_DOUBLE_RESOLUTION
#    define RD_WANT_DOUBLE_RESOLUTION 0
#endif /* RD_WANT_DOUBLE_RESOLUTION */

using rowptr_t = int64_t;
using colidx_t = int;
using val_t    = double;

enum CachePolicy { POLICY_X = 0, POLICY_XY, POLICY_XYA, POLICY_MAX };

void set_buckets_a64fx(void)
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

#if RD_WANT_L1_CACHE_MISSES
    for (int i = 0; i < L1ways; i++)
        Bucket::min_dists.push_back(L1d_capacity_per_way * (i + 1) / MEMBLOCKLEN);
#endif

    for (int i = 0; i < L2ways; i++)
        Bucket::min_dists.push_back(L2_capacity_per_way * (i + 1) / MEMBLOCKLEN);

#if RD_WANT_DOUBLE_RESOLUTION
    // double resolution in most relevant region
    for (int i = 0; i < L2ways; i++)
        Bucket::min_dists.push_back((L2_capacity_per_way * (i + 1) - L2_capacity_per_way / 2) / MEMBLOCKLEN);
#endif /* RD_WANT_DOUBLE_RESOLUTION */

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
        for (rowptr_t i = matrix.row_ptr[r]; i < matrix.row_ptr[r + 1]; ++i) {
            auto cl_x             = cline<val_t, MEMBLOCKLEN>(matrix.col_idx[i]);
            bool first_nnz_in_row = (i == matrix.row_ptr[r]);
#if RD_WANT_L1_CACHE_MISSES
            pc.handle_cline(cl_x, first_nnz_in_row);
#endif
            sc.handle_cline_shared(tid, cl_x, first_nnz_in_row);
        }
    }
}

// sector0 : x / xy / xya
// sector1 : ya / a / -

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

#if !RD_POLICY_X || RD_POLICY_XYA
    // row_ptr[r]
    unsigned first_row = (matrix.nrow / omp_get_num_threads()) * omp_get_thread_num();

    auto cl_row = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(first_row);
    pc.handle_cline(cl_row);
    sc.handle_cline_shared(tid, cl_row);
#endif

#pragma omp for schedule(static)
    for (unsigned r = 0; r < matrix.nrow; ++r) {
        // fprintf(stderr, "row: %d\n", r);
#if !RD_POLICY_X || RD_POLICY_XYA
        // rowptr[r + 1]
        auto cl_row_plus1 = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(r + 1);
        pc.handle_cline(cl_row_plus1);

        // y[r]
        auto cl_y = cl_y_start + cline<val_t, MEMBLOCKLEN>(r);
        pc.handle_cline(cl_y);

        sc.handle_clines_shared(tid, cl_row_plus1, cl_y);
#endif

        for (rowptr_t i = matrix.row_ptr[r]; i < matrix.row_ptr[r + 1]; ++i) {
#if RD_POLICY_XYA
            // a[i]
            auto cl_a = cl_a_start + cline<val_t, MEMBLOCKLEN>(i);
            pc.handle_cline(cl_a);
            // col_idx[i]
            auto cl_col = cl_col_start + cline<colidx_t, MEMBLOCKLEN>(i);
            pc.handle_cline(cl_col);
#endif /* RD_POLICY_XYA */
            // x[col_idx[i]]
            auto cl_x = cline<val_t, MEMBLOCKLEN>(matrix.col_idx[i]);
            pc.handle_cline(cl_x);
#if RD_POLICY_XYA
            sc.handle_clines_shared(tid, cl_a, cl_col, cl_x);
#else
            sc.handle_cline_shared(tid, cl_x);
#endif /* RD_POLICY_XYA */
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

    // row_ptr[r]
    unsigned first_row = (matrix.nrow / omp_get_num_threads()) * omp_get_thread_num();

    auto cl_row = cl_row_start + cline<rowptr_t, MEMBLOCKLEN>(first_row);
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

#if !NDEBUG
    fprintf(stderr, "[!!] running %s in debug mode [!!]\n", argv[0]);
#endif

    if (omp_get_max_threads() > MAX_THREADS) {
        fprintf(stderr, "Error: %s configured for max. %d threads\n", argv[0], MAX_THREADS);
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "reading matrix: %s ...", matrix_path);
    matrix_csr<val_t, rowptr_t, colidx_t> matrix;
    read_matrix(matrix, matrix_path);
    fprintf(stderr, " done!\n");

    // matrix values are not required ==> free to make space for reuse distance algorithm
    free(matrix.val);
    matrix.val = nullptr;

    char *needle = strrchr(matrix_path, '/');
    matrix.name  = needle ? needle + 1 : matrix_path;

    constexpr int threads_per_shared_cache = THREADS_PER_SHARED_CACHE;
    constexpr int num_shared_caches        = NUM_SHARED_CACHES;

    set_buckets_a64fx();

    std::array<SharedCache, num_shared_caches> shared_caches{};

    fprintf(csv_file, Cache::csv_header_);
    double time;
    double time_diff;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        PrivateCache pc{};
        SharedCache &sc = shared_caches[tid / threads_per_shared_cache];

        // TODO: set refmap size depending on policy
        pc.set_refmap_size(matrix.ncol);
        sc.set_refmap_size(matrix.ncol);

#pragma omp barrier
#pragma omp single
        time = omp_get_wtime();

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

#pragma omp barrier
#pragma omp single
        {
            time_diff = omp_get_wtime() - time;
            if (verbose) {
                fprintf(stderr, "matrix: %s, time: %f sec\n", matrix_path, time_diff);
            }
        }

#pragma omp critical
        pc.print_csv(csv_file, matrix, tid, time_diff);
    } /* parallel */

    size_t i = 0u;
    for (auto &sc : shared_caches) {
        sc.print_csv(csv_file, matrix, i, time_diff);
        ++i;
    }
    fclose(csv_file);
}
