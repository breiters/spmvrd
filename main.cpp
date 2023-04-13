#include "cache.h"
#include "read_matrix.h"

#include <array>
#include <omp.h>
#include <string>

void set_buckets_a64fx()
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

    for (int i = 0; i < L1ways; i++)
        Bucket::min_dists.push_back(L1d_capacity_per_way * (i + 1) / MEMBLOCKLEN);

    for (int i = 0; i < L2ways; i++)
        Bucket::min_dists.push_back(L2_capacity_per_way * (i + 1) / MEMBLOCKLEN);

    // TODO: don't use unnecessary buckets for shared / private caches 
    // TODO: think about the scaling factor and the required buckets again
    int scale = (sizeof(double) + sizeof(uint32_t)) / sizeof(uint32_t);
    Bucket::min_dists.push_back(64 * KiB / scale);
    Bucket::min_dists.push_back(8 * MiB / scale);

    // bucket for cold misses (infinite reuse distance)
    Bucket::min_dists.push_back(Bucket::inf_dist);

    // sort buckets in ascending order
    std::sort(Bucket::min_dists.begin(), Bucket::min_dists.end());
}

void set_buckets()
{
    unsigned long KiB = 1024;

    Bucket::min_dists.push_back(0);

    for (int i = 0; i < 20; i++) {
        Bucket::min_dists.push_back(KiB / MEMBLOCKLEN);
        KiB *= 1.5;
    }

    // bucket for cold misses (infinite reuse distance)
    Bucket::min_dists.push_back(Bucket::inf_dist);

    // sort buckets in ascending order
    std::sort(Bucket::min_dists.begin(), Bucket::min_dists.end());
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
            fprintf(stderr, "Usage: %s -f <matrix file>  [-o csv file] [-v]", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (!matrix_path)
        goto usage;

    char *matrix_name = matrix_path; // may want to strip first

    auto matrix = matrix_csr<double, uint32_t, uint32_t>::read_matrix(matrix_path);

    constexpr int threads_per_shared_cache = 12;
    constexpr int num_shared_caches        = 4;

    assert(omp_get_max_threads() <= MAX_THREADS);

    set_buckets();
    // set_buckets_a64fx();
    std::array<SharedCache, num_shared_caches> shared_caches{};

    fprintf(csv_file, "matrix,cache_id,shared,mindist,count\n");
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

#pragma omp for
        for (uint r = 0; r < matrix.nrow; ++r) {
            for (uint i = matrix.row_ptr[r]; i < matrix.row_ptr[r + 1]; ++i) {
                // printf("row: %u, coldix: %u, val: %lf\n", r, matrix.col_idx[i], matrix.val[i]);
                auto cl = cline<uint32_t, MEMBLOCKLEN>(matrix.col_idx[i]);
                pc.handle_cline(cl);
                sc.handle_cline_shared(cl, tid);
            }
        }

        if (verbose) {
#pragma omp barrier
#pragma omp single
            fprintf(stderr, "time: %lf sec\n", omp_get_wtime() - time);
        }

#pragma omp critical
        pc.print_csv(csv_file, matrix_name, tid);
    }

    size_t i = 0u;
    for (auto &c : shared_caches) {
        c.print_csv(csv_file, matrix_name, i);
        ++i;
    }
    fclose(csv_file);
}