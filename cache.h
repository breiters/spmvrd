#pragma once

#include "bucket.h"
#include "config.h"
#include "mcslock.h"

#include <algorithm>
#include <atomic>
#include <list>
#include <unordered_map>
#include <vector>

#include <cstdint>
#include <cstdio>

typedef unsigned Addr;

// make sure that memblocks are powers of two
static constexpr bool is_pow2(int a) { return !(a & (a - 1)); }
static_assert(is_pow2(MEMBLOCKLEN), "is_pow2(MEMBLOCKLEN)");
static_assert(is_pow2(CACHE_LINESIZE), "is_pow2(MEMBLOCKLEN)");

// find first bit set
static constexpr int ffs_constexpr(int x)
{
    int n = 0;
    while ((x & 1) == 0) {
        x >>= 1;
        n++;
    }
    return n;
}

static_assert(ffs_constexpr(256) == 8, "ffs_constexpr(256) == 8"); // sanity check

// gets the cache line number of an index
template <typename T, size_t CLSIZE>
unsigned cline(uint64_t idx)
{
    constexpr static auto first_bit_set = ffs_constexpr(CLSIZE / sizeof(T));
    return idx >> first_bit_set;
}

class Cache
{
public:
    void handle_cline(Addr addr)
    {
        if (addr == last_) {
            incr_access({0u, 0u, 0u});
            stack_.begin()->nnz_count = nnz_count_;
            nnz_count_++;
            return;
        }
        last_ = addr;

        auto map_it = refmap_.find(addr);

        if (map_it == refmap_.end()) {
            refmap_[addr] = on_block_new(MemoryBlock{0u, row_count_, nnz_count_});
            incr_access_inf();
        } else {
            incr_access(on_block_seen(map_it->second));
            refmap_[addr]->nnz_count = nnz_count_;
            refmap_[addr]->row_count = row_count_;
        }

        nnz_count_++;
    }

    StackIterator on_block_new(MemoryBlock &&mb);

    void incr_access(Bucket::Counts &&bucket)
    {
        auto [bx, bxy, bxya] = bucket;
        buckets_[bx].access_counts.count_x++;
        buckets_[bxy].access_counts.count_xy++;
        buckets_[bxya].access_counts.count_xya++;
    }

    void incr_access_inf()
    {
        auto bucket_inf = Bucket::min_dists.size() - 1;
        incr_access({bucket_inf, bucket_inf, bucket_inf});
    }

#if 0
    size_t bounded_distance(StackIterator first, StackIterator last, size_t bound)
    {
        size_t result{0u};
        while (result != bound && first != last) {
            ++first;
            ++result;
        }
        return result;
    }
#endif

    Bucket::Counts on_block_seen(StackIterator &it)
    {
        reuse_count_++;

        // if already on top of stack: do nothing (bucket is zero anyway)
        if (it == stack_.begin()) {
            return {0u, 0u, 0u};
        }

        // compute additional distance due to y and rowptr
        auto ncl_row_y = (row_count_ - it->row_count) * 16 / 256;

        // compute additional distance due to a and colidx
        auto ncl_col_a = (nnz_count_ - it->nnz_count) * 12 / 256;

        // get current bucket
        unsigned bucket = it->bucket;

        Bucket::Counts result = {bucket, bucket, bucket};

        size_t reuse_distance_xy  = ncl_row_y;
        size_t reuse_distance_xya = ncl_row_y + ncl_col_a;

        unsigned bucket_xy  = bucket;
        unsigned bucket_xya = bucket;

        size_t rd_min = Bucket::min_dists[bucket];
        size_t rd_max = Bucket::min_dists[bucket + 1] - 1;
        size_t rd_avg = (rd_min + rd_max) / 2;

        size_t bucket_before_inf = buckets_.size() - 2;
        if (bucket >= bucket_before_inf) {
            result = {bucket_before_inf, bucket_before_inf, bucket_before_inf};
        } else {
            reuse_distance_xy += rd_avg;
            reuse_distance_xya += rd_avg;

            // compute buckets of policy xy and policy xya
            for (size_t b = bucket + 1; b < buckets_.size(); ++b) {
                auto min_distance = Bucket::min_dists[b];
                if (min_distance >= reuse_distance_xy) {
                    ++bucket_xy;
                }
                if (min_distance >= reuse_distance_xya) {
                    ++bucket_xya;
                }
            }

            result = {bucket, bucket_xy, bucket_xya};
        }

#if 0
        auto next_marker = buckets_[bucket + 1].marker;

        // if already in last bucket before infinite reuse distance, return
        auto bucket_before_inf = buckets_.size() - 2;
        if (bucket >= bucket_before_inf) {
            result = {bucket_before_inf, bucket_before_inf, bucket_before_inf};
            goto out;
        } else if (Bucket::min_dists[bucket] + ncl_row_y >= Bucket::min_dists[bucket_before_inf]) {
            // if min reuse distance of xy is in last bucket before infinite reuse distance, return
            result = {bucket, bucket_before_inf, bucket_before_inf};
            goto out;
        }

        // else get exact reuse distance

        // first get distance to next marker in stack
        size_t distance;

        // exact reuse distance is either:
        // - next marker's min. reuse distance minus distance (forward)
        // - current marker's min. reuse distance plus distance (backward)
        size_t reuse_distance_x;

        // stack end marks that marker is not yet set
        if (next_marker != stack_.end()) {
            // forward distance
            distance = bounded_distance(it, next_marker, ncl_row_y + ncl_col_a);

            // fast path if next marker not reached
            if (distance < (ncl_row_y + ncl_col_a)) {
                result = {bucket, bucket, bucket};
                goto out;
            }

            reuse_distance_x = Bucket::min_dists[bucket + 1] - distance;
        } else {
            // backward distance
            // distance = bounded_distance(buckets_[bucket].marker, it, ncl_row_y + ncl_col_a);
            distance = std::distance(buckets_[bucket].marker, it);

            // fast path if next marker not reached
            if (distance >= (ncl_row_y + ncl_col_a)) {
                result = {bucket, bucket, bucket};
                goto out;
            }

            reuse_distance_x = Bucket::min_dists[bucket] + distance;
        }

        reuse_distance_xy += reuse_distance_x;
        reuse_distance_xya += reuse_distance_x;

        // compute buckets of policy xy and policy xya
        for (auto b = bucket + 1; b < buckets_.size(); ++b) {
            auto min_distance = Bucket::min_dists[b];
            if (min_distance <= reuse_distance_xy) {
                ++bucket_xy;
            }
            if (min_distance <= reuse_distance_xya) {
                ++bucket_xya;
            }
        }

        // printf("reuse distances: %zu, %zu, %zu\n", reuse_distance_x, reuse_distance_xy, reuse_distance_xya);
        // printf("buckets: %zu, %zu, %zu\n", bucket_x, bucket_xy, bucket_xya);

out:
#endif

        // then move all markers below current memory block's bucket
        move_markers(bucket);

        // put current memory block on top of stack
        stack_.splice(stack_.begin(), stack_, it);

        // bucket of blockIt is zero now because it is on top of stack
        it->bucket = 0u;

#if RD_DEBUG > 1
        check_consistency();
#endif /* RD_DEBUG */

        return result;
    }

#if 0
    void incr_access(unsigned bucket) { buckets_[bucket].access_count++; }
    // count access with infinite reuse distance
    void incr_access_inf() { incr_access(buckets_.size() - 1); }
#endif

    // void increment_row_count() { row_count_++; }

    static constexpr const char *csv_header_ =
        "matrix,nnz,nrows,cache_id,shared,time,working_set_size,reuses,mindist,count_x,count_xy,count_xya\n";
    void print_csv(FILE *file, const auto &matrix, int id, double time) const
    {
        size_t working_set_size = refmap_.size();
        for (size_t i = 0u; i != Bucket::min_dists.size(); ++i) {
            // matrix name, nnz, nrow, cache id, shared, min bucket, count
            fprintf(file,
                    "%s,%zu,%zu,%d,%d,%f,%zu,%zu,%lu,%lu,%lu,%lu\n",
                    matrix.name,
                    matrix.nnz,
                    matrix.nrow,
                    id,
                    shared_,
                    time,
                    working_set_size,
                    reuse_count_,
                    Bucket::min_dists[i],
                    buckets_[i].access_counts.count_x,
                    buckets_[i].access_counts.count_xy,
                    buckets_[i].access_counts.count_xya);
        }
    }

    void reset_buckets()
    {
        for (auto &b : buckets_) {
            b.access_counts = {0u, 0u, 0u};
        }
    }

    std::atomic<uint32_t> row_count_{0u};

private:
    void move_markers(unsigned);
    void on_next_bucket_gets_active();
    void check_consistency();

    std::list<MemoryBlock>                  stack_{};
    std::unordered_map<Addr, StackIterator> refmap_{};

    uint64_t nnz_count_{0u};
    uint64_t reuse_count_{0u};

    Addr     last_{(Addr)-1};
    unsigned next_bucket_{1u};

    std::vector<Bucket> buckets_{
        std::vector<Bucket>{Bucket::min_dists.size(), stack_.end()}
    };

protected:
    bool shared_{false};
};

class SharedCache : public Cache
{
public:
    SharedCache() { shared_ = true; }

    void handle_cline_shared(int tid, Addr a)
    {
        mcslock_.lock(tid);
        handle_cline(a);
        mcslock_.unlock(tid);
    }

    void handle_clines_shared(int tid, Addr a0, Addr a1)
    {
        mcslock_.lock(tid);
        handle_cline(a0);
        handle_cline(a1);
        mcslock_.unlock(tid);
    }

    void handle_clines_shared(int tid, Addr a0, Addr a1, Addr a2)
    {
        mcslock_.lock(tid);
        handle_cline(a0);
        handle_cline(a1);
        handle_cline(a2);
        mcslock_.unlock(tid);
    }

#if 0
    template<typename... As>
    void handle_clines_shared(int tid, As... addrs)
    {
        mcslock_.lock(tid);
        handle_clines_shared(addrs...);
        mcslock_.unlock(tid);
    }

    template<typename... As>
    void handle_clines_shared(Addr a, As... addrs)
    {
        handle_cline(a);
        handle_clines_shared(addrs...);
    }

    void handle_clines_shared(Addr a)
    {
        handle_cline(a);
    }
#endif

    void reset_buckets_shared(int tid)
    {
        mcslock_.lock(tid);
        reset_buckets();
        mcslock_.unlock(tid);
    }

private:
    MCSLock mcslock_{};
};

class PrivateCache : public Cache
{
};
