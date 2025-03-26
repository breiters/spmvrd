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

// TODO: make bucket list a template parameter

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

    Bucket::Counts on_block_seen(StackIterator &it)
    {
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

        // if already in last bucket before infinite reuse distance, return
        auto bucket_before_inf = buckets_.size() - 2;
        if (bucket >= bucket_before_inf) {
            return {bucket_before_inf, bucket_before_inf, bucket_before_inf};
        }

        // else get exact reuse distance
        auto next_marker = buckets_[bucket + 1].marker;

        // first get distance to next marker in stack
        size_t distance;
        
        // exact reuse distance is either:
        // - next marker's min. reuse distance minus distance (forward)
        // - current marker's min. reuse distance plus distance (backward)
        size_t reuse_distance_x;

        if (next_marker != stack_.end()) {
            // forward distance
            distance = std::distance(it, next_marker);
            reuse_distance_x   = Bucket::min_dists[bucket + 1] - distance;
        } else {
            // backward distance
            distance = std::distance(buckets_[bucket].marker, it);
            reuse_distance_x   = Bucket::min_dists[bucket] + distance;
        }

        size_t reuse_distance_xy  = reuse_distance_x + ncl_row_y;
        size_t reuse_distance_xya = reuse_distance_x + ncl_row_y + ncl_col_a;

        // compute buckets of policy xy and policy xya
        unsigned bucket_xy  = bucket;
        unsigned bucket_xya = bucket;

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

        // then move all markers below current memory block's bucket
        move_markers(bucket);

        // put current memory block on top of stack
        stack_.splice(stack_.begin(), stack_, it);

        // bucket of blockIt is zero now because it is on top of stack
        it->bucket = 0u;

#if RD_DEBUG > 1
        check_consistency();
#endif /* RD_DEBUG */

        return {bucket, bucket_xy, bucket_xya};
    }

#if 0
    void incr_access(unsigned bucket) { buckets_[bucket].access_count++; }
    // count access with infinite reuse distance
    void incr_access_inf() { incr_access(buckets_.size() - 1); }
#endif

    // void increment_row_count() { row_count_++; }

    static constexpr const char *csv_header_ =
        "matrix,nnz,nrows,cache_id,shared,working_set,mindist,count_x,count_xy,count_xya\n";
    void print_csv(FILE *file, const auto &matrix, int id) const
    {
        size_t working_set_size = refmap_.size();
        for (size_t i = 0u; i != Bucket::min_dists.size(); ++i) {
            // matrix name, nnz, nrow, cache id, shared, min bucket, count
            fprintf(file,
                    "%s,%zu,%zu,%d,%d,%zu,%lu,%lu,%lu,%lu\n",
                    matrix.name,
                    matrix.nnz,
                    matrix.nrow,
                    id,
                    shared_,
                    working_set_size,
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

    Addr     last_{(Addr)-1};
    unsigned next_bucket_{1u};

    std::vector<Bucket> buckets_{std::vector<Bucket>{Bucket::min_dists.size(), stack_.end()}};

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
