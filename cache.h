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
            incr_access(std::make_tuple(0, 0, 0));
            return;
        }
        last_ = addr;

        auto map_it = refmap_.find(addr);

        if (map_it == refmap_.end()) {
            refmap_[addr] = on_block_new(MemoryBlock{});
            incr_access_inf();
        } else {
            auto bucket = on_block_seen(map_it->second);
            incr_access(bucket);
        }

        nnz_count_++;
    }

    StackIterator on_block_new(MemoryBlock &&mb);

    template <typename Buckets>
    void incr_access(Buckets &&bucket)
    {
        auto [bx, bxy, bxya] = bucket;
        buckets_x_[bx].access_count++;
        buckets_xy_[bxy].access_count++;
        buckets_xya_[bxya].access_count++;
    }

    void incr_access_inf()
    {
        auto bucket_inf = Bucket::min_dists.size() - 1;
        incr_access(std::make_tuple(bucket_inf, bucket_inf, bucket_inf));
    }

    std::tuple<int, int, int> on_block_seen(StackIterator &it)
    {
        // if already on top of stack: do nothing (bucket is zero anyway)
        if (it == stack_.begin()) {
            return {};
        }

        // compute additional distance due to y and rowptr
        auto ncl_row_y = (row_count_ - it->row_count) * 16 / 256;

        // compute additional distance due to a and colidx
        auto ncl_col_a = (nnz_count_ - it->nnz_count) * 12 / 256;

        // get depth in stack
        auto reuse_distance_x   = std::distance(stack_.begin(), it);
        auto reuse_distance_xy  = reuse_distance_x + ncl_row_y;
        auto reuse_distance_xya = reuse_distance_x + ncl_row_y + ncl_col_a;

        // compute buckets
        int bucket_x   = 0;
        int bucket_xy  = 0;
        int bucket_xya = 0;
        for (auto min_distance : Bucket::min_dists) {
            // skip first bucket
            if(min_distance == 0) continue;

            if (min_distance <= reuse_distance_x) {
                ++bucket_x;
            }
            if (min_distance <= reuse_distance_xy) {
                ++bucket_xy;
            }
            if (min_distance <= reuse_distance_xya) {
                ++bucket_xya;
            }
        }

        printf("reuse distances: %zu, %zu, %zu\n", reuse_distance_x, reuse_distance_xy, reuse_distance_xya);
        printf("buckets: %zu, %zu, %zu\n", bucket_x, bucket_xy, bucket_xya);

        // put current memory block on top of stack
        stack_.splice(stack_.begin(), stack_, it);

        // bucket of blockIt is zero now because it is on top of stack
        it->bucket = 0;

        return {bucket_x, bucket_xy, bucket_xya};
    }

#if 0
    void incr_access(unsigned bucket) { buckets_[bucket].access_count++; }
    // count access with infinite reuse distance
    void incr_access_inf() { incr_access(buckets_.size() - 1); }
#endif

    // void increment_row_count() { row_count_++; }

    void print_csv(FILE *file, const auto &matrix, int id) const
    {
        for (size_t i = 0u; i != Bucket::min_dists.size(); ++i) {
            // matrix name, nnz, nrow, cache id, shared, min bucket, count
            fprintf(file,
                    "%s,%zu,%zu,%d,%d,%lu,%lu,%lu,%lu\n",
                    matrix.name,
                    matrix.nnz,
                    matrix.nrow,
                    id,
                    shared_,
                    Bucket::min_dists[i],
                    buckets_x_[i].access_count,
                    buckets_xy_[i].access_count,
                    buckets_xya_[i].access_count
                );
        }
    }
#if 0
    void print_csv2(FILE *file, const auto &matrix, int id) const
    {
        size_t i = 0u;
        for (auto &b : buckets_) {
            // matrix name, nnz, nrow, cache id, shared, min bucket, count
            fprintf(file,
                    "%s,%zu,%zu,%d,%d,%lu,%lu\n",
                    matrix.name,
                    matrix.nnz,
                    matrix.nrow,
                    id,
                    shared_,
                    Bucket::min_dists[i],
                    b.access_count);
            ++i;
        }
    }
#endif
    void reset_buckets()
    {
        // for (auto &b : buckets_)
        // b.access_count = 0;

        for (auto &b : buckets_x_)
            b.access_count = 0;

        for (auto &b : buckets_xy_)
            b.access_count = 0;

        for (auto &b : buckets_xya_)
            b.access_count = 0;
    }

    std::atomic<uint32_t> row_count_{0u};
    std::atomic<uint64_t> nnz_count_{0u};

private:
    void move_markers(unsigned);
    void on_next_bucket_gets_active();
    void check_consistency();

    std::list<MemoryBlock>                  stack_{};
    std::unordered_map<Addr, StackIterator> refmap_{};

    Addr     last_{(Addr)-1};
    unsigned next_bucket_{1u};

    // TODO: tuple of vectors, or vector of tuples?
    using Buckets = std::tuple<std::vector<Bucket>, std::vector<Bucket>, std::vector<Bucket>>;
    std::vector<Bucket> buckets_x_{std::vector<Bucket>{Bucket::min_dists.size()}};
    std::vector<Bucket> buckets_xy_{std::vector<Bucket>{Bucket::min_dists.size()}};
    std::vector<Bucket> buckets_xya_{std::vector<Bucket>{Bucket::min_dists.size()}};

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
