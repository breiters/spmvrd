#pragma once

#include "bucket.h"
#include "config.h"
#include "mcslock.h"

#include <algorithm>
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
static constexpr int first_set = ffs_constexpr(MEMBLOCKLEN);

// gets the cache line number of an index
template <typename T, size_t CLSIZE>
unsigned cline(uint64_t idx)
{
    constexpr static auto fs = ffs_constexpr(CLSIZE / sizeof(T));
    return idx >> fs;
}

struct alignas(CACHE_LINESIZE) AddrWrapper {
    Addr addr = 0;
};
static_assert(sizeof(AddrWrapper) == CACHE_LINESIZE, "sizeof(AddrWrapper) != CACHE_LINESIZE");

class Cache
{
public:
    void handle_cline(Addr addr)
    {
        if (addr == last_) {
            incr_access(0);
            return;
        }
        last_ = addr;

        auto map_it = refmap_.find(addr);

        if (map_it == refmap_.end()) {
            refmap_[addr] = on_block_new(MemoryBlock{});
            incr_access_inf();
        } else {
            int bucket = on_block_seen(map_it->second);
            incr_access(bucket);
        }
    }

    StackIterator on_block_new(MemoryBlock &&mb);

    int on_block_seen(StackIterator &it);

    void incr_access(unsigned bucket) { buckets_[bucket].access_count++; }
    // count access with infinite reuse distance
    void incr_access_inf() { incr_access(buckets_.size() - 1); }

    const std::vector<Bucket> &buckets() const { return buckets_; }

    void print_csv(FILE *file, const auto &matrix, int id) const
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

    void reset_buckets()
    {
        for (auto &b : buckets_)
            b.access_count = 0;
    }

private:
    void move_markers(unsigned);
    void on_next_bucket_gets_active();
    void check_consistency();

    std::list<MemoryBlock>                  stack_{};
    std::unordered_map<Addr, StackIterator> refmap_{};

    Addr                last_{(Addr)-1};
    unsigned                next_bucket_{1u};                                        //
    std::vector<Bucket> buckets_{std::vector<Bucket>{Bucket::min_dists.size()}}; //

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
