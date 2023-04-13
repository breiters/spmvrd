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

typedef uint Addr;

// make sure that memblocks are powers of two
static constexpr bool is_pow2(int a) { return !(a & (a - 1)); }
static_assert(is_pow2(MEMBLOCKLEN), "is_pow2(MEMBLOCKLEN)");

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

// gets the cache line number of a column index
template <typename T, size_t CLSIZE>
uint cline(T colidx)
{
    constexpr static auto fs = ffs_constexpr(CLSIZE / sizeof(T));
    return colidx >> fs;
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
        static Addr last{(Addr)-1};
        if(addr == last) {
            incr_access(0);
            return;
        }
        last = addr;

        auto map_it = refmap.find(addr);

        if (map_it == refmap.end()) {
            refmap[addr] = on_block_new(MemoryBlock{});
            incr_access_inf();
        } else {
            int bucket = on_block_seen(map_it->second);
            incr_access(bucket);
        }
    }

    StackIterator on_block_new(MemoryBlock &&mb);

    int on_block_seen(StackIterator &it);

    void incr_access(uint bucket) { buckets_[bucket].access_count++; }
    // count access with infinite reuse distance
    void incr_access_inf() { incr_access(buckets_.size() - 1); }

    StackIterator stack_begin() { return stack_.begin(); }

    const std::vector<Bucket> &buckets() const { return buckets_; }

    void print_csv(FILE *file, const char *matrix_name, int id) const
    {
        size_t i = 0u;
        for (auto &b : buckets_) {
            // matrix name, cache id, shared, min bucket, count
            fprintf(file, "%s,%d,%d,%lu,%lu\n", matrix_name, id, shared_, Bucket::min_dists[i], b.access_count);
            ++i;
        }
    }

private:
    void move_markers(uint);
    void on_next_bucket_gets_active();
    void check_consistency();

    std::list<MemoryBlock>                  stack_{}; // Stack structure with MemoryBlock as element
    std::unordered_map<Addr, StackIterator> refmap{};

    uint                next_bucket_{1u};                                        //
    std::vector<Bucket> buckets_{std::vector<Bucket>{Bucket::min_dists.size()}}; //

protected:
    bool shared_{false};
};

class SharedCache : public Cache
{
public:
    SharedCache() { shared_ = true; }

    void handle_cline_shared(Addr addr, int tid)
    {
        mcslock.lock(tid);
        handle_cline(addr);
        mcslock.unlock(tid);
    }

private:
    MCSLock mcslock{};
};

class PrivateCache : public Cache
{
};
