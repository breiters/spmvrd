#pragma once

#include "bucket.h"
#include "config.h"
#include "mcslock.h"

#include <algorithm>
#include <list>
// #include <unordered_map>
#include <vector>

#include <cstdint>
#include <cstdio>

typedef unsigned Addr;

// make sure that memblocks are powers of two
static constexpr bool is_pow2(int a) { return !(a & (a - 1)); }
static_assert(is_pow2(MEMBLOCKLEN), "is_pow2(MEMBLOCKLEN)");
static_assert(is_pow2(CACHE_LINESIZE), "is_pow2(CACHE_LINESIZE)");

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
Addr cline(uint64_t idx)
{
    constexpr static auto first_bit_set = ffs_constexpr(CLSIZE / sizeof(T));
    return static_cast<Addr>(idx >> first_bit_set);
}

class Cache
{
public:
    void handle_cline(Addr addr)
    {
        if (addr == last_) {
            incr_access({0u, 0u, 0u});
            stack_.begin()->nnz_count = nnz_count_;
            stack_.begin()->row_count = row_count_;
            nnz_count_++;
            reuse_count_++;
            return;
        }
        last_ = addr;

        // auto map_it = refmap_.find(addr);
        StackIterator &it = refmap_[addr];

        // if (map_it == refmap_.end()) {
        if (it == stack_.end()) {
            incr_access_inf();
            refmap_[addr] = on_block_new(MemoryBlock{0u, row_count_, nnz_count_});
        } else {
            // incr_access(on_block_seen(map_it->second));
            incr_access(on_block_seen(it));
            refmap_[addr]->nnz_count = nnz_count_;
            refmap_[addr]->row_count = row_count_;
            reuse_count_++;
        }

        nnz_count_++;
    }

    void handle_cline(Addr addr, bool increment_row)
    {
        row_count_++;
        handle_cline(addr);
    }

    Bucket::Counts on_block_seen(StackIterator &it);
    StackIterator  on_block_new(MemoryBlock &&mb);

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

    static constexpr const char *csv_header_ =
        "matrix,nnz,nrows,cache_id,shared,time,nnz_count,working_set_size,reuses,mindist,count_x,count_xy,count_xya\n";
    void print_csv(FILE *file, const auto &matrix, int id, double time) const
    {
        size_t working_set_size = stack_.size();

        for (size_t i{0u}; i != Bucket::min_dists.size(); ++i) {
            // matrix name, nnz, nrow, cache id, shared, min bucket, count
            fprintf(file,
                    "%s,%zu,%zu,%d,%d,%f,%zu,%zu,%zu,%lu,%lu,%lu,%lu\n",
                    matrix.name,
                    matrix.nnz,
                    matrix.nrow,
                    id,
                    shared_,
                    time,
                    nnz_count_ / 2, // TODO: hotfix since we are doing two rounds, thus the nnz count is doubled
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
        reuse_count_ = 0u;
    }

    void set_refmap_size(size_t nelem)
    {
        size_t nlines = nelem * sizeof(double) / MEMBLOCKLEN + 1;
#pragma omp critical
        {
            if (refmap_.empty()) {
                refmap_.reserve(nlines);
                for (size_t i{0u}; i != nlines; ++i) {
                    refmap_[i] = stack_.end();
                }
            }
        }
    }

private:
    void move_markers(unsigned bucket_max);
    void on_next_bucket_gets_active();
    void check_consistency(bool force);

    std::list<MemoryBlock> stack_{};
    // std::unordered_map<Addr, StackIterator> refmap_{};

    /* we use a vector as reference map because it's faster */
    /* this is possible due to the virtual manually-assigned cache line numbers */
    std::vector<StackIterator> refmap_{};

    uint64_t nnz_count_{0u};
    uint64_t reuse_count_{0u};
    uint32_t row_count_{0u};

    unsigned next_bucket_{1u};
    Addr     last_{(Addr)-1};

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

    void handle_cline_shared(int tid, Addr a, bool increment_row)
    {
        mcslock_.lock(tid);
        handle_cline(a, increment_row);
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
