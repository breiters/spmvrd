#include "cache.h"
#include "bucket.h"
#include "config.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <utility>

void Cache::on_next_bucket_gets_active()
{
    // set new buckets marker to end of stack first then set marker to last stack element
    buckets_[next_bucket_].marker = stack_.end();

    --(buckets_[next_bucket_].marker);

    // fprintf(stderr, "stack size: %d, next bucket: %d\n", stack_.size(), next_bucket_);

#if RD_DEBUG
    StackIterator it = stack_.begin();
    for (Bucket::min_type i = 0; i < Bucket::min_dists[next_bucket_]; i++)
        it++;

    assert(it == buckets_[next_bucket_].marker);
#endif /* RD_DEBUG */

    // last stack element is now in the next higher bucket
    (buckets_[next_bucket_].marker)->bucket++;

    assert((buckets_[next_bucket_].marker)->bucket == next_bucket_);

    next_bucket_++;

#if RD_DEBUG
    check_consistency();
#endif /* RD_DEBUG */
}

/**
 * @brief Adds new memory block to top of stack. Moves active bucket markers.
 * Adds next bucket if necessary.
 *
 * @param mb The to memory block.
 * @return The stack begin iterator
 */
StackIterator Cache::on_block_new(MemoryBlock &&mb)
{
    stack_.push_front(std::move(mb));

    // move markers upwards after inserting new block on stack
    move_markers(next_bucket_ - 1);

    // does another bucket get active?
    if (Bucket::min_dists[next_bucket_] != Bucket::INF_DIST && (stack_.size() > Bucket::min_dists[next_bucket_])) {
        on_next_bucket_gets_active();
    }

#if RD_DEBUG > 1
    check_consistency();
#endif /* RD_DEBUG */

    return stack_.begin();
}

/**
 * Sanity check:
 * - every active marker must be found in stack in the right order
 * - distance of bucket marker to stack begin must be equal to the min distance for the
 * bucket
 */
void Cache::check_consistency()
{
#if RD_DEBUG
    const size_t  DO_CHECK = 10;
    static size_t iter     = 0;
    iter++;
    if (iter < DO_CHECK) {
        return;
    }
    iter = 0;

    auto     it       = stack_.begin();
    unsigned distance = 0;
    for (unsigned b = 1; b < next_bucket_; b++) {
        for (; it != buckets_[b].marker; it++) {
            assert(it != stack_.end());
            distance++;
        }
        assert(distance == Bucket::min_dists[b]);
    }
#endif /* RD_DEBUG */
}

void Cache::move_markers(unsigned topBucket)
{
    assert(topBucket < next_bucket_);
    for (unsigned b = 1; b <= topBucket; b++) {
        assert(buckets_[next_bucket_].marker != stack_.begin());

        // decrement marker so it stays always on same distance to stack begin
        --(buckets_[b].marker);

        // increment bucket of memory block where current marker points to
        (buckets_[b].marker)->bucket++;
    }
}

Bucket::Counts Cache::on_block_seen(StackIterator &it)
{
    reuse_count_++;

    // if already on top of stack: do nothing (bucket is zero anyway)
    if (it == stack_.begin()) {
        return {0u, 0u, 0u};
    }

    // compute additional distance due to y and rowptr
    auto ncl_row_y = ((row_count_ - it->row_count) * 16) / 256;

    // compute additional distance due to a and colidx
    auto ncl_col_a = ((nnz_count_ - it->nnz_count) * 12) / 256;

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
            bool done         = true;
            auto min_distance = Bucket::min_dists[b];
            if (reuse_distance_xy >= min_distance) {
                ++bucket_xy;
                done = false;
            }
            if (reuse_distance_xya >= min_distance) {
                ++bucket_xya;
                done = false;
            }
            if (done) {
                break;
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
