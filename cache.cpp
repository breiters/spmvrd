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
 * @brief
 *
 * @param blockIt
 */
int Cache::on_block_seen(StackIterator &blockIt)
{
    // if already on top of stack: do nothing (bucket is zero anyway)
    if (blockIt == stack_.begin()) {
        return 0;
    }

    // move all markers below current memory blocks bucket
    int bucket = blockIt->bucket;
    move_markers(bucket);

    // put current memory block on top of stack
    stack_.splice(stack_.begin(), stack_, blockIt);

    // bucket of blockIt is zero now because it is on top of stack
    blockIt->bucket = 0;

#if RD_DEBUG > 1
    check_consistency();
#endif /* RD_DEBUG */

    return bucket;
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
    iter          = 0;
    auto it_start = stack_.begin();
    for (uint b = 1; b < next_bucket_; b++) {
        unsigned distance = 0;
        for (auto it = it_start; it != buckets_[b].marker; it++) {
            assert(it != stack_.end());
            distance++;
        }
        assert(distance == Bucket::min_dists[b]);
    }
#endif /* RD_DEBUG */
}

void Cache::move_markers(uint topBucket)
{
    assert(topBucket < next_bucket_);

    for (uint b = 1; b <= topBucket; b++) {
        assert(buckets_[next_bucket_].marker != stack_.begin());

        // decrement marker so it stays always on same distance to stack begin
        --(buckets_[b].marker);

        // increment bucket of memory block where current marker points to
        (buckets_[b].marker)->bucket++;
    }
}
