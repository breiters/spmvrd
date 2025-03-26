#pragma once

#include "memoryblock.h"

#include <limits>
#include <vector>

struct Bucket {
    using StackIterator = std::list<MemoryBlock>::iterator;
    using min_type      = unsigned long;

    using count_type = unsigned long;
    // count_type    access_count = count_type{0};

    Bucket(StackIterator it) : marker{it} {}
    
    struct Counts {
        count_type count_x;
        count_type count_xy;
        count_type count_xya;
    };

    Counts        access_counts = Counts{0u, 0u, 0u};
    StackIterator marker;

    static inline std::vector<min_type> min_dists; // minimum reuse distance in bucket
    static inline constexpr min_type    INF_DIST{std::numeric_limits<min_type>::max()};
};
