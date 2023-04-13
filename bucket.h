#pragma once

#include "memoryblock.h"

#include <limits>
#include <vector>

class Bucket
{
public:
    using min_type      = unsigned long;
    using count_type    = unsigned long;
    using StackIterator = std::list<MemoryBlock>::iterator;

    count_type    access_count = count_type{0};
    StackIterator marker;

    static inline std::vector<min_type> min_dists; // minimum reuse distance in bucket
    static inline constexpr min_type    inf_dist{std::numeric_limits<min_type>::max()};
};
