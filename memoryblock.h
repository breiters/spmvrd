#pragma once

#include "config.h"

#include <list>

#if REUSE_DISTANCE_METHOD_NEW
#include <cinttypes>
#endif
struct MemoryBlock {
#if REUSE_DISTANCE_METHOD_NEW
    unsigned reference_time_{0u}; // possibly needs to be uint64_t for very large matrices
#endif
    unsigned bucket{0u};
};

using StackIterator = std::list<MemoryBlock>::iterator;
