#pragma once

#include <cstdint>
#include <list>

struct MemoryBlock {
    uint32_t bucket{0u};
    uint32_t row_count;
    uint64_t nnz_count;
};

using StackIterator = std::list<MemoryBlock>::iterator;
