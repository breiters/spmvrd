#pragma once

#include <list>

struct MemoryBlock {
    unsigned bucket{0u};
};

using StackIterator = std::list<MemoryBlock>::iterator;
