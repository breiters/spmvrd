#pragma once

// cache line size on the machine used for profiling
#ifndef CACHE_LINESIZE
#    define CACHE_LINESIZE 256
#endif

// cache line size of the target architecture (must be a power-of-two)
#define MEMBLOCKLEN 256

// Consistency checks?
#define DEBUG       0

// 2: Huge amount of debug output, 1: checks, 0: silent
#define VERBOSE     0

// Assertions and consistency check?

#if NDEBUG
#    define RD_DEBUG 0
#else
#    define RD_DEBUG 2
#endif

// 2: Huge amount of debug output, 1: checks, 0: silent
#define RD_VERBOSE                    0

// Print only up to N-th marker in debug output
#define RD_PRINT_MARKER_MAX           2

// Do not account for zero distance accesses?
#define RD_DO_NOT_COUNT_ZERO_DISTANCE 0

#define MAX_THREADS                   48

#if RD_VERBOSE
#    define eprintf(...) fprintf(stderr, __VA_ARGS__)
#else
#    define eprintf(...)
#endif
