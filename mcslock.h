#pragma once

#include "config.h"

// C++ atomic library not implemented in PIN C++ STL
// ==> we use C atomics instead
#include <atomic>

#include <cstddef>
#include <cstdint>

// using std::atomic;
using std::atomic_bool;
using std::atomic_exchange;
using std::atomic_store;
using std::atomic_uintptr_t;

struct alignas(CACHE_LINESIZE) qnode {
    atomic_uintptr_t next;
    atomic_bool      wait;
};

class MCSLock
{
public:
    MCSLock() { atomic_store(&_tail, (uintptr_t)NULL); };
    ~MCSLock(){};

    void lock(int tid) { lock(&_nodes[tid]); }

    void unlock(int tid) { unlock(&_nodes[tid]); }

    void lock(qnode *p)
    {
        atomic_store(&(p->next), (uintptr_t)NULL);
        atomic_store(&(p->wait), true);

        // exchange _tail with p, return old value of _tail
        qnode *prev = (qnode *)atomic_exchange(&_tail, (uintptr_t)p);

        if (prev != NULL) {
            atomic_store(&(prev->next), (uintptr_t)p);
            while (atomic_load(&(p->wait))) {
                /* spin */
            }
        }
    }

    // lock can be acquired only when p is head of list
    // => release is only called when p is head of list
    void unlock(qnode *p)
    {
        qnode *succ = (qnode *)atomic_load(&(p->next));

        if (succ == NULL) {
            qnode *old_p = p; // need to keep p
            if (atomic_compare_exchange_strong(&_tail, (uintptr_t *)&old_p, (uintptr_t)NULL)) {
                // _tail was p and is now nullptr
                return;
            }
            // another thread changed _tail (new qnode enters queue)
            // wait until other thread has appended his qnode!
            do {
                succ = (qnode *)atomic_load(&(p->next));
            } while (succ == NULL);
        }
        atomic_store(&(succ->wait), false);
    }

private:
    atomic_uintptr_t _tail;
    qnode            _nodes[MAX_THREADS];
};
