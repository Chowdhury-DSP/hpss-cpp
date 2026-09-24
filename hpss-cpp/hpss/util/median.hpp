#pragma once

#include "memory_arena.hpp"

namespace hpss
{
struct Median
{
    using Idx_Type = int32_t;
    float* window {};
    Idx_Type* idxs {};
    Idx_Type ptr = {};
    Idx_Type window_size = {};

    void init (Memory_Arena& arena, int num)
    {
        window = arena.allocate<float> (num, 16);
        idxs = arena.allocate<Idx_Type> (num, 16);
        ptr = num / 2;
        window_size = num;

        init();
    }

    void init()
    {
        for (Idx_Type i = 0; i < window_size; ++i)
            idxs[i] = i;

        std::fill (window, window + window_size / 2, std::numeric_limits<float>::lowest());
        std::fill (window + window_size / 2, window + window_size, std::numeric_limits<float>::max());
    }

    // includes the worst-case padding for create()'s three allocations
    static size_t bytes_required (int num)
    {
        return sizeof (Median) + alignof (Median) + (size_t) num * (sizeof (float) + sizeof (Idx_Type)) + 2 * 16;
    }

    static Median* create (Memory_Arena& arena, int num)
    {
        auto* median = (Median*) arena.allocate_bytes (sizeof (Median), alignof (Median));
        median->init (arena, num);
        return median;
    }

    float push_and_return (float x)
    {
        auto i = std::distance (idxs, std::find (idxs, idxs + window_size, ptr));
        window[ptr] = x;

        while (i > 0 && window[idxs[i - 1]] > window[idxs[i]])
        {
            std::swap (idxs[i - 1], idxs[i]);
            i--;
        }

        while (i < window_size - 1 && window[idxs[i + 1]] < window[idxs[i]])
        {
            std::swap (idxs[i + 1], idxs[i]);
            i++;
        }

        ptr = (ptr == window_size - 1) ? 0 : ptr + 1;
        return window[idxs[window_size / 2]];
    }
};
}
