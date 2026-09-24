#include <chrono>
#include <iostream>
#include <random>

#include <hpss/util/memory_arena.hpp>
#include <hpss/util/median.hpp>

#include "heap_mediator/mediator.hpp"
#include "heap_mediator/mediator.cpp"

static constexpr int N = 17;
static constexpr size_t M = 10;

template <typename Container>
static void print_span (const Container& buffer)
{
    std::cout << "{";
    for (size_t i = 0; i < buffer.size(); ++i)
    {
        if (i != 0)
            std::cout << ", ";
        std::cout << buffer[i];
    }
    std::cout << "}\n";
}

struct Median_Naive
{
    std::vector<float> window {};
    std::vector<float> window_copy {};
    size_t ptr = 0;

    Median_Naive()
    {
        window.reserve (N);
        window_copy.reserve (N);
    }

    float push_and_return (float x)
    {
        if (window.size() < N)
        {
            window.push_back (x);
            window_copy.resize (window.size());
        }
        else
        {
            window[ptr] = x;
            ptr++;
            if (ptr == N)
                ptr = 0;
        }

        std::copy (window.begin(), window.end(), window_copy.begin());
        std::nth_element (window_copy.begin(), window_copy.begin() + window.size() / 2, window_copy.end());
        return window_copy[window.size() / 2];
    }
};

struct Median_New
{
    using Idx_Type = int32_t;
    void* data {};
    // float* window {};
    // Idx_Type* idxs {};
    Idx_Type ptr = {};
    Idx_Type window_size = {};

    Median_New (hpss::Memory_Arena& arena, int num)
    {
        data = arena.allocate_bytes (4 * 2 * num, 8); // <float> (num)
        // window = arena.allocate<float> (num);
        // idxs = arena.allocate<Idx_Type> (num);
        ptr = num / 2;
        window_size = num;

        auto* window = (float*) data;
        auto* idxs = (Idx_Type*) data + window_size;
        std::iota (idxs, idxs + window_size, 0);
        std::fill (window, window + num / 2, std::numeric_limits<float>::lowest());
        std::fill (window + num / 2, window + window_size, std::numeric_limits<float>::max());
    }

    inline static void swap (Idx_Type* x, Idx_Type* y)
    {
        // trying out different swap methods
        // in theory the xor swap should be slower?
        // in measurement it seems about the same?

        std::swap (*x, *y);

        // auto tmp = *x;
        // *x = *y;
        // *y = tmp;

        // *x ^= *y;
        // *y ^= *x;
        // *x ^= *y;

        // static constexpr int32_t m = 0xFF;
        // *y |= (*x << 8);
        // *x = *y & m;
        // *y >>= 8;
    }

    float push_and_return (float x)
    {
        auto* window = (float*) data;
        auto* idxs = (Idx_Type*) data + window_size;

        auto i = std::distance (idxs, std::find (idxs, idxs + window_size, ptr));
        window[ptr] = x;

        assert (idxs[i] == ptr);
        while (i > 0 && window[idxs[i - 1]] > window[idxs[i]])
        {
            swap (idxs + i - 1, idxs + i);
            // std::swap (idxs[i - 1], idxs[i]);
            i--;
        }

        while (i < window_size - 1 && window[idxs[i + 1]] < window[idxs[i]])
        {
            swap (idxs + i + 1, idxs + i);
            // std::swap (idxs[i + 1], idxs[i]);
            i++;
        }

        ptr = (ptr == window_size - 1) ? 0 : ptr + 1;

        return window[idxs[window_size / 2]];
    }

    void print_sorted (float x) const
    {
        auto* window = (float*) data;
        auto* idxs = (Idx_Type*) data + window_size;

        std::cout << x << " " << ptr << " {";
        for (size_t i = 0; i < N; ++i)
        {
            if (i != 0)
                std::cout << ", ";
            std::cout << window[idxs[i]];
        }
        std::cout << "}\n";
    }
};

int main()
{
    std::cout << "Median test" << std::endl;

    hpss::Memory_Arena arena { 4096 };
    {
        auto* mediator = hpss::MediatorNew (arena, N);
        Median_Naive median_naive {};
        hpss::Median median {};
        median.init (arena, N);

        std::array<float, M> data { 1.0f, 3.0f, 2.0f, -1.0f, 10.0f, 9.0f, 18.0f, -20.0f, 0.0f, 1.0f };
        std::array<float, M> ref {};
        std::array<float, M> naive {};
        std::array<float, M> mine {};

        for (size_t i = 0; i < M; ++i)
        {
            hpss::MediatorInsert (mediator, data[i]);
            ref[i] = hpss::MediatorMedian (mediator);
            naive[i] = median_naive.push_and_return (data[i]);
            mine[i] = median.push_and_return (data[i]);
        }

        print_span (ref);
        print_span (naive);
        print_span (mine);
    }

    arena.clear();

    {
        auto* mediator = hpss::MediatorNew (arena, N);
        Median_Naive median_naive {};
        hpss::Median median {};
        median.init (arena, N);
        Median_New median_new { arena, N };

        std::vector<float> data {};
        data.resize (1'000'000);
        std::random_device rd {};
        std::mt19937 rand_eng { rd() };
        std::uniform_real_distribution<float> dist { -10.0f, 10.0f };
        for (auto& x : data)
            x = dist (rand_eng);

        float dummy {};

        auto start = std::chrono::high_resolution_clock::now();
        for (auto x : data)
        {
            hpss::MediatorInsert (mediator, x);
            dummy = hpss::MediatorMedian (mediator);
        }
        auto duration = std::chrono::high_resolution_clock::now() - start;
        auto duration_seconds = std::chrono::duration<float> (duration).count();
        std::cout << "Heap Median: " << duration_seconds << " seconds" << std::endl;
        std::cout << dummy << std::endl;

        // start = std::chrono::high_resolution_clock::now();
        // for (auto x : data)
        // {
        //     dummy = median_naive.push_and_return (x);
        // }
        // duration = std::chrono::high_resolution_clock::now() - start;
        // duration_seconds = std::chrono::duration<float> (duration).count();
        // std::cout << "Naive Median: " << duration_seconds << " seconds" << std::endl;
        // std::cout << dummy << std::endl;

        start = std::chrono::high_resolution_clock::now();
        for (auto x : data)
        {
            dummy = median.push_and_return (x);
        }
        duration = std::chrono::high_resolution_clock::now() - start;
        duration_seconds = std::chrono::duration<float> (duration).count();
        std::cout << "My Median: " << duration_seconds << " seconds" << std::endl;
        std::cout << dummy << std::endl;

        start = std::chrono::high_resolution_clock::now();
        for (auto x : data)
        {
            dummy = median_new.push_and_return (x);
        }
        duration = std::chrono::high_resolution_clock::now() - start;
        duration_seconds = std::chrono::duration<float> (duration).count();
        std::cout << "My Median New: " << duration_seconds << " seconds" << std::endl;
        std::cout << dummy << std::endl;
    }

    return 0;
}
