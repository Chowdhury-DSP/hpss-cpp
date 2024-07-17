#pragma once

#include "util/memory_arena.hpp"

#if __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic ignored "-Wshadow-field-in-constructor"
#pragma GCC diagnostic ignored "-Wdeprecated-dynamic-exception-spec"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4100)
#endif

#include "pffft/pffft.h"

#if __clang__
#pragma GCC diagnostic pop
#endif

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

namespace hpss
{
struct Params
{
    int window_size = 1 << 11; // integer power of 2
    int hop_factor = 2; // integer power of 2
    int zero_pad = 2; // integer power of 2
    int kernel_size = 17; // must be odd
    float mask_power = 2.0f; // usually either one or two (maybe this can be an int?)
};

struct HPSS_Processor
{
    Memory_Arena<>* arena = nullptr;
    Memory_Arena<>::Frame arena_frame;

    int window_size = 0;
    int hop_size = 0;
    int fft_size = 0;
    int kernel_size = 0;
    float mask_power = 0.0f;

    PFFFT_Setup* fft_setup = nullptr;
    float* fft_io_data = nullptr;

    std::span<std::span<float>> fft_history;
    int fft_history_index = 0;
};

HPSS_Processor init (Params params);
void deinit (HPSS_Processor&);
std::pair<std::span<float>, std::span<float>> process_window (HPSS_Processor& proc, std::span<const float> window_data);
}
