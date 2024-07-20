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

#if defined(__ARM_NEON__)
#define PFFFT_ENABLE_NEON
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
    int hop_factor = 2; // integer power of 2 (right now only 1 and 2 are allowed)
    int zero_pad = 2; // integer power of 2
    int kernel_size = 17; // must be odd
    int mask_power = 2; // [0, 16] (usually either one or two)
};

struct Mediator;
struct HPSS_Processor
{
    Memory_Arena<>* arena = nullptr;
    Memory_Arena<>::Frame arena_frame;

    int window_size = 0;
    int hop_size = 0;
    int fft_size = 0;
    int kernel_size = 0;
    int mask_power = 0;
    bool use_squares = false;

    PFFFT_Setup* fft_setup = nullptr;
    float* fft_io_data = nullptr;

    std::span<float> hann_window;
    std::span<float> window_in;
    std::span<float> last_half_window_harm;
    std::span<float> last_half_window_perc;

    std::span<Mediator*> horizontal_mediators;
};

HPSS_Processor init (Params params);
void deinit (HPSS_Processor&);
std::pair<std::span<float>, std::span<float>> process_window (HPSS_Processor& proc, std::span<const float> hop_data);
}
