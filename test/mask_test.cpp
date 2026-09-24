#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include <hpss/hpss.hpp>

static float naive_median (std::vector<float> values)
{
    const auto mid = values.begin() + (std::ptrdiff_t) values.size() / 2;
    std::nth_element (values.begin(), mid, values.end());
    return *mid;
}

// Sets up the harmonic medians in caller-owned memory (like Tenax does),
// filled with garbage so that anything reset_harmonic_medians() misses will show up.
static hpss::HPSS_Processor init_external (hpss::Memory_Arena& arena, int fft_size, int harmonic_kernel_size, int percussive_kernel_size)
{
    const auto num_bins = fft_size / 2 + 1;
    hpss::HPSS_Processor proc {};
    proc.fft_size = fft_size;
    proc.harmonic_kernel_size = harmonic_kernel_size;
    proc.percussive_kernel_size = percussive_kernel_size;
    proc.mask_power = 2;
    proc.arena = &arena;
    proc.median_windows = arena.make_span<float> (num_bins * harmonic_kernel_size);
    proc.median_idxs = arena.make_span<int32_t> (num_bins * harmonic_kernel_size);
    std::memset (proc.median_windows.data(), 0xAB, proc.median_windows.size_bytes());
    std::memset (proc.median_idxs.data(), 0xAB, proc.median_idxs.size_bytes());
    proc.median_ptr = 12345;
    hpss::reset_harmonic_medians (proc);
    proc.arena_frame = arena.create_frame();
    return proc;
}

static int check_masks (hpss::HPSS_Processor& proc, const char* name)
{
    const auto num_bins = proc.fft_size / 2 + 1;
    const auto half_kernel = proc.percussive_kernel_size / 2;

    std::mt19937 rng { 0x5eed };
    std::uniform_real_distribution<float> dist { 0.0f, 1.0f };

    // reference for the harmonic medians: one ring per bin, starting half-full of
    // lowest() and half-full of max(), with the first write in the middle
    std::vector<std::vector<float>> rings ((size_t) num_bins);
    for (auto& ring : rings)
    {
        ring.resize ((size_t) proc.harmonic_kernel_size, std::numeric_limits<float>::max());
        std::fill_n (ring.begin(), proc.harmonic_kernel_size / 2, std::numeric_limits<float>::lowest());
    }
    auto ring_ptr = (size_t) proc.harmonic_kernel_size / 2;

    int failures = 0;
    for (int frame = 0; frame < 4 * proc.harmonic_kernel_size; ++frame)
    {
        std::vector<float> mags ((size_t) num_bins);
        for (auto& x : mags)
            x = dist (rng);
        for (size_t n = 0; n < rings.size(); ++n)
            rings[n][ring_ptr] = mags[n];
        ring_ptr = (ring_ptr + 1) % (size_t) proc.harmonic_kernel_size;

        proc.arena->reset_to_frame (proc.arena_frame);
        const auto percussive = hpss::generate_percussive_mask (proc, mags);
        const auto harmonic = hpss::generate_harmonic_mask (proc, mags);

        // percussive: median across the neighbouring bins, truncated at the edges
        for (int n = 0; n < num_bins; ++n)
        {
            const auto start = std::max (n - half_kernel, 0);
            const auto end = std::min (n + half_kernel + 1, num_bins);
            failures += percussive[n] != naive_median ({ mags.begin() + start, mags.begin() + end });
        }

        // harmonic: median across the last harmonic_kernel_size frames
        for (int n = 0; n < num_bins; ++n)
            failures += harmonic[n] != naive_median (rings[(size_t) n]);

        hpss::balance_masks (proc, percussive, harmonic);
        for (int n = 0; n < num_bins; ++n)
        {
            if (! std::isfinite (percussive[n]) || ! std::isfinite (harmonic[n]))
            {
                std::cout << name << ": non-finite mask at frame " << frame << ", bin " << n << "\n";
                return 1;
            }
        }
    }

    std::cout << name << ": " << (failures == 0 ? "passed" : "FAILED") << " (" << failures << " mismatches)\n";
    return failures > 0;
}

int main()
{
    int result = 0;
    for (auto [harmonic_kernel_size, percussive_kernel_size] : { std::pair { 17, 17 }, { 33, 9 }, { 5, 65 } })
    {
        for (int fft_size : { 256, 2048 })
        {
            auto proc = hpss::init ({
                .window_size = fft_size,
                .zero_pad = 1,
                .harmonic_kernel_size = harmonic_kernel_size,
                .percussive_kernel_size = percussive_kernel_size,
            });
            result |= check_masks (proc, "hpss::init()");
            hpss::deinit (proc);

            hpss::Memory_Arena arena { 1 << 21 };
            auto external_proc = init_external (arena, fft_size, harmonic_kernel_size, percussive_kernel_size);
            result |= check_masks (external_proc, "reset_harmonic_medians()");
        }
    }
    return result;
}
