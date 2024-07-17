#include "hpss.hpp"

#include <algorithm>
#include <complex>
#include <limits>

namespace hpss
{
// @TODO: maybe it would be faster to work with de-interleaved data?
using complex = std::complex<float>;

HPSS_Processor init (Params params)
{
    HPSS_Processor proc {};

    proc.window_size = params.window_size;
    proc.hop_size = params.window_size / params.hop_factor;
    proc.fft_size = params.window_size * params.zero_pad;
    proc.kernel_size = params.kernel_size;
    proc.mask_power = params.mask_power;

    proc.arena = new Memory_Arena<> {
        2 * proc.fft_size * sizeof (complex)
        + 2 * proc.window_size * sizeof (float)
        + (params.kernel_size + 2) * (proc.fft_size / 2 + 8) * sizeof (float)
        + 2 * params.kernel_size * sizeof (float)
        + 2048
    };

    proc.fft_setup = pffft_new_setup (proc.fft_size, PFFFT_COMPLEX);
    proc.fft_io_data = static_cast<float*> (pffft_aligned_malloc (2 * proc.fft_size * sizeof (complex)));

    proc.fft_history = proc.arena->make_span<std::span<float>> (proc.kernel_size);
    for (auto& fft_data : proc.fft_history)
    {
        fft_data = proc.arena->make_span<float> (proc.fft_size / 2 + 1, 32);
        std::fill (fft_data.begin(), fft_data.end(), 0.0f);
    }

    proc.arena_frame = proc.arena->create_frame();

    return proc;
}

void deinit (HPSS_Processor& proc)
{
    pffft_aligned_free (proc.fft_io_data);
    pffft_destroy_setup (proc.fft_setup);

    delete proc.arena;
}

static std::span<complex> process_forward_fft (HPSS_Processor& proc, std::span<const float> window_data)
{
    // zero-pad and convert real-> split-complex
    std::fill (proc.fft_io_data, proc.fft_io_data + 2 * proc.fft_size, 0.0f);
    for (size_t n = 0; n < window_data.size(); ++n)
        proc.fft_io_data[n * 2] = window_data[n];

    pffft_transform_ordered (proc.fft_setup,
                             proc.fft_io_data,
                             proc.fft_io_data,
                             nullptr, // optional "work buffer"
                             PFFFT_FORWARD);

    const auto fft_out = proc.arena->make_span<complex> (proc.fft_size, 32);
    std::copy (reinterpret_cast<complex*> (proc.fft_io_data),
               reinterpret_cast<complex*> (proc.fft_io_data) + proc.fft_size,
               fft_out.data());

    return fft_out;
}

static std::span<float> process_inverse_fft (HPSS_Processor& proc, std::span<const complex> fft_data)
{
    std::copy (fft_data.begin(), fft_data.end(), reinterpret_cast<complex*> (proc.fft_io_data));

    pffft_transform_ordered (proc.fft_setup,
                             proc.fft_io_data,
                             proc.fft_io_data,
                             nullptr, // optional "work buffer"
                             PFFFT_BACKWARD);

    const auto norm_gain = 1.0f / static_cast<float> (proc.fft_size);
    const auto ifft_out = proc.arena->make_span<float> (proc.window_size, 32);
    for (int n = 0; n < proc.window_size; ++n)
        ifft_out[n] = proc.fft_io_data[n * 2] * norm_gain;

    return ifft_out;
}

static std::span<complex> apply_spectral_mask (Memory_Arena<>& arena, std::span<const complex> spectrum, std::span<const float> mask)
{
    const auto spectrum_out = arena.make_span<complex> (spectrum.size(), 32);

    const auto N = spectrum.size();
    const auto M = mask.size();
    int n = 0;
    for (; n < M; ++n)
        spectrum_out[n] = spectrum[n] * mask[n];
    for (; n < N; ++n)
        spectrum_out[n] = spectrum[n] * mask[N - n];

    return spectrum_out;
}

std::pair<std::span<float>, std::span<float>> process_window (HPSS_Processor& proc, std::span<const float> window_data)
{
    proc.arena->reset_to_frame (proc.arena_frame);

    const auto fft_frame = process_forward_fft (proc, window_data);

    const auto fft_abs_data = proc.fft_history[proc.fft_history_index];
    proc.fft_history_index = (proc.fft_history_index + 1) % proc.kernel_size;

    fft_abs_data[0] = fft_frame[0].real();
    fft_abs_data[proc.fft_size / 2] = fft_frame[proc.fft_size / 2].real();
    for (int n = 1; n < proc.fft_size / 2; ++n)
        fft_abs_data[n] = std::abs (fft_frame[n]);

    const auto percussive_mask = proc.arena->make_span<float> (proc.fft_size / 2 + 1, 32);
    for (int n = 0; n < proc.fft_size / 2 + 1; ++n)
    {
        const auto start_index = std::max (n - proc.kernel_size / 2, 0);
        const auto end_index = std::min (n + proc.kernel_size / 2, proc.fft_size / 2);
        const auto median_count = end_index - start_index;
        const auto median_element = (median_count / 2) + 1;

        const auto _ = proc.arena->create_frame();
        const auto temp_data = proc.arena->make_span<float> (median_count);
        std::copy (fft_abs_data.begin() + start_index, fft_abs_data.begin() + start_index + median_count, temp_data.begin());

        std::nth_element (temp_data.begin(), temp_data.begin() + median_element, temp_data.end());
        percussive_mask[n] = temp_data[median_element];
    }

    const auto harmonic_mask = proc.arena->make_span<float> (proc.fft_size / 2 + 1, 32);
    for (int n = 0; n < proc.fft_size / 2 + 1; ++n)
    {
        const auto median_element = (proc.kernel_size / 2) + 1;

        const auto _ = proc.arena->create_frame();
        const auto temp_data = proc.arena->make_span<float> (proc.kernel_size);

        for (int k = 0; k < proc.kernel_size; ++k)
            temp_data[k] = proc.fft_history[k][n];

        std::nth_element (temp_data.begin(), temp_data.begin() + median_element, temp_data.end());
        harmonic_mask[n] = temp_data[median_element];
    }

    for (int n = 0; n < proc.fft_size / 2 + 1; ++n)
    {
        static constexpr auto eps = std::numeric_limits<float>::epsilon();

        const auto H_p = std::pow (harmonic_mask[n], proc.mask_power);
        const auto P_p = std::pow (percussive_mask[n], proc.mask_power);
        const auto denom = 1.0f / (H_p + P_p + eps);

        harmonic_mask[n] = H_p * denom;
        percussive_mask[n] = P_p * denom;
    }

    std::span<float> harmonic_out {};
    std::span<float> percussive_out {};

    {
        const auto _ = proc.arena->create_frame();
        const auto harmonic_spectrum = apply_spectral_mask (*proc.arena, fft_frame, harmonic_mask);
        harmonic_out = process_inverse_fft (proc, harmonic_spectrum);
    }

    {
        const auto _ = proc.arena->create_frame();
        const auto percussive_spectrum = apply_spectral_mask (*proc.arena, fft_frame, percussive_mask);
        percussive_out = process_inverse_fft (proc, percussive_spectrum);
    }

    return { harmonic_out, percussive_out };
}
}// namespace hpss

#if __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wvla-extension"
#endif

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4005)
#endif

#include "pffft/pffft.c" // NOLINT
#include "pffft/pffft_common.c" // NOLINT

#if __clang__
#pragma GCC diagnostic pop
#endif

#if defined(_MSC_VER)
#pragma warning(pop)
#endif
