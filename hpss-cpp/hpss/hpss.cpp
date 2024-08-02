#include "hpss.hpp"

#include <algorithm>
#include <complex>
#include <limits>
#include <numeric>

#include "util/mediator.hpp"
#include "util/power.hpp"
#include "util/pffft_wrapper.hpp"

namespace hpss
{
HPSS_Processor init (Params params)
{
    HPSS_Processor proc {};

    proc.window_size = params.window_size;
    proc.hop_size = params.window_size / params.hop_factor;
    proc.fft_size = params.window_size * params.zero_pad;
    proc.kernel_size = params.kernel_size;
    proc.mask_power = params.mask_power;
    proc.use_squares = proc.mask_power > 2 && proc.mask_power % 2 == 0;

    const auto mediator_size = MediatorSizeBytes (proc.kernel_size);
    proc.arena = new Memory_Arena<> {
        1 * proc.fft_size * sizeof (complex)
        + 3 * proc.window_size * sizeof (float)
        + 2 * (proc.window_size / 2) * sizeof (float)
        + 2 * proc.hop_size * sizeof (float)
        + 5 * (proc.fft_size / 2 + 8) * sizeof (float)
        + (proc.fft_size / 2 + 1) * (mediator_size + 16)
        + 2048
    };

    proc.fft_setup = pffft_new_setup (proc.fft_size, PFFFT_REAL);
    proc.fft_io_data = static_cast<float*> (pffft_aligned_malloc (2 * proc.fft_size * sizeof (complex)));

    proc.hann_window = proc.arena->make_span<float> (proc.window_size, 32);
    for (int n = 0; n < proc.window_size; ++n)
    {
        const auto sine = std::sin ((float) n * (float) M_PI / (float) proc.window_size);
        proc.hann_window[n] = sine * sine;
    }

    proc.window_in = proc.arena->make_span<float> (proc.window_size, 32);
    std::fill (proc.window_in.begin(), proc.window_in.end(), 0.0f);
    proc.last_half_window_harm = proc.arena->make_span<float> (proc.window_size / 2, 32);
    std::fill (proc.last_half_window_harm.begin(), proc.last_half_window_harm.end(), 0.0f);
    proc.last_half_window_perc = proc.arena->make_span<float> (proc.window_size / 2, 32);
    std::fill (proc.last_half_window_harm.begin(), proc.last_half_window_harm.end(), 0.0f);

    proc.horizontal_mediators = proc.arena->make_span<Mediator*> (proc.fft_size / 2 + 1);
    for (auto& mediator : proc.horizontal_mediators)
        mediator = MediatorNew (*proc.arena, proc.kernel_size);

    proc.arena_frame = proc.arena->create_frame();

    return proc;
}

void deinit (HPSS_Processor& proc)
{
    pffft_aligned_free (proc.fft_io_data);
    pffft_destroy_setup (proc.fft_setup);

    delete proc.arena;
}

std::span<complex> process_forward_fft (HPSS_Processor& proc, std::span<const float> window_data)
{
    // zero-pad
    std::copy (window_data.begin(), window_data.end(), proc.fft_io_data);
    std::fill (proc.fft_io_data + proc.window_size, proc.fft_io_data + proc.fft_size, 0.0f);

    pffft_transform_ordered (proc.fft_setup,
                             proc.fft_io_data,
                             proc.fft_io_data,
                             nullptr, // optional "work buffer"
                             PFFFT_FORWARD);

    const auto fft_out = proc.arena->make_span<complex> (proc.fft_size / 2 + 1, 32);
    std::copy (reinterpret_cast<complex*> (proc.fft_io_data),
               reinterpret_cast<complex*> (proc.fft_io_data) + proc.fft_size,
               fft_out.data());
    fft_out.back() = fft_out.front().imag();
    fft_out.front() = fft_out.front().real();

    return fft_out;
}

std::span<float> compute_fft_magnitudes (HPSS_Processor& proc, std::span<const complex> fft_frame)
{
    const auto fft_abs_data = proc.arena->make_span<float> (proc.fft_size / 2 + 1, 32);

    // If mask_power is a multiple of 2, we compute the squared FFT values instead of absolute FFT values,
    // which saves us a lot of std::sqrt calls. Then later we can use half the mask power when combining the
    // masks, which saves us a bunch of multiplies!
    if (proc.use_squares)
    {
        fft_abs_data[0] = power::ipow<2> (fft_frame[0].real());
        fft_abs_data[proc.fft_size / 2] = power::ipow<2> (fft_frame[proc.fft_size / 2].real());
        for (int n = 1; n < proc.fft_size / 2; ++n)
            fft_abs_data[n] = power::ipow<2> (fft_frame[n].real()) + power::ipow<2> (fft_frame[n].imag());
    }
    else
    {
        fft_abs_data[0] = fft_frame[0].real();
        fft_abs_data[proc.fft_size / 2] = fft_frame[proc.fft_size / 2].real();
        for (int n = 1; n < proc.fft_size / 2; ++n)
            fft_abs_data[n] = std::sqrt (power::ipow<2> (fft_frame[n].real()) + power::ipow<2> (fft_frame[n].imag()));
    }

    return fft_abs_data;
}

std::span<float> process_inverse_fft (HPSS_Processor& proc, std::span<const complex> fft_data)
{
    std::copy (fft_data.begin(), fft_data.end(), reinterpret_cast<complex*> (proc.fft_io_data));
    proc.fft_io_data[1] = proc.fft_io_data[proc.fft_size / 2];
    proc.fft_io_data[proc.fft_size / 2] = 0.0f;

    pffft_transform_ordered (proc.fft_setup,
                             proc.fft_io_data,
                             proc.fft_io_data,
                             nullptr, // optional "work buffer"
                             PFFFT_BACKWARD);

    const auto norm_gain = 1.0f / static_cast<float> (proc.fft_size);
    const auto ifft_out = proc.arena->make_span<float> (proc.window_size, 32);
    for (int n = 0; n < proc.window_size; ++n)
        ifft_out[n] = proc.fft_io_data[n] * norm_gain;

    return ifft_out;
}

std::span<float> generate_percussive_mask (HPSS_Processor& proc, std::span<const float> fft_abs_data)
{
    const auto half_kernel = proc.kernel_size / 2;
    const auto percussive_mask = proc.arena->make_span<float> (proc.fft_size / 2 + 1, 32);

#if 1
    const auto frame = proc.arena->create_frame();
    auto* mediator = MediatorNew (*proc.arena, proc.kernel_size);

    for (int n = 0; n < half_kernel; ++n)
    {
        MediatorInsert (mediator, fft_abs_data[n]);
    }

    int n;
    for (n = 0; n < proc.fft_size / 2 + 1 - half_kernel; ++n)
    {
        MediatorInsert (mediator, fft_abs_data[n + half_kernel]);
        percussive_mask[n] = MediatorMedian (mediator);
    }

    for (; n < proc.fft_size / 2 + 1; ++n)
    {
        MediatorInsert (mediator, n % 2 ? 0.0f : 10000.0f);
        percussive_mask[n] = MediatorMedian (mediator);
    }
#else
    for (int n = 0; n < proc.fft_size / 2 + 1; ++n)
    {
        const auto start_index = std::max (n - half_kernel, 0);
        const auto end_index = std::min (n + half_kernel, proc.fft_size / 2) + 1;
        const auto median_count = end_index - start_index;
        const auto median_element = median_count / 2;

        const auto _ = proc.arena->create_frame();
        const auto temp_data = proc.arena->make_span<float> (median_count);
        std::copy (fft_abs_data.begin() + start_index, fft_abs_data.begin() + start_index + median_count, temp_data.begin());

        std::nth_element (temp_data.begin(), temp_data.begin() + median_element, temp_data.end());
        percussive_mask[n] = temp_data[median_element];
    }
#endif
    return percussive_mask;
}

std::span<float> generate_harmonic_mask (HPSS_Processor& proc, std::span<const float> fft_abs_data)
{
    const auto harmonic_mask = proc.arena->make_span<float> (proc.fft_size / 2 + 1, 32);
    for (int n = 0; n < proc.fft_size / 2 + 1; ++n)
    {
#if 1
        MediatorInsert (proc.horizontal_mediators[n], fft_abs_data[n]);
        harmonic_mask[n] = MediatorMedian (proc.horizontal_mediators[n]);
#else
        const auto median_element = (proc.kernel_size / 2) + 1;

        const auto _ = proc.arena->create_frame();
        const auto temp_data = proc.arena->make_span<float> (proc.kernel_size);

        for (int k = 0; k < proc.kernel_size; ++k)
            temp_data[k] = proc.fft_history[k][n];

        std::nth_element (temp_data.begin(), temp_data.begin() + median_element, temp_data.end());
        harmonic_mask[n] = temp_data[median_element];
#endif
    }
    return harmonic_mask;
}

static void apply_power (int exp, std::span<float> data)
{
#define HPSS_POWER_EXP(exp_val) \
    case (exp_val): \
    for (auto& x : data) \
        x = power::ipow<(exp_val)> (x); \
    return

    switch (exp)
    {
        case 0:
            std::fill (data.begin(), data.end(), 1.0f);
            return;
        case 1:
            return;
        HPSS_POWER_EXP (2);
        HPSS_POWER_EXP (3);
        HPSS_POWER_EXP (4);
        HPSS_POWER_EXP (5);
        HPSS_POWER_EXP (6);
        HPSS_POWER_EXP (7);
        HPSS_POWER_EXP (8);
        HPSS_POWER_EXP (9);
        HPSS_POWER_EXP (10);
        HPSS_POWER_EXP (11);
        HPSS_POWER_EXP (12);
        HPSS_POWER_EXP (13);
        HPSS_POWER_EXP (14);
        HPSS_POWER_EXP (15);
        HPSS_POWER_EXP (16);
        default:
            return;
    }
}

void balance_masks (HPSS_Processor& proc, std::span<float> percussive_mask, std::span<float> harmonic_mask)
{
    // See comment in compute_fft_magnitudes().
    const auto effective_mask_power = proc.use_squares ? proc.mask_power / 2 : proc.mask_power;
    apply_power (effective_mask_power, percussive_mask);
    apply_power (effective_mask_power, harmonic_mask);

    for (int n = 0; n < (int) percussive_mask.size(); ++n)
    {
        static constexpr auto eps = std::numeric_limits<float>::epsilon();

        const auto H_p = harmonic_mask[n];
        const auto P_p = percussive_mask[n];
        const auto denom = 1.0f / (H_p + P_p + eps);

        harmonic_mask[n] = H_p * denom;
        percussive_mask[n] = P_p * denom;
    }
}

std::span<complex> apply_spectral_mask (HPSS_Processor& proc, std::span<const complex> spectrum, std::span<const float> mask)
{
    const auto spectrum_out = proc.arena->make_span<complex> (spectrum.size(), 32);

    const auto M = mask.size();
    for (int n = 0; n < M; ++n)
        spectrum_out[n] = spectrum[n] * mask[n];

    return spectrum_out;
}

std::span<float> overlap_add (HPSS_Processor& proc, std::span<const float> window, std::span<float> last_half_window)
{
    const auto hop_out = proc.arena->make_span<float> (proc.hop_size, 32);
    if (proc.hop_size == proc.window_size)
    {
        std::copy (window.begin(), window.end(), hop_out.begin());
    }
    else if (proc.hop_size == proc.window_size / 2)
    {
        for (int n = 0; n < proc.hop_size; ++n)
        {
            hop_out[n] = window[n] * proc.hann_window[n];
            hop_out[n] += last_half_window[n] * proc.hann_window[n + proc.hop_size];
        }
        std::copy (window.begin() + proc.hop_size, window.end(), last_half_window.begin());
    }
    else
    {
        assert (false); // @TODO!
    }

    return hop_out;
}

void push_new_window (HPSS_Processor& proc, std::span<const float> hop_data)
{
    proc.arena->reset_to_frame (proc.arena_frame);

    if (proc.hop_size == proc.window_size)
    {
        std::copy (hop_data.begin(), hop_data.end(), proc.window_in.begin());
    }
    else
    {
        std::copy (proc.window_in.begin() + proc.hop_size, proc.window_in.end(), proc.window_in.begin());
        std::copy (hop_data.begin(), hop_data.end(), proc.window_in.begin() + (proc.window_size - proc.hop_size));
    }
}

std::pair<std::span<float>, std::span<float>> process_window (HPSS_Processor& proc, std::span<const float> hop_data)
{
    push_new_window (proc, hop_data);

    const auto fft_frame = process_forward_fft (proc, proc.window_in);
    const auto fft_abs_data = compute_fft_magnitudes (proc, fft_frame);

    const auto percussive_mask = generate_percussive_mask (proc, fft_abs_data);
    const auto harmonic_mask = generate_harmonic_mask (proc, fft_abs_data);
    balance_masks (proc, percussive_mask, harmonic_mask);

    const auto harmonic_spectrum = apply_spectral_mask (proc, fft_frame, harmonic_mask);
    const auto harmonic_out = process_inverse_fft (proc, harmonic_spectrum);

    const auto percussive_spectrum = apply_spectral_mask (proc, fft_frame, percussive_mask);
    const auto percussive_out = process_inverse_fft (proc, percussive_spectrum);

    const auto hop_harmonic = overlap_add (proc, harmonic_out, proc.last_half_window_harm);
    const auto hop_percussive = overlap_add (proc, percussive_out, proc.last_half_window_perc);

    return { hop_harmonic, hop_percussive };
}
} // namespace hpss

#include "util/mediator.cpp"
