#pragma once

#include <complex>

#include "util/memory_arena.hpp"

// Forward declarations
struct PFFFT_Setup;

/**
 * Harmonic/Percussive Source Separation (HPSS), based on
 * Harmonic/Percussive Separation Using Median Filtering, by Derry Fitzgerald, published at DAFx10.
 * (https://arrow.tudublin.ie/cgi/viewcontent.cgi?article=1078&context=argcon)
 *
 * There are two APIs that a user can use for processing data, the "Simple" API,
 * and the "Complex" API. The initialization/de-initialization methods
 * are shared by both APIs. The Simple API only contains one method: `process_window()`.
 * All other methods are part of the Complex API, starting with `push_new_window()`.
 */
namespace hpss
{
/**
 * HPSS Parameters
 * Many of these parameters are common for frequency-domain signal processing.
 * Please refer to the reference paper for more information about the kernel
 * size and mask power.
 */
struct Params
{
    int window_size = 1 << 11; // integer power of 2
    int hop_factor = 2; // integer power of 2 (right now only 1 and 2 are allowed)
    int zero_pad = 2; // integer power of 2
    int kernel_size = 17; // must be odd
    int mask_power = 2; // [0, 16] (usually either one or two)
};

// Forward declaration
struct Mediator;

/** HPSS computation state and pre-computed data. */
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

using complex = std::complex<float>;

/**
 * Creates an HPSS Processor for a given set of parameters.
 * All required memory for the HPSS processing is allocated here.
 */
HPSS_Processor init (Params params);

/** De-initializes an HPSS Processor and frees all allocate memory. */
void deinit (HPSS_Processor&);

/**
 * Generates a harmonic and percussive signal for a window of audio data.
 * The returned data is guaranteed to be valid until the next call to this
 * method.
 *
 * @param hop_data A span[proc.hop_size] of audio data.
 * @return A pair of span[proc.hop_size], with the order { harmonic, percussive }.
 */
std::pair<std::span<float>, std::span<float>> process_window (HPSS_Processor& proc, std::span<const float> hop_data);

/**
 * Pushes a new window of audio data into the processor.
 * Calling this method will reset the memory arena, thereby
 * invalidating the results of all the following methods.
 *
 * @param hop_data A span[proc.hop_size] of audio data.
 */
void push_new_window (HPSS_Processor& proc, std::span<const float> hop_data);

/**
 * Computes an FFT frame from the given window.
 *
 * @param window_data A span[proc.window_size] of audio data.
 * @return A span[proc.fft_size] of interleaved-complex FFT data.
 */
std::span<complex> process_forward_fft (HPSS_Processor& proc, std::span<const float> window_data);

/**
 * Computes the magnitudes of the FFT frame.
 * This will either be the absolute value or the squared value,
 * depending on the value of proc.mask_power.
 *
 * @param fft_frame A span[proc.fft_size] of interleaved-complex FFT data.
 * @return A span[proc.fft_size / 2 + 1] of magnitude data.
 */
std::span<float> compute_fft_magnitudes (HPSS_Processor& proc, std::span<const complex> fft_frame);

/**
 * Computes the inverse FFT of the given FFT frame.
 *
 * @param fft_data A span[proc.fft_size] of interleaved-complex FFT data.
 * @return A span[proc.window_size] of audio data.
 */
std::span<float> process_inverse_fft (HPSS_Processor& proc, std::span<const complex> fft_data);

/**
 * Generates a percussive (vertical) mask for some FFT magnitude data.
 * Note that unlike generate_harmonic_mask(), this method is "stateless".
 *
 * @param fft_abs_data A span[proc.fft_size / 2 + 1] of magnitude data.
 * @return A span[proc.fft_size / 2 + 1] of mask data.
 */
std::span<float> generate_percussive_mask (HPSS_Processor& proc, std::span<const float> fft_abs_data);

/**
 * Generates a harmonic (horizontal) mask for some FFT magnitude data.
 * Note that unlike generate_percussive_mask(), this method is "stateful".
 *
 * @param fft_abs_data A span[proc.fft_size / 2 + 1] of magnitude data.
 * @return A span[proc.fft_size / 2 + 1] of mask data.
 */
std::span<float> generate_harmonic_mask (HPSS_Processor& proc, std::span<const float> fft_abs_data);

/**
 * Balances the percussive and harmonic masks.
 *
 * @param percussive_mask A span[proc.fft_size / 2 + 1] of mask data.
 * @param harmonic_mask A span[proc.fft_size / 2 + 1] of mask data.
 */
void balance_masks (HPSS_Processor& proc, std::span<float> percussive_mask, std::span<float> harmonic_mask);

/**
 * Applies a spectral mask to an FFT frame.
 *
 * @param spectrum A span[proc.fft_size] of interleaved-complex FFT data.
 * @param mask A span[proc.fft_size / 2 + 1] of mask data.
 * @return A span[proc.fft_size] of interleaved-complex FFT data.
 */
std::span<complex> apply_spectral_mask (HPSS_Processor& proc, std::span<const complex> spectrum, std::span<const float> mask);

/**
 * Computes overlap-add processing for consecutive windows.
 *
 * @param window A span[proc.window_size] of audio data.
 * @param last_half_window A span[proc.window_size / 2] of audio data (e.g. proc.last_harmonic_half_window).
 * @return A span[proc.hop_size] of audio data.
 */
std::span<float> overlap_add (HPSS_Processor& proc, std::span<const float> window, std::span<float> last_half_window);
}
