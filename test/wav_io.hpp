#pragma once

#pragma once

#include <hpss/util/memory_arena.hpp>
#include <iostream>
#include <sndfile.h>
#include <span>

/**
 * Utility functions for loading .wav files
 * into 2D vectors using libsndfile.
 */
namespace wav_io
{
using Buffer = std::span<std::span<float>>;
inline Buffer make_buffer (hpss::Memory_Arena& arena, int num_channels, int num_samples)
{
    auto* data = arena.allocate<float> (num_channels * num_samples);
    auto buffer = arena.make_span<std::span<float>> (num_channels);
    for (int channel_index = 0; auto& channel_data : buffer)
        channel_data = std::span { data + (channel_index++) * num_samples, (size_t) num_samples };
    return buffer;
}

using SND_PTR = std::unique_ptr<SNDFILE, decltype (&sf_close)>;

inline Buffer load_file (const char* file, SF_INFO& sf_info, hpss::Memory_Arena& arena)
{
    std::cout << "Loading file: " << file << std::endl;

    SND_PTR wav_file { sf_open (file, SFM_READ, &sf_info), &sf_close };

    if (sf_info.frames == 0)
    {
        std::cout << "File could not be opened!" << std::endl;
        std::exit (1);
    }

    const auto buffer = make_buffer (arena, sf_info.channels, (int) sf_info.frames);

    const auto frame = arena.create_frame();
    auto* interleaved_data = arena.allocate<float> (sf_info.channels * sf_info.frames);
    sf_readf_float (wav_file.get(), interleaved_data, sf_info.frames);

    // de-interleave channels
    for (int i = 0; i < sf_info.frames; ++i)
    {
        int interleaved_ptr = i * sf_info.channels;
        for (size_t ch = 0; ch < sf_info.channels; ++ch)
            buffer[ch][i] = interleaved_data[interleaved_ptr + ch];
    }

    return buffer;
}

inline void normalize (const Buffer& buffer)
{
    float max = 0.0f;
    for (const auto& channel : buffer)
    {
        for (auto x : channel)
            max = std::max (std::abs (x), max);
    }

    if (max > 1.0f)
    {
        const auto gain = 1.0f / max;
        for (const auto& channel : buffer)
        {
            for (auto& x : channel)
                x *= gain;
        }
    }
}

inline void write_file (const char* file, const Buffer& audio, SF_INFO& sf_info, hpss::Memory_Arena& arena)
{
    std::cout << "Writing to file: " << file << std::endl;

    normalize (audio);

    const auto channels = (int) audio.size();
    const auto frames = (sf_count_t) audio[0].size();
    sf_info.frames = frames;

    SND_PTR wav_file { sf_open (file, SFM_WRITE, &sf_info), &sf_close };

    const auto frame = arena.create_frame();

    auto* interleaved_data = arena.allocate<float> (channels * frames);

    // de-interleave channels
    for (int i = 0; i < frames; ++i)
    {
        int interleaved_ptr = i * channels;
        for (int ch = 0; ch < channels; ++ch)
            interleaved_data[interleaved_ptr + ch] = audio[ch][i];
    }

    sf_writef_float (wav_file.get(), interleaved_data, frames);
}

} // namespace wav_io
