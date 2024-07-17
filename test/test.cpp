#include "hpss/hpss.hpp"
#include "wav_io.hpp"

void help()
{
    std::cout << "Utility to separate harmonic and percussive signals from a .wav file" << std::endl;
    std::cout << "Usage: hpss <wav_file> [<num_seconds> <start_seconds>]" << std::endl;
}

int main (int argc, char* argv[])
{
    if (argc < 2 || argc > 4)
    {
        help();
        return 1;
    }

    if (argc == 2 && std::string (argv[1]) == "--help")
    {
        help();
        return 1;
    }

    std::string test_file = std::string (argv[1]);

    float num_seconds = 10.0f;
    if (argc >= 3)
        num_seconds = (float) std::atof (argv[2]);

    float start_seconds = 0.0f;
    if (argc >= 4)
        start_seconds = (float) std::atof (argv[3]);

    hpss::Memory_Arena<> arena { 1 << 29 };
    SF_INFO sf_info;
    auto ref_signal = wav_io::load_file (test_file.c_str(), sf_info, arena);
    const auto fs = sf_info.samplerate;

    // trim signal
    const auto num_channels = (int) ref_signal.size();
    const auto start_sample = int ((float) fs * start_seconds);
    const auto num_samples = std::min (int ((float) fs * num_seconds), (int) sf_info.frames - start_sample);
    for (auto& channel : ref_signal)
        channel = channel.subspan (start_sample, num_samples);

    const auto harmonic_signal = wav_io::make_buffer (arena, num_channels, num_samples);
    const auto percussive_signal = wav_io::make_buffer (arena, num_channels, num_samples);
    const auto sum_signal = wav_io::make_buffer (arena, num_channels, num_samples);

    const hpss::Params params
    {
        .window_size = 1 << 12,
        .hop_factor = 2,
        .zero_pad = 2,
    };
    const auto hpss_procs = arena.make_span<hpss::HPSS_Processor> (num_channels);
    for (auto& proc : hpss_procs)
        proc = hpss::init (params);

    int sample_count = 0;
    while (sample_count + hpss_procs[0].hop_size <= num_samples)
    {
        for (int channel = 0; channel < num_channels; ++channel)
        {
            const auto window = ref_signal[channel].subspan (sample_count, hpss_procs[channel].hop_size);

            const auto [harmonic_out, percussive_out] = hpss::process_window (hpss_procs[channel], window);

            std::copy (harmonic_out.begin(), harmonic_out.end(), harmonic_signal[channel].begin() + sample_count);
            std::copy (percussive_out.begin(), percussive_out.end(), percussive_signal[channel].begin() + sample_count);

            for (int n = sample_count; n < sample_count + hpss_procs[0].hop_size; ++n)
                sum_signal[channel][n] = harmonic_signal[channel][n] + percussive_signal[channel][n];
        }

        sample_count += hpss_procs[0].hop_size;
    }

    for (auto& proc : hpss_procs)
        hpss::deinit (proc);

    wav_io::write_file ("ref.wav", ref_signal, sf_info, arena);
    wav_io::write_file("harmonic.wav", harmonic_signal, sf_info, arena);
    wav_io::write_file("percussive.wav", percussive_signal, sf_info, arena);
    wav_io::write_file("sum.wav", sum_signal, sf_info, arena);

    return 0;
}
