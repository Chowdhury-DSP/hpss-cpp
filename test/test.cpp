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
    while (sample_count + params.window_size <= num_samples)
    {
        for (int channel = 0; channel < num_channels; ++channel)
        {
            const auto window = ref_signal[channel].subspan (sample_count, params.window_size);

            hpss::process_window (hpss_procs[channel], window);
        }

        sample_count += hpss_procs[0].hop_size;
    }

    for (auto& proc : hpss_procs)
        hpss::deinit (proc);

    // HPSS::HPSS_PARAMS params;
    // params.sample_rate = fs;
    // params.debug = true;
    //
    // std::vector<std::vector<float>> h_signal;
    // std::vector<std::vector<float>> p_signal;
    // std::vector<std::vector<float>> sum_signal;
    // for(int ch = 0; ch < (int) ref_signal.size(); ++ch)
    // {
    //     auto[h_ch, p_ch] = HPSS::hpss(ref_signal[ch], params);
    //     h_signal.push_back(h_ch);
    //     p_signal.push_back(p_ch);
    //
    //     std::vector<float> sum_ch (h_ch.size(), 0.0f);
    //     for(int i = 0; i < (int) h_ch.size(); ++i)
    //         sum_ch[i] = h_ch[i] + p_ch[i];
    //
    //     sum_signal.push_back(sum_ch);
    // }

    wav_io::write_file ("ref.wav", ref_signal, sf_info, arena);
    // WavIO::write_file("harmonic.wav", h_signal, sf_info);
    // WavIO::write_file("percussive.wav", p_signal, sf_info);
    // WavIO::write_file("sum.wav", sum_signal, sf_info);

    return 0;
}
