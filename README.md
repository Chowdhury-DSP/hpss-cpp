# Harmonic-Percussive Source Separation

This repository contains a C++ library for computing
[harmonic-percussive source separation](https://www.audiolabs-erlangen.de/resources/MIR/FMP/C8/C8S1_HPS.html).
The library intends to provide fast and high-quality source
separation, while also remaining suitable for real-time use-cases.

## Compilation

### Compiling with CMake

The library exports a single static library (`hpss`) as a CMake target.
```cpp
add_subdirectory(hpss-cpp)
target_link_libraries(MyCoolProject PRIVATE hpss)
```

## Usage

### Using the Simple API

After initializing a processor, you may process your signal
by repeatedly calling `hpss::process_window()`.

```cpp
#include <hpss/hpss.cpp>

// Set up parameters
const hpss::Params params {
    .window_size = 1 << 12,
    .hop_factor = 2,
    .zero_pad = 2,
    .mask_power = 2,
};

// Initialize an HPSS processor
auto processor = hpss::init (params);

// Process one "hop" of audio data
const auto [harmonic_out, percussive_out] = hpss::process_window (processor, window);

// De-initialize the processor
hpss::deinit (processor);
```

### Using the Complete API

Most of the "intermediate" steps of the algorithm are exposed as
their own standalone methods. These can be called independently,
which is useful in some cases.

For example, in real-time processing, the audio device buffer size
may be smaller than the hop size used  by the HPSS algorithm. When
using the Simple API, the user would  have wait until enough buffers
have accumulated to process a single "hop", at which point the entire
algorithm would run. In this scenario, the entire algorithm would
need to be fast enough to meet the real-time deadline determined by
the device buffer size.

With the Complete API, the work being done by the algorithm can
be split up into multiple steps. This allows an implementation
where a little bit of work is done with each incoming buffer,
while waiting for the next hop to accumulate, thereby reducing
the worst-case computation time for any one buffer.

Please review the implementation of `hpss::process_window()` to
see how the individual methods of the Complete API should be used.

## Testing

1. Configure and compile the test program
```bash
$ cmake -Bbuild -G<generator> -DHPSS_BUILD_TESTS=ON
$ cmake --build build --target hpss_test --parallel
```

2. Run the test program
```bash
$ .build/.../hpss_test audio_file.wav <num_seconds> <start_trim_seconds>
```

The test program will generate 4 files:
- ref.wav: A trimmed copy of the input file.
- harmonic.wav: Just the harmonic part of the signal.
- percussive.wav: Just the percussive part of the signal.
- sum.wav: The sum of the harmonic and percussive signals.

## References

The main algorithm being implemented is based on Derry Fitzgerald's
paper from DAFx10 ["Harmonic/Percussive Separation Using Median Filtering"](https://arrow.tudublin.ie/cgi/viewcontent.cgi?article=1078&context=argcon).

The median computation is implemented using a min-max heap
algorithm adapted from https://ideone.com/8VVEa (MIT license).

The FFT is computed using [`chowdsp_fft`](https://github.com/Chowdhury-DSP/chowdsp_fft).

## Future Improvements

- Currently the slowest part of the algorithm is the median
  computation. There's probably some ways that we can improve on it.
- It would be nice to either support more FFT backends (FFTW,
  Accelerate, etc), or allow the user to provide their own FFT backend.
- I've experimented with some SIMD optimizations, but have abandoned
  them for now, as I was not seeing a measurable performance improvement.
  However, it would be nice to take another shot at some point!

## LICENSE

`hpss-cpp` is published under the BSD 3-clause license.
