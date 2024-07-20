#pragma once

#include <cstddef>

namespace hpss::simd
{
namespace sse_or_neon
{
    void multiply_4 (const float* a, const float* b, float* y, size_t count);

    void multiply_add_4 (const float* a, const float* b, float* y, size_t count);
}

namespace avx
{
    void multiply_8 (const float* a, const float* b, float* y, size_t count);

    void multiply_add_8 (const float* a, const float* b, float* y, size_t count);
}
} // namespace hpss::simd
