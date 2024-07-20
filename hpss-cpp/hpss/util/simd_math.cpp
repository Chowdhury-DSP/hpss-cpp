#include "simd_math.hpp"

#include <cassert>

#if defined(__AVX2__)
#include <immintrin.h>
#elif defined(__SSE2__)
#include <immintrin.h>
#elif defined(__ARM_NEON__)
#include <arm_neon.h>
#define PFFFT_ENABLE_NEON
#else
static_assert (false, "Compiling for an un-supported architecture!");
#endif

namespace hpss::simd
{
#if COMPILING_WITH_AVX
namespace avx
{
    void multiply_4 (const float* a, const float* b, float* y, size_t count)
    {
        assert (count % 8 == 0);
        for (size_t i = 0; i < count; i += 8)
        {
            auto a_reg = _mm256_load_ps (a + i);
            auto b_reg = _mm256_load_ps (b + i);
            auto y_reg = _mm256_mul_ps (a_reg, b_reg);
            _mm256_store_ps (y + i, y_reg);
        }
    }

    void multiply_add_4 (const float* a, const float* b, float* y, size_t count)
    {
        assert (count % 8 == 0);
        for (size_t i = 0; i < count; i += 8)
        {
            auto a_reg = _mm256_load_ps (a + i);
            auto b_reg = _mm256_load_ps (b + i);
            auto y_reg = _mm256_load_ps (y + i);
            y_reg = _mm256_add_ps (_mm256_mul_ps (a_reg, b_reg), y_reg);
            _mm256_store_ps (y + i, y_reg);
        }
    }
}
#else
namespace sse_or_neon
{
    void multiply_4 (const float* a, const float* b, float* y, size_t count)
    {
        assert (count % 4 == 0);
        for (size_t i = 0; i < count; i += 4)
        {
            auto a_reg = _mm_load_ps (a + i);
            auto b_reg = _mm_load_ps (b + i);
            auto y_reg = _mm_mul_ps (a_reg, b_reg);
            _mm_store_ps (y + i, y_reg);
        }
    }

    void multiply_add_4 (const float* a, const float* b, float* y, size_t count)
    {
        assert (count % 4 == 0);
        for (size_t i = 0; i < count; i += 4)
        {
            auto a_reg = _mm_load_ps (a + i);
            auto b_reg = _mm_load_ps (b + i);
            auto y_reg = _mm_load_ps (y + i);
            y_reg = _mm_add_ps (_mm_mul_ps (a_reg, b_reg), y_reg);
            _mm_store_ps (y + i, y_reg);
        }
    }
} // namespace sse_or_neon
#endif
} // namespace hpss::simd
