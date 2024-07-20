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
    void multiply_8 (const float* a, const float* b, float* y, size_t count)
    {
        assert (count % 8 == 0);
#if defined(__AVX2__)
        for (size_t i = 0; i < count; i += 8)
        {
            auto a_reg = _mm256_load_ps (a + i);
            auto b_reg = _mm256_load_ps (b + i);
            auto y_reg = _mm256_mul_ps (a_reg, b_reg);
            _mm256_store_ps (y + i, y_reg);
        }
#endif
    }

    void multiply_add_8 (const float* a, const float* b, float* y, size_t count)
    {
        assert (count % 8 == 0);
#if defined(__AVX2__)
        for (size_t i = 0; i < count; i += 8)
        {
            auto a_reg = _mm256_load_ps (a + i);
            auto b_reg = _mm256_load_ps (b + i);
            auto y_reg = _mm256_load_ps (y + i);
            y_reg = _mm256_add_ps (_mm256_mul_ps (a_reg, b_reg), y_reg);
            _mm256_store_ps (y + i, y_reg);
        }
#endif
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
#if defined(__SSE2__)
            auto a_reg = _mm_load_ps (a + i);
            auto b_reg = _mm_load_ps (b + i);
            auto y_reg = _mm_mul_ps (a_reg, b_reg);
            _mm_store_ps (y + i, y_reg);
#elif defined(__ARM_NEON__)
            auto a_reg = vld1q_f32(a + i);
            auto b_reg = vld1q_f32(b + i);
            auto y_reg = vmulq_f32(a_reg, b_reg);
            vst1q_f32(y + i, y_reg);
#endif
        }
    }

    void multiply_add_4 (const float* a, const float* b, float* y, size_t count)
    {
        assert (count % 4 == 0);
        for (size_t i = 0; i < count; i += 4)
        {
#if defined(__SSE2__)
            auto a_reg = _mm_load_ps (a + i);
            auto b_reg = _mm_load_ps (b + i);
            auto y_reg = _mm_load_ps (y + i);
            y_reg = _mm_add_ps (_mm_mul_ps (a_reg, b_reg), y_reg);
            _mm_store_ps (y + i, y_reg);
#elif defined(__ARM_NEON__)
            auto a_reg = vld1q_f32(a + i);
            auto b_reg = vld1q_f32(b + i);
            auto y_reg = vld1q_f32(y + i);
            y_reg = vmlaq_f32(y_reg, a_reg, b_reg);
            vst1q_f32(y + i, y_reg);
#endif
        }
    }
} // namespace sse_or_neon
#endif
} // namespace hpss::simd
