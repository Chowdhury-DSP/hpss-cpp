#pragma once

#if __GNUC__ || __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic ignored "-Wshadow-field-in-constructor"
#pragma GCC diagnostic ignored "-Wdeprecated-dynamic-exception-spec"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wvla-extension"
#endif

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4005 4100)
#endif

#if defined(__ARM_NEON__)
#define PFFFT_ENABLE_NEON
#endif

#include "../pffft/pffft.h"
#include "../pffft/pffft.c" // NOLINT
#include "../pffft/pffft_common.c" // NOLINT

#if __GNUC__ || __clang__
#pragma GCC diagnostic pop
#endif

#if defined(_MSC_VER)
#pragma warning(pop)
#endif
