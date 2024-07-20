#pragma once

#include <type_traits>

namespace hpss::power
{
template <int exp, typename T>
constexpr std::enable_if_t<exp == 0, T>
    ipow (T) noexcept
{
    return (T) 1;
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 1, T>
    ipow (T a) noexcept
{
    return a;
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 2, T>
    ipow (T a) noexcept
{
    return a * a;
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 3, T>
    ipow (T a) noexcept
{
    return ipow<2> (a) * a;
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 4, T>
    ipow (T a) noexcept
{
    const auto b = ipow<2> (a);
    return ipow<2> (b);
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 5, T>
    ipow (T a) noexcept
{
    return ipow<4> (a) * a;
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 6, T>
    ipow (T a) noexcept
{
    const auto b = ipow<2> (a);
    return ipow<3> (b);
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 7, T>
    ipow (T a) noexcept
{
    return ipow<6> (a) * a;
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 8, T>
    ipow (T a) noexcept
{
    const auto d = ipow<4> (a);
    return ipow<2> (d);
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 9, T>
    ipow (T a) noexcept
{
    const auto c = ipow<3> (a);
    return ipow<3> (c);
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 10, T>
    ipow (T a) noexcept
{
    const auto e = ipow<5> (a);
    return ipow<2> (e);
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 11, T>
    ipow (T a) noexcept
{
    return ipow<10> (a) * a;
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 12, T>
    ipow (T a) noexcept
{
    const auto d = ipow<4> (a);
    return ipow<3> (d);
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 13, T>
    ipow (T a) noexcept
{
    return ipow<12> (a) * a;
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 14, T>
    ipow (T a) noexcept
{
    const auto g = ipow<7> (a);
    return ipow<2> (g);
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 15, T>
    ipow (T a) noexcept
{
    const auto e = ipow<5> (a);
    return ipow<3> (e);
}

template <int exp, typename T>
constexpr std::enable_if_t<exp == 16, T>
    ipow (T a) noexcept
{
    const auto h = ipow<8> (a);
    return ipow<2> (h);
}
} // namespace hpss::power
