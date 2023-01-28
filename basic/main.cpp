/*

References:
- Affine Arithmetic and its Applications to Computer Graphics
- Affine Arithmeticについて

*/

// Platform defines
#if defined(_WIN32) || defined(_WIN64)
#   define Platform_Windows
#   if defined(_MSC_VER)
#       define Platform_Windows_MSVC
#       if defined(__INTELLISENSE__)
#           define Platform_CodeCompletion
#       endif
#   endif
#elif defined(__APPLE__)
#   define Platform_macOS
#endif

#if defined(Platform_Windows_MSVC)
#   define WIN32_LEAN_AND_MEAN
#   define NOMINMAX
#   define _USE_MATH_DEFINES
#   include <Windows.h>
#   undef WIN32_LEAN_AND_MEAN
#   undef NOMINMAX
#   undef near
#   undef far
#   undef RGB
#endif

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cfloat>
#include <cfenv>
#include <cmath>
#include <stdexcept>
#include <numbers>
#include <numeric>
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>

#if defined(Platform_Windows_MSVC)
static void devPrintf(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char str[4096];
    vsnprintf_s(str, sizeof(str), _TRUNCATE, fmt, args);
    va_end(args);
    OutputDebugString(str);
}
#else
#   define devPrintf printf
#endif



template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;



template <std::floating_point FloatType>
class Interval {
    FloatType m_lo;
    FloatType m_hi;

public:
    Interval(FloatType x = 0) : m_lo(x), m_hi(x) {}
    Interval(FloatType lo, FloatType hi) : m_lo(lo), m_hi(hi) {}

    FloatType &lo() {
        return m_lo;
    }
    FloatType &hi() {
        return m_hi;
    }
    const FloatType &lo() const {
        return m_lo;
    }
    const FloatType &hi() const {
        return m_hi;
    }
    FloatType center() const {
        return (m_lo + m_hi) / 2;
    }
    FloatType radius() const {
        FloatType c = center();
        std::fesetround(FE_UPWARD);
        FloatType d = std::max(c - m_lo, m_hi - c);
        std::fesetround(FE_TONEAREST);
        return d;
    }

    Interval operator+() const {
        return *this;
    }
    Interval operator-() const {
        Interval ret;
        ret.m_lo = -m_hi;
        ret.m_hi = -m_lo;
        return ret;
    }

    Interval &operator+=(const Interval &r) {
        std::fesetround(FE_DOWNWARD);
        m_lo += r.m_lo;
        std::fesetround(FE_UPWARD);
        m_hi += r.m_hi;
        std::fesetround(FE_TONEAREST);

        return *this;
    }
    Interval &operator+=(FloatType r) {
        return *this += Interval(r);
    }
    Interval &operator-=(const Interval &r) {
        std::fesetround(FE_DOWNWARD);
        m_lo -= r.m_lo;
        std::fesetround(FE_UPWARD);
        m_hi -= r.m_hi;
        std::fesetround(FE_TONEAREST);

        return *this;
    }
    Interval &operator-=(FloatType r) {
        return *this -= Interval(r);
    }
    Interval &operator*=(const Interval &r) {
        Interval l = *this;
        if (l.m_lo >= 0.0f) {
            if (r.m_lo >= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_lo * r.m_lo;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_hi * r.m_hi;
            }
            else if (r.m_hi <= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_hi * r.m_lo;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_lo * r.m_hi;
            }
            else {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_hi * r.m_lo;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_hi * r.m_hi;
            }
        }
        else if (l.m_hi <= 0.0f) {
            if (r.m_lo >= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_lo * r.m_hi;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_hi * r.m_lo;
            }
            else if (r.m_hi <= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_hi * r.m_hi;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_lo * r.m_lo;
            }
            else {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_lo * r.m_hi;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_lo * r.m_lo;
            }
        }
        else {
            if (r.m_lo >= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_lo * r.m_hi;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_hi * r.m_hi;
            }
            else if (r.m_hi <= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_hi * r.m_lo;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_lo * r.m_lo;
            }
            else {
                std::fesetround(FE_DOWNWARD);
                m_lo = std::min(l.m_lo * r.m_hi, l.m_hi * r.m_lo);
                std::fesetround(FE_UPWARD);
                m_hi = std::max(l.m_lo * r.m_lo, l.m_hi * r.m_hi);
            }
        }
        std::fesetround(FE_TONEAREST);

        return *this;
    }
    Interval &operator*=(FloatType r) {
        Interval l = *this;
        if (r >= 0.0f) {
            std::fesetround(FE_DOWNWARD);
            m_lo = l.m_lo * r;
            std::fesetround(FE_UPWARD);
            m_hi = l.m_hi * r;
        }
        else {
            std::fesetround(FE_DOWNWARD);
            m_lo = l.m_hi * r;
            std::fesetround(FE_UPWARD);
            m_hi = l.m_lo * r;
        }
        std::fesetround(FE_TONEAREST);
        
        return *this;
    }
    Interval &operator/=(const Interval &r) {
        Interval l = *this;
        if (r.m_lo > 0.0f) {
            if (l.m_lo >= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_lo / r.m_hi;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_hi / r.m_lo;
            }
            else if (l.m_hi <= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_lo / r.m_lo;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_hi / r.m_hi;
            }
            else {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_lo / r.m_lo;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_hi / r.m_lo;
            }
        }
        else if (r.m_hi < 0.0f) {
            if (l.m_lo >= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_hi / r.m_hi;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_lo / r.m_lo;
            }
            else if (l.m_hi <= 0.0f) {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_hi / r.m_lo;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_lo / r.m_hi;
            }
            else {
                std::fesetround(FE_DOWNWARD);
                m_lo = l.m_hi / r.m_hi;
                std::fesetround(FE_UPWARD);
                m_hi = l.m_lo / r.m_hi;
            }
        }
        else {
            std::fesetround(FE_TONEAREST);
            throw std::domain_error("Interval: division by 0.");
        }
        std::fesetround(FE_TONEAREST);

        return *this;
    }
    Interval &operator/=(FloatType r) {
        Interval l = *this;
        if (r > 0.0f) {
            std::fesetround(FE_DOWNWARD);
            m_lo = l.m_lo / r;
            std::fesetround(FE_UPWARD);
            m_hi = l.m_hi / r;
        }
        else if (r < 0.0f) {
            std::fesetround(FE_DOWNWARD);
            m_lo = l.m_hi / r;
            std::fesetround(FE_UPWARD);
            m_hi = l.m_lo / r;
        }
        else {
            std::fesetround(FE_TONEAREST);
            throw std::domain_error("Interval: division by 0.");
        }
        std::fesetround(FE_TONEAREST);

        return *this;
    }

    template <std::floating_point FloatType>
    friend inline Interval<FloatType> sqrt(const Interval<FloatType> &v) {
        Interval<FloatType> ret;
        if (v.m_lo < 0.0f)
            throw std::domain_error("Interval: sqrt of a negative value.");

        std::fesetround(FE_DOWNWARD);
        ret.m_lo = std::sqrt(v.m_lo);
        std::fesetround(FE_UPWARD);
        ret.m_hi = std::sqrt(v.m_hi);
        std::fesetround(FE_TONEAREST);

        return ret;
    }
};



template <std::floating_point FloatType>
inline Interval<FloatType> operator+(const Interval<FloatType> &a, const Interval<FloatType> &b) {
    Interval<FloatType> ret = a;
    ret += b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Interval<FloatType> operator+(const Interval<FloatType> &a, N b) {
    Interval<FloatType> ret = a;
    ret += static_cast<FloatType>(b);
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Interval<FloatType> operator+(N a, const Interval<FloatType> &b) {
    Interval<FloatType> ret = b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline Interval<FloatType> operator-(const Interval<FloatType> &a, const Interval<FloatType> &b) {
    Interval<FloatType> ret = a;
    ret -= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Interval<FloatType> operator-(const Interval<FloatType> &a, N b) {
    Interval<FloatType> ret = a;
    ret -= static_cast<FloatType>(b);
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Interval<FloatType> operator-(N a, const Interval<FloatType> &b) {
    Interval<FloatType> ret = -b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline Interval<FloatType> operator*(const Interval<FloatType> &a, const Interval<FloatType> &b) {
    Interval<FloatType> ret = a;
    ret *= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Interval<FloatType> operator*(const Interval<FloatType> &a, N b) {
    Interval<FloatType> ret = a;
    ret *= static_cast<FloatType>(b);
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Interval<FloatType> operator*(N a, const Interval<FloatType> &b) {
    Interval<FloatType> ret = b;
    ret *= static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline Interval<FloatType> operator/(const Interval<FloatType> &a, const Interval<FloatType> &b) {
    Interval<FloatType> ret = a;
    ret /= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Interval<FloatType> operator/(const Interval<FloatType> &a, N b) {
    Interval<FloatType> ret = a;
    ret /= static_cast<FloatType>(b);
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Interval<FloatType> operator/(N a, const Interval<FloatType> &b) {
    Interval<FloatType> ret(static_cast<FloatType>(a));
    ret /= b;
    return ret;
}



template <std::floating_point FloatType>
class Affine {
    using IntervalF = Interval<FloatType>;

    static uint32_t s_ID;

    static uint32_t getNewNoiseSymbol() {
        return s_ID++;
    }

    std::unordered_map<uint32_t, FloatType> m_coeffs;
    FloatType m_roundOffCoeff;

    void cleanUp() {
        for (auto it = m_coeffs.begin(); it != m_coeffs.end();) {
            if (it->second == 0.0f && it->first > 0)
                it = m_coeffs.erase(it);
            else
                ++it;
        }
    }

public:
    Affine(FloatType x = 0) {
        m_coeffs[0] = x;
        m_roundOffCoeff = 0;
    }
    Affine(FloatType lo, FloatType hi) {
        IntervalF i(lo, hi);
        m_coeffs[0] = i.center();
        m_coeffs[getNewNoiseSymbol()] = i.radius();
        m_roundOffCoeff = 0;
    }

    operator Interval<FloatType>() const {
        FloatType c = m_coeffs.at(0);
        FloatType d = 0;
        for (auto it : m_coeffs) {
            if (it.first == 0)
                continue;
            d += std::fabs(it.second);
        }
        return Interval<FloatType>(c - d, c + d);
    }

    Affine operator+() const {
        return *this;
    }
    Affine operator-() const {
        Affine ret = *this;
        for (auto &it : ret.m_coeffs)
            it.second *= -1;
        return ret;
    }

    Affine &operator+=(const Affine &r) {
        IntervalF accRoundOff;
        for (auto &it : m_coeffs) {
            if (!r.m_coeffs.contains(it.first))
                continue;
            IntervalF ic(it.second);
            ic += r.m_coeffs.at(it.first);
            it.second = ic.center();
            accRoundOff += ic.radius();
        }
        for (auto &it : r.m_coeffs) {
            if (m_coeffs.contains(it.first))
                continue;
            m_coeffs[it.first] = it.second;
        }
        accRoundOff += m_roundOffCoeff;
        accRoundOff += r.m_roundOffCoeff;
        m_roundOffCoeff = accRoundOff.hi();

        cleanUp();

        return *this;
    }
    Affine &operator+=(FloatType r) {
        IntervalF accRoundOff;
        FloatType &c0 = m_coeffs.at(0);
        IntervalF ic(c0);
        ic += r;
        c0 = ic.center();
        accRoundOff += ic.radius();
        accRoundOff += m_roundOffCoeff;
        m_roundOffCoeff = accRoundOff.hi();
        return *this;
    }
    Affine &operator-=(const Affine &r) {
        IntervalF accRoundOff;
        for (auto &it : m_coeffs) {
            if (!r.m_coeffs.contains(it.first))
                continue;
            IntervalF ic(it.second);
            ic -= r.m_coeffs.at(it.first);
            it.second = ic.center();
            accRoundOff += ic.radius();
        }
        for (auto &it : r.m_coeffs) {
            if (m_coeffs.contains(it.first))
                continue;
            m_coeffs[it.first] = it.second;
        }
        accRoundOff += m_roundOffCoeff;
        accRoundOff += r.m_roundOffCoeff;
        m_roundOffCoeff = accRoundOff.hi();

        cleanUp();

        return *this;
    }
    Affine &operator-=(FloatType r) {
        IntervalF accRoundOff;
        FloatType &c0 = m_coeffs.at(0);
        IntervalF ic(c0);
        ic -= r;
        c0 = ic.center();
        accRoundOff += ic.radius();
        accRoundOff += m_roundOffCoeff;
        m_roundOffCoeff = accRoundOff.hi();
        return *this;
    }
    Affine &operator*=(const Affine &r) {
        IntervalF accRoundOff;
        IntervalF u = 0;
        IntervalF v = 0;
        IntervalF c0 = m_coeffs.at(0);
        IntervalF r_c0 = r.m_coeffs.at(0);

        // Affine terms
        {
            IntervalF ic = c0 * r_c0;
            m_coeffs.at(0) = ic.center();
            accRoundOff += ic.radius();
        }
        for (auto &it : m_coeffs) {
            if (it.first == 0)
                continue;

            u += std::fabs(it.second);

            FloatType r_ci = 0;
            if (r.m_coeffs.contains(it.first))
                r_ci = r.m_coeffs.at(it.first);
            IntervalF ic = c0 * r_ci + r_c0 * it.second;
            it.second = ic.center();
            accRoundOff += ic.radius();
        }
        {
            u += m_roundOffCoeff;
            accRoundOff += r_c0 * m_roundOffCoeff;
        }
        for (auto &it : r.m_coeffs) {
            if (it.first == 0)
                continue;

            v += std::fabs(it.second);

            if (m_coeffs.contains(it.first))
                continue;
            IntervalF ic = c0 * it.second;
            m_coeffs[it.first] = ic.center();
            accRoundOff += ic.radius();
        }
        {
            v += r.m_roundOffCoeff;
            accRoundOff += c0 * r.m_roundOffCoeff;
        }

        // Non-Affine term
        {
            // Quick conservative estimate
            accRoundOff += u * v;
            // New noise symbol for the non-affine term handles round-off errors also.
            m_coeffs[getNewNoiseSymbol()] = accRoundOff.hi();
        }

        m_roundOffCoeff = 0;

        cleanUp();

        return *this;
    }
    Affine &operator*=(FloatType r) {
        if (r == 0) {
            m_coeffs.clear();
            m_coeffs[0] = 0;
            m_roundOffCoeff = 0;
        }
        else {
            IntervalF accRoundOff;
            for (auto &it : m_coeffs) {
                IntervalF ic(it.second);
                ic *= r;
                it.second = ic.center();
                accRoundOff += ic.radius();
            }
            accRoundOff += IntervalF(m_roundOffCoeff) * std::fabs(r);
            m_roundOffCoeff = accRoundOff.hi();
        }
        return *this;
    }
};

template <std::floating_point FloatType> uint32_t Affine<FloatType>::s_ID = 1;



template <std::floating_point FloatType>
inline Affine<FloatType> operator+(const Affine<FloatType> &a, const Affine<FloatType> &b) {
    Affine<FloatType> ret = a;
    ret += b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Affine<FloatType> operator+(const Affine<FloatType> &a, N b) {
    Affine<FloatType> ret = a;
    ret += static_cast<FloatType>(b);
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Affine<FloatType> operator+(N a, const Affine<FloatType> &b) {
    Affine<FloatType> ret = b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline Affine<FloatType> operator-(const Affine<FloatType> &a, const Affine<FloatType> &b) {
    Affine<FloatType> ret = a;
    ret -= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Affine<FloatType> operator-(const Affine<FloatType> &a, N b) {
    Affine<FloatType> ret = a;
    ret -= static_cast<FloatType>(b);
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Affine<FloatType> operator-(N a, const Affine<FloatType> &b) {
    Affine<FloatType> ret = -b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline Affine<FloatType> operator*(const Affine<FloatType> &a, const Affine<FloatType> &b) {
    Affine<FloatType> ret = a;
    ret *= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Affine<FloatType> operator*(const Affine<FloatType> &a, N b) {
    Affine<FloatType> ret = a;
    ret *= static_cast<FloatType>(b);
    return ret;
}

template <std::floating_point FloatType, Number N>
inline Affine<FloatType> operator*(N a, const Affine<FloatType> &b) {
    Affine<FloatType> ret = b;
    ret *= static_cast<FloatType>(a);
    return ret;
}



using FP32Interval = Interval<float>;
using FP32Affine = Affine<float>;

int32_t main(int32_t argc, const char* argv[]) {
    {
        FP32Interval x(-2.0f, 2.0f);
        FP32Interval r(-1.0f, 1.0f);
        FP32Interval s(-1.0f, 1.0f);
        FP32Interval z = (10 + x + r) * (10 - x + s);
        printf("");
    }
    {
        FP32Affine x(-2.0f, 2.0f);
        FP32Affine r(-1.0f, 1.0f);
        FP32Affine s(-1.0f, 1.0f);
        FP32Affine z = (10 + x + r) * (10 - x + s);
        auto iz = static_cast<FP32Interval>(z);
        printf("");
    }

    return 0;
}
