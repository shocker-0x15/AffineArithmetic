/*

References:
- Affine Arithmetic and its Applications to Computer Graphics
- Affine Arithmeticについて
- Fast reliable interrogation of procedurally defined implicit surfaces
  using extended revised affine arithmetic

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
class IAFloat {
    static bool s_truncateImaginaryValue;

    FloatType m_lo;
    FloatType m_hi;

public:
    static void setImaginaryValueHandling(bool truncate) {
        s_truncateImaginaryValue = truncate;
    }
    static bool truncateImaginaryValue() {
        return s_truncateImaginaryValue;
    }

    IAFloat(FloatType x = 0) : m_lo(x), m_hi(x) {}
    IAFloat(FloatType lo, FloatType hi) : m_lo(lo), m_hi(hi) {}

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

    IAFloat operator+() const {
        return *this;
    }
    IAFloat operator-() const {
        IAFloat ret;
        ret.m_lo = -m_hi;
        ret.m_hi = -m_lo;
        return ret;
    }

    IAFloat &operator+=(const IAFloat &r) {
        std::fesetround(FE_DOWNWARD);
        m_lo += r.m_lo;
        std::fesetround(FE_UPWARD);
        m_hi += r.m_hi;
        std::fesetround(FE_TONEAREST);

        return *this;
    }
    IAFloat &operator+=(FloatType r) {
        return *this += IAFloat(r);
    }
    IAFloat &operator-=(const IAFloat &r) {
        std::fesetround(FE_DOWNWARD);
        m_lo -= r.m_hi;
        std::fesetround(FE_UPWARD);
        m_hi -= r.m_lo;
        std::fesetround(FE_TONEAREST);

        return *this;
    }
    IAFloat &operator-=(FloatType r) {
        return *this -= IAFloat(r);
    }
    IAFloat &operator*=(const IAFloat &r) {
        IAFloat l = *this;
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
    IAFloat &operator*=(FloatType r) {
        IAFloat l = *this;
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
    IAFloat &operator/=(const IAFloat &r) {
        IAFloat l = *this;
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
            throw std::domain_error("IAFloat: division by 0.");
        }
        std::fesetround(FE_TONEAREST);

        return *this;
    }
    IAFloat &operator/=(FloatType r) {
        IAFloat l = *this;
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
            throw std::domain_error("IAFloat: division by 0.");
        }
        std::fesetround(FE_TONEAREST);

        return *this;
    }

    friend inline IAFloat abs(const IAFloat &v) {
        IAFloat ret;
        FloatType absLo = std::fabs(v.m_lo);
        FloatType absHi = std::fabs(v.m_hi);
        if (absLo > absHi)
            std::swap(absLo, absHi);
        ret.m_lo = absLo;
        ret.m_hi = absHi;

        return ret;
    }
    friend inline IAFloat pow2(const IAFloat &v) {
        IAFloat ret;
        FloatType absLo = std::fabs(v.m_lo);
        FloatType absHi = std::fabs(v.m_hi);
        if (absLo > absHi)
            std::swap(absLo, absHi);
        std::fesetround(FE_DOWNWARD);
        ret.m_lo = absLo * absLo;
        std::fesetround(FE_UPWARD);
        ret.m_hi = absHi * absHi;
        std::fesetround(FE_TONEAREST);

        return ret;
    }
    friend inline IAFloat sqrt(const IAFloat &v) {
        IAFloat ret;
        IAFloat mv = v;
        if (mv.m_lo < 0.0f) {
            if (IAFloat::truncateImaginaryValue())
                mv.m_lo = 0.0f;
            else
                throw std::domain_error("IAFloat: sqrt of a negative value.");
        }

        std::fesetround(FE_DOWNWARD);
        ret.m_lo = std::sqrt(mv.m_lo);
        std::fesetround(FE_UPWARD);
        ret.m_hi = std::sqrt(mv.m_hi);
        std::fesetround(FE_TONEAREST);

        return ret;
    }
};

template <std::floating_point FloatType>
bool IAFloat<FloatType>::s_truncateImaginaryValue = false;



template <std::floating_point FloatType>
inline IAFloat<FloatType> operator+(const IAFloat<FloatType> &a, const IAFloat<FloatType> &b) {
    IAFloat<FloatType> ret = a;
    ret += b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline IAFloat<FloatType> operator+(const IAFloat<FloatType> &a, N b) {
    IAFloat<FloatType> ret = a;
    ret += static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType>
inline IAFloat<FloatType> operator+(N a, const IAFloat<FloatType> &b) {
    IAFloat<FloatType> ret = b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline IAFloat<FloatType> operator-(const IAFloat<FloatType> &a, const IAFloat<FloatType> &b) {
    IAFloat<FloatType> ret = a;
    ret -= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline IAFloat<FloatType> operator-(const IAFloat<FloatType> &a, N b) {
    IAFloat<FloatType> ret = a;
    ret -= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType>
inline IAFloat<FloatType> operator-(N a, const IAFloat<FloatType> &b) {
    IAFloat<FloatType> ret = -b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline IAFloat<FloatType> operator*(const IAFloat<FloatType> &a, const IAFloat<FloatType> &b) {
    IAFloat<FloatType> ret = a;
    ret *= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline IAFloat<FloatType> operator*(const IAFloat<FloatType> &a, N b) {
    IAFloat<FloatType> ret = a;
    ret *= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType>
inline IAFloat<FloatType> operator*(N a, const IAFloat<FloatType> &b) {
    IAFloat<FloatType> ret = b;
    ret *= static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline IAFloat<FloatType> operator/(const IAFloat<FloatType> &a, const IAFloat<FloatType> &b) {
    IAFloat<FloatType> ret = a;
    ret /= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline IAFloat<FloatType> operator/(const IAFloat<FloatType> &a, N b) {
    IAFloat<FloatType> ret = a;
    ret /= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType>
inline IAFloat<FloatType> operator/(N a, const IAFloat<FloatType> &b) {
    IAFloat<FloatType> ret(static_cast<FloatType>(a));
    ret /= b;
    return ret;
}



template <std::floating_point FloatType>
class AAFloat {
    using ia_fp_t = IAFloat<FloatType>;

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
    AAFloat(FloatType x = 0) {
        m_coeffs[0] = x;
        m_roundOffCoeff = 0;
    }
    AAFloat(const ia_fp_t &x) {
        m_coeffs[0] = x.center();
        m_coeffs[getNewNoiseSymbol()] = x.radius();
        m_roundOffCoeff = 0;
    }
    AAFloat(FloatType lo, FloatType hi) : AAFloat(ia_fp_t(lo, hi)) {}

    operator IAFloat<FloatType>() const {
        ia_fp_t ret(m_coeffs.at(0));
        for (auto it : m_coeffs) {
            if (it.first == 0)
                continue;
            FloatType d = std::fabs(it.second);
            ret += ia_fp_t(-d, d);
        }
        {
            ret += ia_fp_t(-m_roundOffCoeff, m_roundOffCoeff);
        }
        return ret;
    }

    AAFloat operator+() const {
        return *this;
    }
    AAFloat operator-() const {
        AAFloat ret = *this;
        for (auto &it : ret.m_coeffs)
            it.second *= -1;
        return ret;
    }

    AAFloat &operator+=(const AAFloat &r) {
        ia_fp_t accRoundOff;
        for (auto &it : m_coeffs) {
            if (!r.m_coeffs.contains(it.first))
                continue;
            ia_fp_t ic(it.second);
            ic += r.m_coeffs.at(it.first);
            it.second = ic.center();
            accRoundOff += ic.radius();
        }
        for (auto &it : r.m_coeffs) {
            if (m_coeffs.contains(it.first))
                continue;
            m_coeffs[it.first] = it.second;
        }
        {
            accRoundOff += m_roundOffCoeff;
            accRoundOff += r.m_roundOffCoeff;
            m_roundOffCoeff = accRoundOff.hi();
        }

        cleanUp();

        return *this;
    }
    AAFloat &operator+=(FloatType r) {
        ia_fp_t accRoundOff;
        {
            FloatType &c0 = m_coeffs.at(0);
            ia_fp_t ic(c0);
            ic += r;
            c0 = ic.center();
            accRoundOff += ic.radius();
        }
        {
            accRoundOff += m_roundOffCoeff;
            m_roundOffCoeff = accRoundOff.hi();
        }
        return *this;
    }
    AAFloat &operator-=(const AAFloat &r) {
        ia_fp_t accRoundOff;
        for (auto &it : m_coeffs) {
            if (!r.m_coeffs.contains(it.first))
                continue;
            ia_fp_t ic(it.second);
            ic -= r.m_coeffs.at(it.first);
            it.second = ic.center();
            accRoundOff += ic.radius();
        }
        for (auto &it : r.m_coeffs) {
            if (m_coeffs.contains(it.first))
                continue;
            m_coeffs[it.first] = -it.second;
        }
        {
            // The round-off coefficients don't cancel each other.
            accRoundOff += m_roundOffCoeff;
            accRoundOff += r.m_roundOffCoeff;
            m_roundOffCoeff = accRoundOff.hi();
        }

        cleanUp();

        return *this;
    }
    AAFloat &operator-=(FloatType r) {
        ia_fp_t accRoundOff;
        {
            FloatType &c0 = m_coeffs.at(0);
            ia_fp_t ic(c0);
            ic -= r;
            c0 = ic.center();
            accRoundOff += ic.radius();
        }
        {
            accRoundOff += m_roundOffCoeff;
            m_roundOffCoeff = accRoundOff.hi();
        }
        return *this;
    }
    AAFloat &operator*=(const AAFloat &r) {
        ia_fp_t accRoundOff;
        ia_fp_t u = 0;
        ia_fp_t v = 0;
        ia_fp_t c0 = m_coeffs.at(0);
        ia_fp_t r_c0 = r.m_coeffs.at(0);

        // Affine terms
        {
            ia_fp_t ic = c0 * r_c0;
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
            ia_fp_t ic = c0 * r_ci + r_c0 * it.second;
            it.second = ic.center();
            accRoundOff += ic.radius();
        }
        {
            u += m_roundOffCoeff;
            ia_fp_t ic = abs(r_c0) * m_roundOffCoeff;
            accRoundOff += ic;
        }
        for (auto &it : r.m_coeffs) {
            if (it.first == 0)
                continue;

            v += std::fabs(it.second);

            if (m_coeffs.contains(it.first))
                continue;
            ia_fp_t ic = c0 * it.second;
            m_coeffs[it.first] = ic.center();
            accRoundOff += ic.radius();
        }
        {
            v += r.m_roundOffCoeff;
            ia_fp_t ic = abs(c0) * r.m_roundOffCoeff;
            accRoundOff += ic;
        }

        // Non-Affine term
        {
            // Quick conservative estimate
            ia_fp_t ic = u * v;
            if constexpr (false) {
                m_coeffs[getNewNoiseSymbol()] = ic.center();
                accRoundOff += ic.radius();
                m_roundOffCoeff = accRoundOff.hi();
            }
            else {
                ic += accRoundOff;
                // New noise symbol for the non-affine term handles round-off errors also.
                m_coeffs[getNewNoiseSymbol()] = ic.hi();
                m_roundOffCoeff = 0;
            }
        }

        cleanUp();

        return *this;
    }
    AAFloat &operator*=(FloatType r) {
        ia_fp_t accRoundOff;
        for (auto &it : m_coeffs) {
            ia_fp_t ic(it.second);
            ic *= r;
            it.second = ic.center();
            accRoundOff += ic.radius();
        }
        {
            ia_fp_t ic(m_roundOffCoeff);
            ic *= std::fabs(r);
            accRoundOff += ic;
            m_roundOffCoeff = accRoundOff.hi();
        }

        cleanUp();

        return *this;
    }
    friend inline AAFloat reciprocal(const AAFloat &v) {
        ia_fp_t interval = static_cast<ia_fp_t>(v);
        if (interval.lo() <= 0.0f && interval.hi() >= 0.0f)
            throw std::domain_error("AAFloat: division by 0.");

        ia_fp_t ia(interval.lo());
        ia_fp_t ib(interval.hi());
        ia_fp_t iab = ia * ib;
        ia_fp_t ialpha = -1 / iab;
        ia_fp_t isqrtab = (interval.lo() > 0.0f ? 1 : -1) * sqrt(iab);

        ia_fp_t ibeta = (ia + 2 * isqrtab + ib) / (2 * iab);
        ia_fp_t idelta = (ia - 2 * isqrtab + ib) / (2 * iab);

        ia_fp_t accRoundOff;

        ia_fp_t ic0(ialpha * v.m_coeffs.at(0) + ibeta);
        AAFloat ret(ic0.center());
        accRoundOff += ic0.radius();

        for (auto it : v.m_coeffs) {
            if (it.first == 0)
                continue;
            ia_fp_t ic(ialpha * it.second);
            ret.m_coeffs[it.first] = ic.center();
            accRoundOff += ic.radius();
        }
        accRoundOff += ia_fp_t(abs(ialpha) * v.m_roundOffCoeff);
        accRoundOff += idelta;

        // New noise symbol for the non-affine term handles round-off errors also.
        ret.m_coeffs[getNewNoiseSymbol()] = accRoundOff.hi();
        ret.m_roundOffCoeff = 0;

        return ret;
    }
    AAFloat &operator/=(const AAFloat &r) {
        return *this *= reciprocal(r);
    }
    AAFloat &operator/=(FloatType r) {
        if (std::isinf(r)) {
            m_coeffs.clear();
            m_coeffs[0] = 0;
            m_roundOffCoeff = 0;
        }
        else {
            *this *= AAFloat(1.0f / ia_fp_t(r));
        }
        return *this;
    }

    friend inline AAFloat pow2(const AAFloat &v) {
        return v * v;
    }
    friend inline AAFloat sqrt(const AAFloat &v) {
        ia_fp_t interval = static_cast<ia_fp_t>(v);
        if (interval.lo() < 0.0f)
            throw std::domain_error("IAFloat: sqrt of a negative value.");

        ia_fp_t isqrt = sqrt(interval);
        ia_fp_t isqrt_a(isqrt.lo());
        ia_fp_t isqrt_b(isqrt.hi());
        ia_fp_t ialpha = 1 / (isqrt_a + isqrt_b);
        ia_fp_t ibeta = (isqrt_a + isqrt_b) / 8 + (isqrt_a * isqrt_b) / (isqrt_a + isqrt_b) / 2;
        ia_fp_t idelta = pow2(isqrt_b - isqrt_a) / (8 * (isqrt_a + isqrt_b));

        ia_fp_t accRoundOff;

        ia_fp_t ic0(ialpha * v.m_coeffs.at(0) + ibeta);
        AAFloat ret(ic0.center());
        accRoundOff += ic0.radius();

        for (auto it : v.m_coeffs) {
            if (it.first == 0)
                continue;
            ia_fp_t ic(ialpha * it.second);
            ret.m_coeffs[it.first] = ic.center();
            accRoundOff += ic.radius();
        }
        {
            ia_fp_t ic(ialpha * v.m_roundOffCoeff);
            accRoundOff += ic.radius();
        }

        accRoundOff += idelta;
        //ret.m_roundOffCoeff = accRoundOff.hi();
        // New noise symbol for the non-affine term handles round-off errors also.
        ret.m_coeffs[getNewNoiseSymbol()] = accRoundOff.hi();
        ret.m_roundOffCoeff = 0;

        return ret;
    }
};

template <std::floating_point FloatType> uint32_t AAFloat<FloatType>::s_ID = 1;



template <std::floating_point FloatType>
inline AAFloat<FloatType> operator+(const AAFloat<FloatType> &a, const AAFloat<FloatType> &b) {
    AAFloat<FloatType> ret = a;
    ret += b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline AAFloat<FloatType> operator+(const AAFloat<FloatType> &a, N b) {
    AAFloat<FloatType> ret = a;
    ret += static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType>
inline AAFloat<FloatType> operator+(N a, const AAFloat<FloatType> &b) {
    AAFloat<FloatType> ret = b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline AAFloat<FloatType> operator-(const AAFloat<FloatType> &a, const AAFloat<FloatType> &b) {
    AAFloat<FloatType> ret = a;
    ret -= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline AAFloat<FloatType> operator-(const AAFloat<FloatType> &a, N b) {
    AAFloat<FloatType> ret = a;
    ret -= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType>
inline AAFloat<FloatType> operator-(N a, const AAFloat<FloatType> &b) {
    AAFloat<FloatType> ret = -b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline AAFloat<FloatType> operator*(const AAFloat<FloatType> &a, const AAFloat<FloatType> &b) {
    AAFloat<FloatType> ret = a;
    ret *= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline AAFloat<FloatType> operator*(const AAFloat<FloatType> &a, N b) {
    AAFloat<FloatType> ret = a;
    ret *= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType>
inline AAFloat<FloatType> operator*(N a, const AAFloat<FloatType> &b) {
    AAFloat<FloatType> ret = b;
    ret *= static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType>
inline AAFloat<FloatType> operator/(const AAFloat<FloatType> &a, const AAFloat<FloatType> &b) {
    AAFloat<FloatType> ret = a;
    ret /= b;
    return ret;
}

template <std::floating_point FloatType, Number N>
inline AAFloat<FloatType> operator/(const AAFloat<FloatType> &a, N b) {
    AAFloat<FloatType> ret = a;
    ret /= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType>
inline AAFloat<FloatType> operator/(N a, const AAFloat<FloatType> &b) {
    AAFloat<FloatType> ret(a);
    ret /= b;
    return ret;
}



using ia_fp32_t = IAFloat<float>;
using aa_fp32_t = AAFloat<float>;



template <typename ValueType>
inline ValueType g(const ValueType &x) {
    ValueType a = pow2(x);
    return sqrt(a - x + 0.5f) / sqrt(a + 0.5f);
}

int32_t main(int32_t argc, const char* argv[]) {
    //{
    //    ia_fp32_t x(-2.0f, 2.0f);
    //    ia_fp32_t r(-1.0f, 1.0f);
    //    ia_fp32_t s(-1.0f, 1.0f);
    //    ia_fp32_t z = (10 + x + r) * (10 - x + s);
    //    //auto sqx = sqrt(x);
    //    printf("");
    //}
    //{
    //    FP32Affine x(-2.0f, 2.0f);
    //    FP32Affine r(-1.0f, 1.0f);
    //    FP32Affine s(-1.0f, 1.0f);
    //    FP32Affine z = (10 + x + r) * (10 - x + s);
    //    auto iz = static_cast<ia_fp32_t>(z);
    //    auto sqz = sqrt(z);
    //    printf("");
    //}

    {
        ia_fp32_t::setImaginaryValueHandling(true);

        constexpr uint32_t N = 16;
        for (int i = 0; i < N; ++i) {
            //ia_fp32_t iintv(-2 + 4.0f * i / N, -2 + 4.0f * (i + 1) / N);
            //ia_fp32_t iv1 = g(iintv);
            //ia_fp32_t iv2 = g(iv1);
            //devPrintf(
            //    "%g, %g, %g, %g\n",
            //    iv1.lo(), iv1.hi(), iv2.lo(), iv2.hi());

            aa_fp32_t aaintv(-2 + 4.0f * i / N, -2 + 4.0f * (i + 1) / N);
            aa_fp32_t aav1 = g(aaintv);
            aa_fp32_t aav2 = g(aav1);
            auto iaav1 = static_cast<ia_fp32_t>(aav1);
            auto iaav2 = static_cast<ia_fp32_t>(aav2);
            devPrintf(
                "%g, %g, %g, %g\n",
                iaav1.lo(), iaav1.hi(), iaav2.lo(), iaav2.hi());
        }
    }

    return 0;
}
