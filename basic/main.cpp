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



enum class AARoundOffMode {
    Everytime = 0,
    Dedicated,
    IncludesNonAffine
};

template <std::floating_point FloatType, AARoundOffMode roundOffMode = AARoundOffMode::Dedicated>
class AAFloat {
    using ia_fp_t = IAFloat<FloatType>;

    static uint32_t s_ID;

    static uint32_t getNewNoiseSymbol() {
        return s_ID++;
    }

    struct Empty {};

    std::unordered_map<uint32_t, FloatType> m_coeffs;
    [[no_unique_address]] std::conditional_t<
        roundOffMode != AARoundOffMode::Everytime,
        FloatType,
        Empty> m_roundOffCoeff;

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
        if constexpr (roundOffMode != AARoundOffMode::Everytime)
            m_roundOffCoeff = 0;
    }
    AAFloat(const ia_fp_t &x) {
        m_coeffs[0] = x.center();
        m_coeffs[getNewNoiseSymbol()] = x.radius();
        if constexpr (roundOffMode != AARoundOffMode::Everytime)
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
        if constexpr (roundOffMode != AARoundOffMode::Everytime) {
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
        if constexpr (roundOffMode == AARoundOffMode::Everytime) {
            m_coeffs[getNewNoiseSymbol()] = accRoundOff.hi();
        }
        else {
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
        if constexpr (roundOffMode == AARoundOffMode::Everytime) {
            m_coeffs[getNewNoiseSymbol()] = accRoundOff.hi();
        }
        else {
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
        if constexpr (roundOffMode == AARoundOffMode::Everytime) {
            m_coeffs[getNewNoiseSymbol()] = accRoundOff.hi();
        }
        else {
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
        if constexpr (roundOffMode == AARoundOffMode::Everytime) {
            m_coeffs[getNewNoiseSymbol()] = accRoundOff.hi();
        }
        else {
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
        if constexpr (roundOffMode != AARoundOffMode::Everytime) {
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
        if constexpr (roundOffMode != AARoundOffMode::Everytime) {
            v += r.m_roundOffCoeff;
            ia_fp_t ic = abs(c0) * r.m_roundOffCoeff;
            accRoundOff += ic;
        }

        // Non-affine term
        {
            // Quick conservative estimate
            ia_fp_t ic = u * v;
            ic += accRoundOff;
            if constexpr (roundOffMode != AARoundOffMode::IncludesNonAffine) {
                // New noise symbol for the non-affine term handles round-off errors also.
                m_coeffs[getNewNoiseSymbol()] = ic.hi();
                if constexpr (roundOffMode == AARoundOffMode::Dedicated)
                    m_roundOffCoeff = 0;
            }
            else {
                m_roundOffCoeff = ic.hi();
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
        if constexpr (roundOffMode == AARoundOffMode::Everytime) {
            m_coeffs[getNewNoiseSymbol()] = accRoundOff.hi();
        }
        else {
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

        // Affine terms
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
        if constexpr (roundOffMode != AARoundOffMode::Everytime)
            idelta += abs(ialpha) * v.m_roundOffCoeff;

        // Non-affine term
        {
            idelta += accRoundOff;
            if constexpr (roundOffMode != AARoundOffMode::IncludesNonAffine) {
                // New noise symbol for the non-affine term handles round-off errors also.
                ret.m_coeffs[getNewNoiseSymbol()] = idelta.hi();
            }
            else {
                ret.m_roundOffCoeff = idelta.hi();
            }
        }

        return ret;
    }
    AAFloat &operator/=(const AAFloat &r) {
        return *this *= reciprocal(r);
    }
    AAFloat &operator/=(FloatType r) {
        if (std::isinf(r)) {
            m_coeffs.clear();
            m_coeffs[0] = 0;
            if constexpr (roundOffMode != AARoundOffMode::Everytime)
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

        // Affine terms
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
        if constexpr (roundOffMode != AARoundOffMode::Everytime)
            idelta += abs(ialpha) * v.m_roundOffCoeff;

        // Non-affine term
        {
            idelta += accRoundOff;
            if constexpr (roundOffMode != AARoundOffMode::IncludesNonAffine) {
                // New noise symbol for the non-affine term handles round-off errors also.
                ret.m_coeffs[getNewNoiseSymbol()] = idelta.hi();
            }
            else {
                ret.m_roundOffCoeff = idelta.hi();
            }
        }

        return ret;
    }
};

template <std::floating_point FloatType, AARoundOffMode roundOffMode> uint32_t AAFloat<FloatType, roundOffMode>::s_ID = 1;



template <std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator+(const AAFloat<FloatType, roundOffMode> &a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret += b;
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode, Number N>
inline AAFloat<FloatType, roundOffMode> operator+(const AAFloat<FloatType, roundOffMode> &a, N b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret += static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator+(N a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator-(const AAFloat<FloatType, roundOffMode> &a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret -= b;
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode, Number N>
inline AAFloat<FloatType, roundOffMode> operator-(const AAFloat<FloatType, roundOffMode> &a, N b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret -= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator-(N a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = -b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator*(const AAFloat<FloatType, roundOffMode> &a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret *= b;
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode, Number N>
inline AAFloat<FloatType, roundOffMode> operator*(const AAFloat<FloatType, roundOffMode> &a, N b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret *= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator*(N a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = b;
    ret *= static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator/(const AAFloat<FloatType, roundOffMode> &a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret /= b;
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode, Number N>
inline AAFloat<FloatType, roundOffMode> operator/(const AAFloat<FloatType, roundOffMode> &a, N b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret /= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator/(N a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret(a);
    ret /= b;
    return ret;
}



using ia_fp32_t = IAFloat<float>;
using aa_fp32_t = AAFloat<float, AARoundOffMode::Dedicated>;
using ia_fp64_t = IAFloat<double>;
using aa_fp64_t = AAFloat<double, AARoundOffMode::Dedicated>;



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

    //{
    //    ia_fp32_t::setImaginaryValueHandling(true);

    //    constexpr uint32_t N = 16;
    //    for (int i = 0; i < N; ++i) {
    //        //ia_fp32_t iintv(-2 + 4.0f * i / N, -2 + 4.0f * (i + 1) / N);
    //        //ia_fp32_t iv1 = g(iintv);
    //        //ia_fp32_t iv2 = g(iv1);
    //        //devPrintf(
    //        //    "%g, %g, %g, %g\n",
    //        //    iv1.lo(), iv1.hi(), iv2.lo(), iv2.hi());

    //        aa_fp32_t aaintv(-2 + 4.0f * i / N, -2 + 4.0f * (i + 1) / N);
    //        aa_fp32_t aav1 = g(aaintv);
    //        aa_fp32_t aav2 = g(aav1);
    //        auto iaav1 = static_cast<ia_fp32_t>(aav1);
    //        auto iaav2 = static_cast<ia_fp32_t>(aav2);
    //        devPrintf(
    //            "%g, %g, %g, %g\n",
    //            iaav1.lo(), iaav1.hi(), iaav2.lo(), iaav2.hi());
    //    }
    //}

    {
        using FloatType = double;
        using ia_fp_t = IAFloat<FloatType>;
        using aa1_fp_t = AAFloat<FloatType, AARoundOffMode::Everytime>;
        using aa2_fp_t = AAFloat<FloatType, AARoundOffMode::Dedicated>;
        using aa3_fp_t = AAFloat<FloatType, AARoundOffMode::IncludesNonAffine>;
        constexpr size_t s1 = sizeof(aa1_fp_t);
        constexpr size_t s2 = sizeof(aa2_fp_t);
        constexpr size_t s3 = sizeof(aa3_fp_t);

        ia_fp_t ia_x(-1e-5, 1e-5);
        ia_fp_t ia_y(-1e-5, 1e-5);
        aa1_fp_t aa1_x(ia_x);
        aa1_fp_t aa1_y(ia_y);
        aa2_fp_t aa2_x(ia_x);
        aa2_fp_t aa2_y(ia_y);
        aa3_fp_t aa3_x(ia_x);
        aa3_fp_t aa3_y(ia_y);
        {
            auto ia_aa1_x = static_cast<ia_fp_t>(aa1_x);
            auto ia_aa1_y = static_cast<ia_fp_t>(aa1_y);
            auto ia_aa2_x = static_cast<ia_fp_t>(aa2_x);
            auto ia_aa2_y = static_cast<ia_fp_t>(aa2_y);
            auto ia_aa3_x = static_cast<ia_fp_t>(aa3_x);
            auto ia_aa3_y = static_cast<ia_fp_t>(aa3_y);
            devPrintf(
                "%g, %g, %g, %g\n",
                std::max(ia_x.radius(), ia_y.radius()),
                std::max(ia_aa1_x.radius(), ia_aa1_y.radius()),
                std::max(ia_aa2_x.radius(), ia_aa2_y.radius()),
                std::max(ia_aa3_x.radius(), ia_aa3_y.radius()));
        }
        constexpr FloatType a = 1.05;
        constexpr FloatType b = 0.3;
        for (int i = 0; i < 200; ++i) {
            ia_fp_t temp_ia_x = ia_x;
            ia_x = 1 - a * (temp_ia_x * temp_ia_x) + ia_y;
            ia_y = b * temp_ia_x;

            ia_fp_t ia_aa1_x;
            ia_fp_t ia_aa1_y;
            {
                aa1_fp_t temp_aa1_x = aa1_x;
                aa1_x = 1 - a * temp_aa1_x * temp_aa1_x + aa1_y;
                aa1_y = b * temp_aa1_x;
                ia_aa1_x = static_cast<ia_fp_t>(aa1_x);
                ia_aa1_y = static_cast<ia_fp_t>(aa1_y);
            }

            ia_fp_t ia_aa2_x;
            ia_fp_t ia_aa2_y;
            {
                aa2_fp_t temp_aa2_x = aa2_x;
                aa2_x = 1 - a * temp_aa2_x * temp_aa2_x + aa2_y;
                aa2_y = b * temp_aa2_x;
                ia_aa2_x = static_cast<ia_fp_t>(aa2_x);
                ia_aa2_y = static_cast<ia_fp_t>(aa2_y);
            }

            ia_fp_t ia_aa3_x;
            ia_fp_t ia_aa3_y;
            {
                aa3_fp_t temp_aa3_x = aa3_x;
                aa3_x = 1 - a * temp_aa3_x * temp_aa3_x + aa3_y;
                aa3_y = b * temp_aa3_x;
                ia_aa3_x = static_cast<ia_fp_t>(aa3_x);
                ia_aa3_y = static_cast<ia_fp_t>(aa3_y);
            }

            devPrintf(
                "%g, %g, %g, %g\n",
                std::max(ia_x.radius(), ia_y.radius()),
                std::max(ia_aa1_x.radius(), ia_aa1_y.radius()),
                std::max(ia_aa2_x.radius(), ia_aa2_y.radius()),
                std::max(ia_aa3_x.radius(), ia_aa3_y.radius()));
        }
        devPrintf("");
    }

    return 0;
}
