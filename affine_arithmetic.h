#pragma once

#include "interval_arithmetic.h"
#include <unordered_map>



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

    static AAFloat affineApproximation(
        const AAFloat &v,
        const ia_fp_t &ialpha, const ia_fp_t &ibeta, ia_fp_t idelta) {
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
                if (idelta.hi() > 0)
                    ret.m_coeffs[getNewNoiseSymbol()] = idelta.hi();
            }
            else {
                ret.m_roundOffCoeff = idelta.hi();
            }
        }

        return ret;
    }

public:
    AAFloat(FloatType x = 0) {
        m_coeffs[0] = x;
        if constexpr (roundOffMode != AARoundOffMode::Everytime)
            m_roundOffCoeff = 0;
    }
    AAFloat(const ia_fp_t &x) {
        m_coeffs[0] = x.center();
        FloatType r = x.radius();
        if (r > 0)
            m_coeffs[getNewNoiseSymbol()] = r;
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
    IAFloat<FloatType> toIAFloat() const {
        return static_cast<ia_fp_t>(*this);
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
                if (ic.hi() > 0)
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
        ia_fp_t isqrtab = (interval.lo() > 0.0f ? 1 : -1) * sqrt(iab);

        ia_fp_t ialpha = -1 / iab;
        ia_fp_t ibeta = (ia + 2 * isqrtab + ib) / (2 * iab);
        ia_fp_t idelta = (ia - 2 * isqrtab + ib) / (2 * iab);

        return affineApproximation(v, ialpha, ibeta, idelta);
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

        return affineApproximation(v, ialpha, ibeta, idelta);
    }
};

template <std::floating_point FloatType, AARoundOffMode roundOffMode>
uint32_t AAFloat<FloatType, roundOffMode>::s_ID = 1;



template <std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator+(
    const AAFloat<FloatType, roundOffMode> &a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret += b;
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode, Number N>
inline AAFloat<FloatType, roundOffMode> operator+(
    const AAFloat<FloatType, roundOffMode> &a, N b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret += static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator+(
    N a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator-(
    const AAFloat<FloatType, roundOffMode> &a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret -= b;
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode, Number N>
inline AAFloat<FloatType, roundOffMode> operator-(
    const AAFloat<FloatType, roundOffMode> &a, N b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret -= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator-(
    N a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = -b;
    ret += static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator*(
    const AAFloat<FloatType, roundOffMode> &a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret *= b;
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode, Number N>
inline AAFloat<FloatType, roundOffMode> operator*(
    const AAFloat<FloatType, roundOffMode> &a, N b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret *= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator*(
    N a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = b;
    ret *= static_cast<FloatType>(a);
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator/(
    const AAFloat<FloatType, roundOffMode> &a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret /= b;
    return ret;
}

template <std::floating_point FloatType, AARoundOffMode roundOffMode, Number N>
inline AAFloat<FloatType, roundOffMode> operator/(
    const AAFloat<FloatType, roundOffMode> &a, N b) {
    AAFloat<FloatType, roundOffMode> ret = a;
    ret /= static_cast<FloatType>(b);
    return ret;
}

template <Number N, std::floating_point FloatType, AARoundOffMode roundOffMode>
inline AAFloat<FloatType, roundOffMode> operator/(
    N a, const AAFloat<FloatType, roundOffMode> &b) {
    AAFloat<FloatType, roundOffMode> ret(a);
    ret /= b;
    return ret;
}
