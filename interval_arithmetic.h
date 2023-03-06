#pragma once

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
            if (truncateImaginaryValue())
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
