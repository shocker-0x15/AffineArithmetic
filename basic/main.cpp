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
#include <string>
#include <vector>

#include "../affine_arithmetic.h"

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
    {
        ia_fp32_t x(-2.0f, 2.0f);
        ia_fp32_t r(-1.0f, 1.0f);
        ia_fp32_t s(-1.0f, 1.0f);
        ia_fp32_t z = (10 + x + r) * (10 - x + s);
        //auto sqx = sqrt(x);
        printf("");
    }
    {
        aa_fp32_t x(-2.0f, 2.0f);
        aa_fp32_t r(-1.0f, 1.0f);
        aa_fp32_t s(-1.0f, 1.0f);
        aa_fp32_t z = (10 + x + r) * (10 - x + s);
        auto iz = static_cast<ia_fp32_t>(z);
        auto sqz = sqrt(z);
        printf("");
    }

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
        devPrintf("\n");
    }

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
        devPrintf("\n");
    }

    return 0;
}
