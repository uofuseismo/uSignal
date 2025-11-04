#include <iostream>
#include <cmath>
#include <memory>
#include <complex>
#include <algorithm>
#include <limits>
#include "polynomial.hpp"
#include "src/alignment.hpp"

namespace
{
#ifdef WITH_MKL
#include <mkl_lapacke.h>
#ifndef lapack_int
#define lapack_int MKL_INT
#endif
#ifndef lapack_complex_float
#define lapack_complex_float   MKL_Complex8
#endif
#ifndef lapack_complex_double
#define lapack_complex_double   MKL_Complex16
#endif
template<typename T>
USignal::Vector<std::complex<T>> computeEigenvalues(
    USignal::Vector<double> &A, const int nIn, const int ldaIn)
{
    auto n = static_cast<lapack_int> (nIn);
    auto lda = static_cast<lapack_int> (ldaIn);
    USignal::Vector<std::complex<T>> eigenvalues(n, 0 + 0i);
    std::vector<double> wr(n);
    std::vector<double> wi(n);
    double vl;
    constexpr lapack_int ldvl{1};
    double vr;
    constexpr lapack_int ldvr{1};
    auto info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, A.data(), lda,
                              wr.data(), wi.data(), &vl, ldvl, &vr, ldvr);
    if (info > 0)
    {
        throw std::runtime_error("Failed to compute eigenvalues");
    }
    for (int i = 0; i < n; ++i)
    {
        auto wri = static_cast<T> (wr[i]);
        auto wii = static_cast<T> (wi[i]);
        eigenvalues[i] = std::complex<T> (wri, wii);
    }
    return eigenvalues;
}
#endif
}

//using namespace USignal::Utilities::Math::Polynomial;

/// @brief Evalutes the polynomial:
///           p(x) = p_0 x^n + p_1 x^{n-1} + \cdots + p_n  
///        at the given evaluation points.
/// @param[in] polynomialCoefficients  The polynomial coefficients
///                                    ordered as shown above.
/// @param[in] evaluationPoints The x_i's at which to evaluate the polynomial.
/// @result The polynomial evaluated at each evaluation point.
template<typename T>
USignal::Vector<T> 
USignal::Utilities::Math::Polynomial::evaluate(
    const USignal::Vector<T> &polynomialCoefficients,
    const USignal::Vector<T> &evaluationPoints) 
{
    if (polynomialCoefficients.empty())
    {
        throw std::invalid_argument(
           "No coefficients in polynomialCoeffiecients");
    }
    auto order = static_cast<int> (polynomialCoefficients.size()) - 1;
    auto nEvaluationPoints = static_cast<int> (evaluationPoints.size());
    USignal::Vector<T> y;
    if (nEvaluationPoints < 1){return y;}
    y.resize(nEvaluationPoints, 0);
    // Expand the constant case 
    auto yPtr = std::assume_aligned<ALIGNMENT> (y.data());
    if (order == 0) // Constant
    {
        const auto p0 = polynomialCoefficients[0];
        std::fill(yPtr, yPtr + nEvaluationPoints, p0);
    }
    else if (order == 1) // Linear
    {
        const auto xPtr
            = std::assume_aligned<ALIGNMENT> (evaluationPoints.data());
        const auto p0 = polynomialCoefficients[0];
        const auto p1 = polynomialCoefficients[1];
        for (int i = 0; i < nEvaluationPoints; ++i)
        {
            yPtr[i] = p0*xPtr[i] + p1;
        }
    }
    else if (order == 2) // Quadratic
    {
        const auto xPtr
            = std::assume_aligned<ALIGNMENT> (evaluationPoints.data());
        const auto p0 = polynomialCoefficients[0];
        const auto p1 = polynomialCoefficients[1];
        const auto p2 = polynomialCoefficients[2];
        for (int i = 0; i < nEvaluationPoints; ++i)
        {
            auto xi = xPtr[i];
            yPtr[i] = p2 + xi*(p1 + xi*p0);
        }
    }
    else // General case (Horner's rule)
    {
        const auto xPtr
            = std::assume_aligned<ALIGNMENT> (evaluationPoints.data());
        const auto p0 = polynomialCoefficients.at(0);
        for (auto i = 0; i < nEvaluationPoints; ++i)
        {
            yPtr[i] = p0*xPtr[i];
        }
        for (auto j = 1; j < order; j++)
        {
            const auto pj = polynomialCoefficients[j];
            for (auto i = 0; i < nEvaluationPoints; ++i)
            {
                yPtr[i] = (pj + yPtr[i])*xPtr[i];
            }
        }
        const auto pn = polynomialCoefficients[order];
        for (auto i = 0; i < nEvaluationPoints; ++i)
        {
            yPtr[i] = pn + yPtr[i];
        }
    }
    return y;
}

/// @brief Computes the roots of a polynomial:
///        p(x) = p_0 x^n + p_1 x^{n-1} + \cdots + p_n
/// @param[in] coefficients  The coefficients of the polynomial whose
///                          order is shown above.
/// @result The roots of the polynomial.
template<typename T>
USignal::Vector<std::complex<T>>
USignal::Utilities::Math::Polynomial::computeRoots(
    const USignal::Vector<T> &coefficients)
{
    auto nCoefficients = static_cast<int> (coefficients.size());
    if (nCoefficients < 1)
    {
        throw std::invalid_argument("No coefficeints");
    }
    auto order = nCoefficients - 1;
    constexpr T zero{0};
    if (coefficients[0] == zero)
    {
        throw std::invalid_argument("Highest order coefficient is zero");
    }
    // Set space for companion matrix
    int n   = order;
    int lda = std::max(8, order);
    USignal::Vector<double> A(lda*lda, 0.0);
    double ami = 1.0/static_cast<double > (coefficients[0]); //coefficient on highest order term
    // Fill out the non-zeros of the companion matrix
    for (int i = 1; i < order + 1; ++i)
    {
        int indx = (i - 1)*lda + 0;
        A[indx] =-ami*static_cast<double> (coefficients[i]);
        // One in subdiagonal 
        if (i < order)
        {
            indx = (i - 1)*lda + i;
            A[indx] = 1.0;
        }
    }
    return ::computeEigenvalues<T>(A, n, lda);
}

/// Expand a polynomial
template<typename T>
USignal::Vector<std::complex<T>>
USignal::Utilities::Math::Polynomial::expand(
    const USignal::Vector<std::complex<T>> &roots)
{
    constexpr std::complex<T> one{1 + 0i};
    USignal::Vector<std::complex<T>> result;
    auto order = static_cast<int> (roots.size());
    result.resize(order + 1);
    // Special case is order 0
    if (order == 0)
    {
        result[0] = one;
        return result;
    }
    // Linear case: 1*x - r_0
    if (order == 1)
    {
        result[0] = one;
        result[1] =-roots[0];
        return result;
    }
    constexpr std::complex<T> zero{0 + 0i};
    USignal::Vector<std::complex<T>> temp1(roots.size() + 1);
    USignal::Vector<std::complex<T>> temp2(roots.size() + 1);
    // Initialize
    result[0] =-roots[0];
    result[1] = one;
    for (int i = 2; i <= order; ++i)
    {
        // x*(a_0 + a_1 x + ... + a_n x^{n-1}) = a_0 x + a_1 x^2 + ... + a_n x^i 
        // Shift right
        temp1[0] = zero;
        for (int j = 1; j <= i; ++j)
        {
            temp1[j] = result[j - 1];
        }
        // -p_i*(a_0 + .... + a_n x^{i-1}) =-p_i a_0 - ... p_i a_n x^{i-1}
        const auto p_i = roots[i - 1];
        for (int j = 0; j < i; ++j)
        {
            temp2[j] = -p_i*result[j];
        }
        temp2[i] = zero;
        // -p_i a_0 + p_i a_1 x + ... + p_i a_n x^{i-1}
        //          -     a_0 x - ...                   - a_n x^i
        for (int j = 0; j < i + 1; ++j)
        {
            result[j] = temp1[j] + temp2[j];
        }
    }
    /*
    // Clean up some numerical junk
    constexpr T tolerance{std::numeric_limits<T>::epsilon()};
    for (auto &ci : result)
    {
        if (std::abs(ci) < tolerance){ci = std::real(ci) + 0i;}
    }
    */
    // Reverse for consistency with Matlab
    std::reverse(result.data(), result.data() + result.size());
    return result;
}

template<typename T>
USignal::Vector<T>
USignal::Utilities::Math::Polynomial::expandToRealCoefficients(
    const USignal::Vector<std::complex<T>> &roots)
{
    auto complexPolynomialCoefficients = expand(roots);
    USignal::Vector<T> polynomialCoefficients(
        complexPolynomialCoefficients.size());
    const auto zc
        = std::assume_aligned<ALIGNMENT>
          (complexPolynomialCoefficients.data());
    auto c = std::assume_aligned<ALIGNMENT> (polynomialCoefficients.data());
    for (int i = 0; i < static_cast<int> (polynomialCoefficients.size()); ++i)
    {
        c[i] = std::real(zc[i]);
    }
    return polynomialCoefficients;
}

///--------------------------------------------------------------------------///
///                             Instantiation                                ///
///--------------------------------------------------------------------------///
template USignal::Vector<double>
USignal::Utilities::Math::Polynomial::evaluate(
    const USignal::Vector<double> &polynomialCoefficients,
    const USignal::Vector<double> &evaluationPoints);
template USignal::Vector<float>
USignal::Utilities::Math::Polynomial::evaluate(
    const USignal::Vector<float> &polynomialCoefficients,
    const USignal::Vector<float> &evaluationPoints);
template USignal::Vector<std::complex<double>>
USignal::Utilities::Math::Polynomial::evaluate(
    const USignal::Vector<std::complex<double>> &polynomialCoefficients,
    const USignal::Vector<std::complex<double>> &evaluationPoints);
template USignal::Vector<std::complex<float>>
USignal::Utilities::Math::Polynomial::evaluate(
    const USignal::Vector<std::complex<float>> &polynomialCoefficients,
    const USignal::Vector<std::complex<float>> &evaluationPoints);

template USignal::Vector<std::complex<double>>
   USignal::Utilities::Math::Polynomial::computeRoots(const USignal::Vector<double> &coefficients);
template USignal::Vector<std::complex<float>>
   USignal::Utilities::Math::Polynomial::computeRoots(const USignal::Vector<float> &coefficients);

template USignal::Vector<double>
   USignal::Utilities::Math::Polynomial::expandToRealCoefficients(const USignal::Vector<std::complex<double>> &);
template USignal::Vector<float>
   USignal::Utilities::Math::Polynomial::expandToRealCoefficients(const USignal::Vector<std::complex<float>> &);

