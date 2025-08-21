#ifndef USIGNAL_UTILITIES_MATH_POLYNOMIAL_HPP
#define USIGNAL_UTILITIES_MATH_POLYNOMIAL_HPP
#include <complex>
#include <uSignal/vector.hpp>
namespace USignal::Utilities::Math::Polynomial
{

/// @brief Evalutes the polynomial:
///           p(x) = p_0 x^n + p_1 x^{n-1} + \cdots + p_n  
///        at the given evaluation points.
/// @param[in] polynomialCoefficients  The polynomial coefficients
///                                    ordered as shown above.
/// @param[in] evaluationPoints The x_i's at which to evaluate the polynomial.
/// @result The polynomial evaluated at each evaluation point.
template<typename T>
USignal::Vector<T> 
evaluate(const USignal::Vector<T> &polynomialCoefficients,
            const USignal::Vector<T> &evaluationPoints);

/// @brief Computes the roots of a polynomial:
///        p(x) = p_0 x^n + p_1 x^{n-1} + \cdots + p_n
/// @param[in] coefficients  The coefficients of the polynomial whose
///                          order is shown above.
/// @result The roots of the polynomial.
template<typename T>
USignal::Vector<std::complex<T>>
   computeRoots(const USignal::Vector<T> &coefficients);

}
#endif
