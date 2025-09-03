#include <iostream>
#include <cmath>
#include <algorithm>
#include <execution>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uSignal/filterDesign/response.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include "src/utilities/math/polynomial.hpp"
#include "src/alignment.hpp"


namespace
{
template<typename T>
USignal::Vector<std::complex<T>> makeComplex(const USignal::Vector<T> &x)
{
    USignal::Vector<std::complex<T>> cx(x.size());
    const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
    auto cxPtr = std::assume_aligned<ALIGNMENT> (cx.data());
    auto n = static_cast<int> (x.size());
    for (int i = 0; i < n; ++i)
    {
        cxPtr[i] = std::complex<double> (xPtr[i], 0);
    }
    return cx;
}

template<typename T>
USignal::Vector<std::complex<T>> makeComplexAndReverse(
    const USignal::Vector<T> &x)
{
    USignal::Vector<std::complex<T>> cx(x.size());
    const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
    auto cxPtr = std::assume_aligned<ALIGNMENT> (cx.data());
    auto n = static_cast<int> (x.size());
    for (int i = 0; i < n; ++i)
    {
        cxPtr[i] = std::complex<double> (xPtr[n - 1 - i], 0);
    }
    return cx;
}

template<typename T>
void computeTransferFunction(
    const USignal::Vector<std::complex<T>> &numeratorCoefficients,
    const USignal::Vector<std::complex<T>> &denominatorCoefficients,
    const USignal::Vector<std::complex<T>> &abscissa,
    USignal::Vector<std::complex<T>> *response)
{
#ifndef NDEBUG
    assert(hNumerator.size() == hDenominator.size());
    assert(response != nullptr);
#endif
    // Evaluate the numerator and denominator polynomials
    auto hNumerator 
        = USignal::Utilities::Math::Polynomial::evaluate(
            numeratorCoefficients, abscissa); 
    auto hDenominator
        = USignal::Utilities::Math::Polynomial::evaluate(
            denominatorCoefficients, abscissa);
    // Compute the transfer function
    if (response->size() != hNumerator.size())
    {
        response->resize(hNumerator.size());
    }
    response->resize(hNumerator.size());
#ifndef NDEBUG
    assert(hsNumerator.size() == response.size());
    assert(hsDenominator.size() == response.size());
#endif
    const auto hnPtr = std::assume_aligned<ALIGNMENT> (hNumerator.data());
    const auto hdPtr = std::assume_aligned<ALIGNMENT> (hDenominator.data());
    auto rPtr = std::assume_aligned<ALIGNMENT> (response->data());
    std::transform(std::execution::unseq,
                   hnPtr, hnPtr + hNumerator.size(), 
                   hdPtr,
                   rPtr,
                   std::divides< std::complex<T> > ());
}

}


template<typename T>
USignal::Vector<std::complex<T>>
USignal::FilterDesign::Response::computeAnalog(
   const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filter,
   const USignal::Vector<T> &frequencies)
{
    USignal::Vector<std::complex<T>> response;
    if (frequencies.empty()){return response;}
    // Get numerator and denominator coefficients
    const auto &numeratorCoefficients
        = filter.getNumeratorFilterCoefficientsReference();
    const auto &denominatorCoefficients
        = filter.getDenominatorFilterCoefficientsReference();
    auto allZero = std::all_of(denominatorCoefficients.cbegin(),
                               denominatorCoefficients.cend(),
                               [](const auto &a)
                               {
                                   constexpr T zero{0};
                                   return a == zero;
                               });
    if (allZero)
    {
        throw std::invalid_argument("All denominator coefficients are zero");
    }
    // Compute s = i \omega
    auto nFrequencies = static_cast<int> (frequencies.size());
    USignal::Vector<std::complex<T>> s(nFrequencies);
    auto sPtr = std::assume_aligned<ALIGNMENT> (s.data());
    const auto frequenciesPtr
        = std::assume_aligned<ALIGNMENT> (frequencies.data());
    for (int i = 0; i < nFrequencies; ++i)
    {
        sPtr[i] = std::complex<T> (0, frequenciesPtr[i]);
    }
    // Evaluate the numerator and denominator polynomials
    auto complexNumeratorCoefficients
        = ::makeComplex(numeratorCoefficients);
    auto complexDenominatorCoefficients
        = ::makeComplex(denominatorCoefficients);
    ::computeTransferFunction(complexNumeratorCoefficients,
                              complexDenominatorCoefficients,
                              s,
                              &response);
/*
    auto hsNumerator 
        = USignal::Utilities::Math::Polynomial::evaluate(
            ::makeComplex(numeratorCoefficients), s);
    auto hsDenominator
        = USignal::Utilities::Math::Polynomial::evaluate(
            ::makeComplex(denominatorCoefficients), s); 
    // Compute the transfer function
    response.resize(nFrequencies);
#ifndef NDEBUG
    assert(hsNumerator.size() == response.size());
    assert(hsDenominator.size() == response.size());
#endif
    std::transform(//std::execution::unseq,
                   hsNumerator.cbegin(),
                   hsDenominator.cend(),
                   hsDenominator.cbegin(),
                   response.begin(),
                   std::divides< std::complex<T> > ());
std::cout << "return " << std::endl;
*/
    return response;
}

template<typename T>
USignal::Vector<std::complex<T>>
USignal::FilterDesign::Response::computeDigital(
   const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filter,
   const USignal::Vector<T> &frequencies)
{
    USignal::Vector<std::complex<T>> response;
    if (frequencies.empty()){return response;}
    // Get numerator and denominator coefficients
    const auto &numeratorCoefficients
        = filter.getNumeratorFilterCoefficientsReference();
    const auto &denominatorCoefficients
        = filter.getDenominatorFilterCoefficientsReference();
    auto allZero = std::all_of(denominatorCoefficients.cbegin(),
                               denominatorCoefficients.cend(),
                               [](const auto &a) 
                               {
                                   constexpr T zero{0};
                                   return a == zero;
                               }); 
    if (allZero)
    {   
        throw std::invalid_argument("All denominator coefficients are zero");
    }
    // Copy numerators and denominators in reverse order for consistency with
    // evaluatePolynomial and then compute z.  The evaluatePolynomial convention
    // requires us to compute z=e^{-i\omega} to compensate.  
    // Normally, this is would be z^{i \omega}.
    auto reversedComplexNumeratorCoefficients
        = ::makeComplexAndReverse(numeratorCoefficients);
    auto reversedComplexDenominatorCoefficients
        = ::makeComplexAndReverse(denominatorCoefficients);
    // Now compute z = exp(-i w) = cos(w) - i sin(w)
    auto nFrequencies = static_cast<int> (frequencies.size());
    USignal::Vector<std::complex<T>> z(nFrequencies);
    auto zPtr = std::assume_aligned<ALIGNMENT> (z.data());
    const auto frequenciesPtr
        = std::assume_aligned<ALIGNMENT> (frequencies.data());
    for (int i = 0; i < nFrequencies; ++i)
    {
        const T w{frequenciesPtr[i]};
        zPtr[i] = std::complex<T> (std::cos(w), -std::sin(w));
    }
    // Now compute the transfer function
    ::computeTransferFunction(reversedComplexNumeratorCoefficients,
                              reversedComplexDenominatorCoefficients,
                              z,
                              &response);
    return response;
}

///--------------------------------------------------------------------------///
///                             Instantiation                                ///
///--------------------------------------------------------------------------///

template USignal::Vector<std::complex<double>>
USignal::FilterDesign::Response::computeAnalog(
    const USignal::FilterRepresentations::InfiniteImpulseResponse<double> &,
    const USignal::Vector<double> &);
template USignal::Vector<std::complex<float>>
USignal::FilterDesign::Response::computeAnalog(
    const USignal::FilterRepresentations::InfiniteImpulseResponse<float> &,
    const USignal::Vector<float> &); 


template USignal::Vector<std::complex<double>>
USignal::FilterDesign::Response::computeDigital(
    const USignal::FilterRepresentations::InfiniteImpulseResponse<double> &,
    const USignal::Vector<double> &); 
template USignal::Vector<std::complex<float>>
USignal::FilterDesign::Response::computeDigital(
    const USignal::FilterRepresentations::InfiniteImpulseResponse<float> &,
    const USignal::Vector<float> &); 

