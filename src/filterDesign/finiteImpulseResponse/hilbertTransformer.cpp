#include <complex>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uSignal/filterDesign/finiteImpulseResponse/hilbertTransformer.hpp"
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include "src/alignment.hpp"
#include "src/filterDesign/windowFunctions.hpp"

namespace
{

template<typename T>
USignal::Vector<T> sinc(const int n)
{

    if (n < 1){throw std::invalid_argument("No samples in sinc");}
    // Compute the sinc function where fc = 1 and fc/2 = 0.5
    constexpr T one{1};
    constexpr T zero{0};
    USignal::Vector<T> sinc(n, one);
    auto sincPtr = std::assume_aligned<ALIGNMENT> (sinc.data());
    const double term1 = 0.25*static_cast<double> (1 - n);
    for (int i = 0; i < n; ++i)
    {
        auto t = term1 + 0.5*static_cast<double> (i);
        //sincPtr[i] = one;
        if (t != 0)
        {
            double pix = M_PI*t;
            sincPtr[i] = static_cast<T> (std::sin(pix)/pix);
        }
    }
    return sinc;
}

}

/// Hilbert transformer
template<typename T>
USignal::FilterRepresentations::FiniteImpulseResponse<std::complex<T>>
hilbertTransform( //USignal::FilterDesign::FiniteImpulseResponse::hilbertTransformer<T>(
    const int order, const double beta)
{
    if (order < 0)
    {
        throw std::invalid_argument("order cannot be negative");
    }
    int n = order + 1;
    // Special case
    if (n == 1)
    {
        std::vector<double> hfiltR(1, 1);
        std::vector<double> hfiltI(1, 0);
        USignal::Vector<std::complex<T>>
           taps( std::vector<std::complex<T>> {1 + 0i, 0 + 0i} );
        USignal::FilterRepresentations::FiniteImpulseResponse<std::complex<T>> fir{taps};
        return taps;
    }
    // Create a kaiser window
    USignal::Vector<double> kaiser(n, 0);
    auto kaiserPtr = std::assume_aligned<ALIGNMENT> (kaiser.data());
    ::kaiser(n, &kaiserPtr, beta);

    auto sinct = ::sinc<double> (n);
    // Create the complex valued FIR coefficients of 12.66 of Oppenheim and 
    // Schafer.  Note, that hfilt = sinc(t)*exp(i*pi*t) is what Matlab uses.
    // However Oppenheim and Schafer use sin*sinc in 12.67.  This uses the
    // Matlab implementation.
    const double term1 = 0.25*static_cast<double> (1 - n);
    USignal::Vector<std::complex<T>> taps(n, std::complex<T> (0, 0));
    auto sinctPtr = std::assume_aligned<ALIGNMENT> (sinct.data()); 
    auto tapsPtr = std::assume_aligned<ALIGNMENT> (taps.data());   
    double gain{0};
    if (n%2 == 1)
    {
        // Type III has many zeros
        gain = 1; 
        taps.at(n/2) = std::complex<T> (1, 0);
        // Two cases
        int istart = 0;
        if ((n/2)%2 == 0){istart = 1;}
        for (int i = istart; i < n; i = i + 2)
        {
            auto t = term1 + 0.5*static_cast<double> (i);
            double ks = kaiserPtr[i]*sinctPtr[i];
            auto hfiltI = static_cast<T> (ks*(std::sin(M_PI*t)));
            tapsPtr[i] = std::complex<T> (0, hfiltI);
        }
    }
    else
    {
        auto sinctPtr = std::assume_aligned<ALIGNMENT> (sinct.data());
        const double term1 = 0.25*static_cast<double> (1 - n);
        // Type IV is, in general, non-zero
        for (int i = 0; i < n; ++i)
        {
            auto t = term1 + 0.5*static_cast<double> (i);
            double ks = kaiserPtr[i]*sinctPtr[i];
            auto hfiltR = ks*std::cos(M_PI*t);
            gain = gain + hfiltR;
            auto hfiltI = ks*std::sin(M_PI*t);
            tapsPtr[i] = std::complex<T> (hfiltR, hfiltI);
        }
#ifndef NDEBUG
        assert(gain != 0);
#endif
        // Normalize
//        ippsDivC_64f_I(gain, hfiltR.data(), n); 
//        ippsDivC_64f_I(gain, hfiltI.data(), n);
    }
    // Set the filter
    USignal::FilterRepresentations::FiniteImpulseResponse<std::complex<T>>
        firFilter{taps};
    return firFilter;
}
