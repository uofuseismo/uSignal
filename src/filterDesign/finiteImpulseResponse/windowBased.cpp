#include <iostream>
#include <cmath>
#include <string>
#include "uSignal/filterDesign/finiteImpulseResponse/windowBased.hpp"
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#ifndef WITH_IPP
static_assert(false, "FIR window design requires IPP");
#else
#include <ipp/ipps.h>
#endif

using namespace USignal::FilterDesign::FiniteImpulseResponse::WindowBased;

namespace
{

/*
void sinc(const int n, const T x[], T sincx[])
{
    double pix;
    for (int i = 0; i < n; ++i)
    {
        sincx[i] = 1;
        if (x[i] != 0)
        {
            pix = M_PI*x[i];
            sincx[i] = std::sin(pix)/pix;
        }
    }
    return;
}
*/

IppWinType classifyWindow(const Window window)
{
    IppWinType winType = ippWinHamming; // Matlab default
    if (window == Window::Hamming)
    {
        winType = ippWinHamming;
    }
    else if (window == Window::Bartlett)
    {
        winType = ippWinBartlett;
    }
    else if (window == Window::Hanning)
    {
        winType = ippWinHann;
    }
    else if (window == Window::OptimalBlackman)
    {
        winType = ippWinBlackman;
    }
    return winType;
}

}

/// Lowpass
template<typename T>
USignal::FilterRepresentations::FiniteImpulseResponse<T>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::lowPass(
    const int order, const T r, const Window window)
{
    // Check inputs
    if (order < 4 || r <= 0 || r >= 1)
    {
        if (order < 4)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                     + " must be at least 4");
        }
        throw std::invalid_argument("normalized frequency = "
                                  + std::to_string(r)
                                  + " must be in range (0,1)");
    }
    auto windowType = ::classifyWindow(window);
    const int nTaps{order + 1};    // Length of filter is order + 1
    int bufferSize;
    auto status = ippsFIRGenGetBufferSize(nTaps, &bufferSize);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("Error getting buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufferSize);
    if (pBuffer == nullptr)
    {
        throw std::runtime_error("Failed to allocate buffer");
    }
    USignal::Vector<T> pTaps(nTaps);
    constexpr IppBool doNormal{ippTrue}; // Normalize filter coefficients
    const double rFreq{r/2};             // IPP uses 0.5 as Nyquist
    if constexpr (sizeof(T) == sizeof(double))
    {
        status = ippsFIRGenLowpass_64f(rFreq, pTaps.data(), nTaps,
                                       windowType, doNormal, pBuffer);
        ippsFree(pBuffer);
    }
    else if constexpr (sizeof(T) == sizeof(float))
    {
        USignal::Vector<double> pTaps64(nTaps);
        status = ippsFIRGenLowpass_64f(rFreq, pTaps64.data(), nTaps,
                                       windowType, doNormal, pBuffer);
        ippsFree(pBuffer);
        std::copy(pTaps64.begin(), pTaps64.end(), pTaps.begin());
    }
    else
    {
        ippsFree(pBuffer);
        static_assert(false, "Unhandled lowpass precision");
    }
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("FIR lowpass window design failed");
    }
    USignal::FilterRepresentations::FiniteImpulseResponse<T> fir{pTaps};
    return fir;
}

/// Highpass
template<typename T>
USignal::FilterRepresentations::FiniteImpulseResponse<T>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::highPass(
    const int order, const T r, const Window window)
{
    // Check inputs
    if (order < 4 || r <= 0 || r >= 1)
    {   
        if (order < 4)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                     + " must be at least 4");
        }
        throw std::invalid_argument("normalized frequency = "
                                  + std::to_string(r)
                                  + " must be in range (0,1)");
    }
    auto windowType = ::classifyWindow(window);
    const int nTaps{order + 1};    // Length of filter is order + 1
    int bufferSize;
    auto status = ippsFIRGenGetBufferSize(nTaps, &bufferSize);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("Error getting buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufferSize);
    if (pBuffer == nullptr)
    {
        throw std::runtime_error("Failed to allocate buffer");
    }
    USignal::Vector<T> pTaps(nTaps);
    constexpr IppBool doNormal{ippTrue}; // Normalize filter coefficients
    const double rFreq{r/2};             // IPP uses 0.5 as Nyquist
    if constexpr (sizeof(T) == sizeof(double))
    {
        status = ippsFIRGenHighpass_64f(rFreq, pTaps.data(), nTaps,
                                        windowType, doNormal, pBuffer);
        ippsFree(pBuffer);
    }
    else if constexpr (sizeof(T) == sizeof(float))
    {
        USignal::Vector<double> pTaps64(nTaps);
        status = ippsFIRGenHighpass_64f(rFreq, pTaps64.data(), nTaps,
                                        windowType, doNormal, pBuffer);
        ippsFree(pBuffer);
        std::copy(pTaps64.begin(), pTaps64.end(), pTaps.begin());
    }
    else
    {
        ippsFree(pBuffer);
        static_assert(false, "Unhandled highpass precision");
    }
    if (status != ippStsNoErr) 
    {
        throw std::runtime_error("FIR highpass window design failed");
    }
    USignal::FilterRepresentations::FiniteImpulseResponse<T> fir{pTaps};
    return fir;
}

/// Bandpass
template<typename T>
USignal::FilterRepresentations::FiniteImpulseResponse<T>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::bandPass(
    const int order, const std::pair<T, T> &r, const Window window)
{
    // Check inputs
    if (order < 4 ||
        r.first <= 0 || r.first >= 1 ||
        r.second <= 0 || r.second >= 1)
    {   
        if (order < 4)
        {   
            throw std::invalid_argument("order = " + std::to_string(order)
                                     + " must be at least 4");
        }
        if (r.first <= 0 || r.first >= 1)
        {
            throw std::invalid_argument("normalized low pass frequency = "
                                      + std::to_string(r.first)
                                      + " must be in range (0,1)");
        }
        throw std::invalid_argument("normalized high stop frequency = "
                                  + std::to_string(r.second)
                                  + " must be in range (0,1)");
    }   
    auto windowType = ::classifyWindow(window);
    const int nTaps{order + 1};    // Length of filter is order + 1
    int bufferSize;
    auto status = ippsFIRGenGetBufferSize(nTaps, &bufferSize);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("Error getting buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufferSize);
    if (pBuffer == nullptr)
    {   
        throw std::runtime_error("Failed to allocate buffer");
    }
    USignal::Vector<T> pTaps(nTaps);
    constexpr IppBool doNormal{ippTrue}; // Normalize filter coefficients
    const double rFreq0{r.first/2};      // IPP uses 0.5 as Nyquist
    const double rFreq1{r.second/2};     // IPP uses 0.5 as Nyquist
    if constexpr (sizeof(T) == sizeof(double))
    {
        status = ippsFIRGenBandpass_64f(rFreq0, rFreq1, pTaps.data(), nTaps,
                                        windowType, doNormal, pBuffer);
        ippsFree(pBuffer);
    }
    else if constexpr (sizeof(T) == sizeof(float))
    {
        USignal::Vector<double> pTaps64(nTaps);
        status = ippsFIRGenBandpass_64f(rFreq0, rFreq1, pTaps64.data(), nTaps,
                                        windowType, doNormal, pBuffer);
        ippsFree(pBuffer);
        std::copy(pTaps64.begin(), pTaps64.end(), pTaps.begin());
    }
    else
    {
        ippsFree(pBuffer);
        static_assert(false, "Unhandled bandpass precision");
    }
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("FIR bandpass window design failed");
    }
    USignal::FilterRepresentations::FiniteImpulseResponse<T> fir{pTaps};
    return fir;
}

/// Bandstop
template<typename T>
USignal::FilterRepresentations::FiniteImpulseResponse<T>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::bandStop(
    const int order, const std::pair<T, T> &r, const Window window)
{
    // Check inputs
    if (order < 4 ||
        r.first <= 0 || r.first >= 1 ||
        r.second <= 0 || r.second >= 1)
    {
        if (order < 4)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                     + " must be at least 4");
        }
        if (r.first <= 0 || r.first >= 1)
        {
            throw std::invalid_argument("normalized low stop frequency = "
                                      + std::to_string(r.first)
                                      + " must be in range (0,1)");
        }
        throw std::invalid_argument("normalized high pass frequency = "
                                  + std::to_string(r.second)
                                  + " must be in range (0,1)");
    }
    auto windowType = ::classifyWindow(window);
    const int nTaps{order + 1};    // Length of filter is order + 1
    int bufferSize;
    auto status = ippsFIRGenGetBufferSize(nTaps, &bufferSize);
    if (status != ippStsNoErr)
    {   
        throw std::runtime_error("Error getting buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufferSize);
    if (pBuffer == nullptr)
    {   
        throw std::runtime_error("Failed to allocate buffer");
    }
    USignal::Vector<T> pTaps(nTaps);
    constexpr IppBool doNormal{ippTrue}; // Normalize filter coefficients
    const double rFreq0{r.first/2};      // IPP uses 0.5 as Nyquist
    const double rFreq1{r.second/2};     // IPP uses 0.5 as Nyquist
    if constexpr (sizeof(T) == sizeof(double))
    {
        status = ippsFIRGenBandstop_64f(rFreq0, rFreq1, pTaps.data(), nTaps,
                                        windowType, doNormal, pBuffer);
        ippsFree(pBuffer);
    }
    else if constexpr (sizeof(T) == sizeof(float))
    {
        USignal::Vector<double> pTaps64(nTaps);
        status = ippsFIRGenBandstop_64f(rFreq0, rFreq1, pTaps64.data(), nTaps,
                                        windowType, doNormal, pBuffer);
        ippsFree(pBuffer);
        std::copy(pTaps64.begin(), pTaps64.end(), pTaps.begin());
    }
    else
    {
        ippsFree(pBuffer);
        static_assert(false, "Unhandled bandstop precision");
    }
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("FIR bandstop window design failed");
    }
    USignal::FilterRepresentations::FiniteImpulseResponse<T> fir{pTaps};
    return fir;
}

///--------------------------------------------------------------------------///
///                              Instantiation                               ///
///--------------------------------------------------------------------------///
template
USignal::FilterRepresentations::FiniteImpulseResponse<double>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::lowPass(
    int n,
    double,
    USignal::FilterDesign::FiniteImpulseResponse::WindowBased::Window );
template
USignal::FilterRepresentations::FiniteImpulseResponse<float>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::lowPass(
    int n,
    float,
    USignal::FilterDesign::FiniteImpulseResponse::WindowBased::Window );

template
USignal::FilterRepresentations::FiniteImpulseResponse<double>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::highPass(
    int n,
    double,
    USignal::FilterDesign::FiniteImpulseResponse::WindowBased::Window );
template
USignal::FilterRepresentations::FiniteImpulseResponse<float>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::highPass(
    int n,
    float,
    USignal::FilterDesign::FiniteImpulseResponse::WindowBased::Window );

template
USignal::FilterRepresentations::FiniteImpulseResponse<double>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::bandPass(
    int n,
    const std::pair<double, double> &,
    USignal::FilterDesign::FiniteImpulseResponse::WindowBased::Window );
template
USignal::FilterRepresentations::FiniteImpulseResponse<float>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::bandPass(
    int n,
    const std::pair<float, float> &,
    USignal::FilterDesign::FiniteImpulseResponse::WindowBased::Window );

template
USignal::FilterRepresentations::FiniteImpulseResponse<double>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::bandStop(
    int n,
    const std::pair<double, double> &,
    USignal::FilterDesign::FiniteImpulseResponse::WindowBased::Window );
template
USignal::FilterRepresentations::FiniteImpulseResponse<float>
USignal::FilterDesign::FiniteImpulseResponse::WindowBased::bandStop(
    int n,
    const std::pair<float, float> &,
    USignal::FilterDesign::FiniteImpulseResponse::WindowBased::Window );

