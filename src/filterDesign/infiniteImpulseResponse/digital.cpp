#include <iostream>
#include <cmath>
#include <complex>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uSignal/filterDesign/infiniteImpulseResponse/digital.hpp"
#include "uSignal/filterDesign/infiniteImpulseResponse/bilinearTransform.hpp"
#include "uSignal/filterDesign/infiniteImpulseResponse/analogPrototype.hpp"
#include "uSignal/filterDesign/infiniteImpulseResponse/convertBand.hpp"
#include "uSignal/vector.hpp"

using namespace USignal::FilterDesign::InfiniteImpulseResponse;
namespace UFR = USignal::FilterRepresentations;

namespace
{

enum class Prototype
{
    Bessel,
    Butterworth,
    ChebyshevI,
    ChebyshevII
};

[[nodiscard]] double warpDigital(const double normalizedFrequency)
{
    constexpr double fs{2};
    return 2*fs*std::tan(M_PI*normalizedFrequency/fs);
}
[[nodiscard]] std::pair<double, double> warpDigital(
     const std::pair<double, double> &normalizedCutoffFrequencies)
{
    return std::pair {::warpDigital(normalizedCutoffFrequencies.first),
                      ::warpDigital(normalizedCutoffFrequencies.second)}; 
}
void checkFrequencies(
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    auto w0 = normalizedCutoffFrequencies.first;
    if (w0 < 0)
    {
        throw std::invalid_argument("normalizedCutoffFrequencies.first = "
                                  +  std::to_string(w0)
                                  + " must be non-negative");
    }
    auto w1 = normalizedCutoffFrequencies.second;
    if (w1 > 1)
    {   
        throw std::invalid_argument("normalizedCutoffFrequencies.second = "
                                  +  std::to_string(w1)
                                  + " cannot exceed 1");
    }   
    if (w0 >= w1)
    {
        throw std::invalid_argument(
            "normalizedCutoffFrequencies.first = "
          + std::to_string(w0)
          + " must be less than normalizedCutoffFrequencies.second = "
          + std::to_string(w1));
    }
}

[[nodiscard]] UFR::ZerosPolesGain<double>
createDigitalLowpassOrHighpass(
    const int order,
    const double normalizedCutoffFrequency,
    const double ripple,
    const ::Prototype prototype,
    const bool isLowpass)
{
    if (order < 1)
    {
        throw std::invalid_argument("order = " + std::to_string(order)
                                  + " must be positive");
    }
    if (normalizedCutoffFrequency < 0)
    {
        throw std::invalid_argument("normalizedCutoffFrequency = "
                                  + std::to_string(normalizedCutoffFrequency)
                                  + " must be non-negative");
    }
    if (normalizedCutoffFrequency > 1)
    {   
        throw std::invalid_argument("normalizedCutoffFrequency = "
                                  + std::to_string(normalizedCutoffFrequency)
                                  + " cannot exceed 1");
    }   
    if (ripple <= 0 &&
        (prototype == ::Prototype::ChebyshevI ||
         prototype == ::Prototype::ChebyshevII))
    {   
        throw std::invalid_argument("ripple = " + std::to_string(ripple)
                                  + " must be positive");
    }
    std::unique_ptr<UFR::ZerosPolesGain<double>> analogPrototype{nullptr};
    if (prototype == ::Prototype::Bessel)
    {
        analogPrototype
           = std::make_unique<UFR::ZerosPolesGain<double>>
             (AnalogPrototype::bessel(order));
    }
    else if (prototype == ::Prototype::Butterworth)
    {
        analogPrototype
           = std::make_unique<UFR::ZerosPolesGain<double>>
             (AnalogPrototype::butterworth(order));
    }
    else if (prototype == ::Prototype::ChebyshevI)
    {
        analogPrototype
           = std::make_unique<UFR::ZerosPolesGain<double>>
             (AnalogPrototype::chebyshevTypeI(order, ripple));
    }
    else if (prototype == ::Prototype::ChebyshevII)
    {
        analogPrototype
           = std::make_unique<UFR::ZerosPolesGain<double>>
             (AnalogPrototype::chebyshevTypeII(order, ripple));
    }
    else
    {
#ifndef NDEBUG
        assert(false);
#else
        throw std::runtime_error("Unhandled prototype");
#endif
    }
    auto warpedFrequency = ::warpDigital(normalizedCutoffFrequency);
    if (isLowpass)
    {
        auto zpkAnalogLowpass
            = convertLowpassAnalogPrototypeToLowpass(*analogPrototype,
                                                      warpedFrequency);
        return bilinearTransform<double>(zpkAnalogLowpass, 2);
    }
    else
    {
        auto zpkAnalogHighpass
            = convertLowpassAnalogPrototypeToHighpass(*analogPrototype,
                                                      warpedFrequency);
        return bilinearTransform<double>(zpkAnalogHighpass, 2); 
    }
}

[[nodiscard]] UFR::ZerosPolesGain<double>
createDigitalBandpassOrBandstop(
    const int order,
    const std::pair<double, double> &normalizedCutoffFrequencies,
    const double ripple,
    const ::Prototype prototype,
    const bool isBandpass)
{
    if (order < 1)
    {
        throw std::invalid_argument("order = " + std::to_string(order)
                                  + " must be positive");
    }
    ::checkFrequencies(normalizedCutoffFrequencies);
    if (ripple <= 0 &&
        (prototype == ::Prototype::ChebyshevI ||
         prototype == ::Prototype::ChebyshevII))
    {   
        throw std::invalid_argument("ripple = " + std::to_string(ripple)
                                  + " must be positive");
    }
    std::unique_ptr<UFR::ZerosPolesGain<double>> analogPrototype{nullptr};
    if (prototype == ::Prototype::Bessel)
    {
        analogPrototype
           = std::make_unique<UFR::ZerosPolesGain<double>>
             (AnalogPrototype::bessel(order));
    }
    else if (prototype == ::Prototype::Butterworth)
    {
        analogPrototype
           = std::make_unique<UFR::ZerosPolesGain<double>>
             (AnalogPrototype::butterworth(order));
    }
    else if (prototype == ::Prototype::ChebyshevI)
    {
        analogPrototype
           = std::make_unique<UFR::ZerosPolesGain<double>>
             (AnalogPrototype::chebyshevTypeI(order, ripple));
    }
    else if (prototype == ::Prototype::ChebyshevII)
    {
        analogPrototype
           = std::make_unique<UFR::ZerosPolesGain<double>>
             (AnalogPrototype::chebyshevTypeII(order, ripple));
    }
    else
    {
#ifndef NDEBUG
        assert(false);
#else
        throw std::runtime_error("Unhandled prototype");
#endif
    }
    auto warpedFrequencies = ::warpDigital(normalizedCutoffFrequencies);
    double w0 = std::sqrt(warpedFrequencies.first*warpedFrequencies.second);
    double bandwidth = warpedFrequencies.second - warpedFrequencies.first;
    std::pair<double, double> passband(w0, w0 + bandwidth);
    if (isBandpass)
    {
        auto zpkAnalogBandpass
            = convertLowpassAnalogPrototypeToBandpass(*analogPrototype,
                                                      passband);
        return bilinearTransform<double>(zpkAnalogBandpass, 2);
    }
    else
    {
        auto zpkAnalogBandstop
            = convertLowpassAnalogPrototypeToBandstop(*analogPrototype,
                                                      passband);
        return bilinearTransform<double>(zpkAnalogBandstop, 2); 
    }
}

}

///--------------------------------------------------------------------------///
///                                Lowpass                                   ///
///--------------------------------------------------------------------------///
UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createButterworthLowpass(
    int order,
    const double normalizedCutoffFrequency)
{
    constexpr bool isLowpass{true};
    return ::createDigitalLowpassOrHighpass(order,
                                            normalizedCutoffFrequency,
                                            3.0,
                                            ::Prototype::Butterworth,
                                            isLowpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createBesselLowpass(
    int order,
    const double normalizedCutoffFrequency)
{
    constexpr bool isLowpass{true};
    return ::createDigitalLowpassOrHighpass(order,
                                            normalizedCutoffFrequency,
                                            3.0,
                                            ::Prototype::Bessel,
                                            isLowpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createChebyshevTypeILowpass(
    int order,
    const double ripple,
    const double normalizedCutoffFrequency)
{
    constexpr bool isLowpass{true};
    return ::createDigitalLowpassOrHighpass(order,
                                            normalizedCutoffFrequency,
                                            ripple,
                                            ::Prototype::ChebyshevI,
                                            isLowpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createChebyshevTypeIILowpass(
    int order,
    const double ripple,
    const double normalizedCutoffFrequency)
{
    constexpr bool isLowpass{true};
    return ::createDigitalLowpassOrHighpass(order,
                                            normalizedCutoffFrequency,
                                            ripple,
                                            ::Prototype::ChebyshevII,
                                            isLowpass);
}

///--------------------------------------------------------------------------///
///                                Highpass                                  ///
///--------------------------------------------------------------------------///
UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createButterworthHighpass(
    int order,
    const double normalizedCutoffFrequency)
{
    constexpr bool isLowpass{false};
    return ::createDigitalLowpassOrHighpass(order,
                                            normalizedCutoffFrequency,
                                            3.0,
                                            ::Prototype::Butterworth,
                                            isLowpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createBesselHighpass(
    int order,
    const double normalizedCutoffFrequency)
{
    constexpr bool isLowpass{false};
    return ::createDigitalLowpassOrHighpass(order,
                                            normalizedCutoffFrequency,
                                            3.0,
                                            ::Prototype::Bessel,
                                            isLowpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createChebyshevTypeIHighpass(
    int order,
    const double ripple,
    const double normalizedCutoffFrequency)
{
    constexpr bool isLowpass{false};
    return ::createDigitalLowpassOrHighpass(order,
                                            normalizedCutoffFrequency,
                                            ripple,
                                            ::Prototype::ChebyshevI,
                                            isLowpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createChebyshevTypeIIHighpass(
    int order,
    const double ripple,
    const double normalizedCutoffFrequency)
{
    constexpr bool isLowpass{false};
    return ::createDigitalLowpassOrHighpass(order,
                                            normalizedCutoffFrequency,
                                            ripple,
                                            ::Prototype::ChebyshevII,
                                            isLowpass);
}

///--------------------------------------------------------------------------///
///                                Bandpass                                  ///
///--------------------------------------------------------------------------///
UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createButterworthBandpass(
    int order,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    constexpr bool isBandpass{true};
    return ::createDigitalBandpassOrBandstop(order,
                                             normalizedCutoffFrequencies,
                                             3.0,
                                             ::Prototype::Butterworth,
                                             isBandpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createBesselBandpass(
    int order,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    constexpr bool isBandpass{true};
    return ::createDigitalBandpassOrBandstop(order,
                                             normalizedCutoffFrequencies,
                                             3.0,
                                             ::Prototype::Bessel,
                                             isBandpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createChebyshevTypeIBandpass(
    int order,
    const double ripple,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    constexpr bool isBandpass{true};
    return ::createDigitalBandpassOrBandstop(order,
                                             normalizedCutoffFrequencies,
                                             ripple,
                                             ::Prototype::ChebyshevI,
                                             isBandpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createChebyshevTypeIIBandpass(
    int order,
    const double ripple,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    constexpr bool isBandpass{true};
    return ::createDigitalBandpassOrBandstop(order,
                                             normalizedCutoffFrequencies,
                                             ripple,
                                             ::Prototype::ChebyshevII,
                                             isBandpass);
}

///--------------------------------------------------------------------------///
///                                Bandstop                                  ///
///--------------------------------------------------------------------------///
UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createButterworthBandstop(
    int order,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    constexpr bool isBandpass{false};
    return ::createDigitalBandpassOrBandstop(order,
                                             normalizedCutoffFrequencies,
                                             3.0,
                                             ::Prototype::Butterworth,
                                             isBandpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createBesselBandstop(
    int order,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    constexpr bool isBandpass{false};
    return ::createDigitalBandpassOrBandstop(order,
                                             normalizedCutoffFrequencies,
                                             3.0,
                                             ::Prototype::Bessel,
                                             isBandpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createChebyshevTypeIBandstop(
    int order,
    const double ripple,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    constexpr bool isBandpass{false};
    return ::createDigitalBandpassOrBandstop(order,
                                             normalizedCutoffFrequencies,
                                             ripple,
                                             ::Prototype::ChebyshevI,
                                             isBandpass);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createChebyshevTypeIIBandstop(
    int order,
    const double ripple,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    constexpr bool isBandpass{false};
    return ::createDigitalBandpassOrBandstop(order,
                                             normalizedCutoffFrequencies,
                                             ripple,
                                             ::Prototype::ChebyshevII,
                                             isBandpass);
}

