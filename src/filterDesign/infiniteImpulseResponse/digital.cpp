#include <cmath>
#include "uSignal/filterDesign/infiniteImpulseResponse/digital.hpp"
#include "uSignal/filterDesign/infiniteImpulseResponse/bilinearTransform.hpp"
#include "uSignal/filterDesign/infiniteImpulseResponse/analogPrototype.hpp"
#include "uSignal/filterDesign/infiniteImpulseResponse/convertBand.hpp"

using namespace USignal::FilterDesign::InfiniteImpulseResponse;
namespace UFR = USignal::FilterRepresentations;

namespace
{
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
                                  + " must be positive");
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
    
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createButterworthBandpass(
    int order,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    if (order < 1)
    {
        throw std::invalid_argument("order = " + std::to_string(order)
                                  + " must be positive");
    }
    ::checkFrequencies(normalizedCutoffFrequencies);
    auto zpkAnalogPrototype = AnalogPrototype::butterworth(order);
    auto warpedFrequencies = ::warpDigital(normalizedCutoffFrequencies);
    auto bandwidth = warpedFrequencies.second - warpedFrequencies.first;
    auto w0 = std::sqrt(warpedFrequencies.first*warpedFrequencies.second);
    std::pair<double, double> passband(w0, w0 + bandwidth);
    auto zpkAnalogBandpass
        = convertLowpassAnalogPrototypeToBandpass(zpkAnalogPrototype,
                                                  warpedFrequencies);
    
    return bilinearTransform<double>(zpkAnalogBandpass, 2);
}

UFR::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::Digital::createBesselBandpass(
    int order,
    const std::pair<double, double> &normalizedCutoffFrequencies)
{
    if (order < 1)
    {
        throw std::invalid_argument("order = " + std::to_string(order)
                                  + " must be positive");
    }
    ::checkFrequencies(normalizedCutoffFrequencies);
    auto zpkAnalogPrototype = AnalogPrototype::bessel(order);
    auto warpedFrequencies = ::warpDigital(normalizedCutoffFrequencies);
    auto bandwidth = warpedFrequencies.second - warpedFrequencies.first;
    auto w0 = std::sqrt(warpedFrequencies.first*warpedFrequencies.second);
    std::pair<double, double> passband(w0, w0 + bandwidth);
    auto zpkAnalogBandpass
        = convertLowpassAnalogPrototypeToBandpass(zpkAnalogPrototype,
                                                  warpedFrequencies);
    return bilinearTransform<double>(zpkAnalogBandpass, 2);
}

