#ifndef USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_DIGITAL_HPP
#define USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_DIGITAL_HPP
#include <uSignal/filterRepresentations/zerosPolesGain.hpp>
namespace USignal::FilterDesign::InfiniteImpulseResponse::Digital
{

[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createButterworthBandpass(
       int order,
       const std::pair<double, double> &normalizedCutoffFrequencies);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createBesselBandpass(
       int order,
       const std::pair<double, double> &normalizedCutoffFrequencies);

}
#endif
