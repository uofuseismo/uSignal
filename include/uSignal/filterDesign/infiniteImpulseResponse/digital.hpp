#ifndef USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_DIGITAL_HPP
#define USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_DIGITAL_HPP
#include <uSignal/filterRepresentations/zerosPolesGain.hpp>
namespace USignal::FilterDesign::InfiniteImpulseResponse::Digital
{

/// Lowpass
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createButterworthLowpass(
       int order,
       double normalizedCutoffFrequency);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createBesselLowpass(
       int order,
       double normalizedCutoffFrequency);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createChebyshevTypeILowpass(
       int order,
       double ripple,
       double normalizedCutoffFrequency);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createChebyshevTypeIILowpass(
       int order,
       double ripple,
       double normalizedCutoffFrequency);

/// Highpass
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createButterworthHighpass(
       int order,
       double normalizedCutoffFrequency);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createBesselHighpass(
       int order,
       double normalizedCutoffFrequency);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createChebyshevTypeIHighpass(
       int order,
       double ripple,
       double normalizedCutoffFrequency);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createChebyshevTypeIIHighpass(
       int order,
       double ripple,
       double normalizedCutoffFrequency);

/// Bandpass
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createButterworthBandpass(
       int order,
       const std::pair<double, double> &normalizedCutoffFrequencies);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createBesselBandpass(
       int order,
       const std::pair<double, double> &normalizedCutoffFrequencies);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createChebyshevTypeIBandpass(
       int order,
       double ripple,
       const std::pair<double, double> &normalizedCutoffFrequencies);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createChebyshevTypeIIBandpass(
       int order,
       double ripple,
       const std::pair<double, double> &normalizedCutoffFrequencies);

/// Bandstop
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createButterworthBandstop(
       int order,
       const std::pair<double, double> &normalizedCutoffFrequencies);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createBesselBandstop(
       int order,
       const std::pair<double, double> &normalizedCutoffFrequencies);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createChebyshevTypeIBandstop(
       int order,
       double ripple,
       const std::pair<double, double> &normalizedCutoffFrequencies);
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double> 
    createChebyshevTypeIIBandstop(
       int order,
       double ripple,
       const std::pair<double, double> &normalizedCutoffFrequencies);

}
#endif
