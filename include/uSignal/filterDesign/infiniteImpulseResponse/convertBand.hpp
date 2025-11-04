#ifndef USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_CONVERT_BAND_HPP
#define USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_CONVERT_BAND_HPP
#include <uSignal/filterRepresentations/zerosPolesGain.hpp>
namespace USignal::FilterDesign::InfiniteImpulseResponse
{
/// @brief Converts a lowpass analog prototype filter to a lowpass filter
///        with a different cutoff frequency.
/// @param[in] zpk    The analog lowpass prototype with unity cutoff frequency.
/// @param[in] omega0 The desired angular frequency in rad/s.
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double>
   convertLowpassAnalogPrototypeToLowpass(
      const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
      double omega0);
/// @brief Converts a lowpass analog prototype filter to a highpass filter
///        with a different cutoff frequency.
/// @param[in] zpk    The analog lowpass prototype with unity cutoff frequency.
/// @param[in] omega0 The desired angular frequency in rad/s.
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double>
   convertLowpassAnalogPrototypeToHighpass(
      const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
      double omega0);
/// @brief Converts a lowpass analog prototype filter to a bandpass filter
///        with a different cutoff frequency.
/// @param[in] zpk    The analog lowpass prototype with unity cutoff frequency.
/// @param[in] omega0 The desired passband angular frequencies in rad/s.
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double>
   convertLowpassAnalogPrototypeToBandpass(
      const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
      const std::pair<double, double> &omega0);
/// @brief Converts a lowpass analog prototype filter to a bandstop filter
///        with a different cutoff frequency.
/// @param[in] zpk    The analog lowpass prototype with unity cutoff frequency.
/// @param[in] omega0 The desired stopband angular frequencies in rad/s.
[[nodiscard]] USignal::FilterRepresentations::ZerosPolesGain<double>
   convertLowpassAnalogPrototypeToBandstop(
      const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
      const std::pair<double, double> &omega0);
}
#endif
