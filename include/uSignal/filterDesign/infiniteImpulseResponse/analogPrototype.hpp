#ifndef USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_ANALOG_PROTOTYPE_HPP
#define USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_ANALOG_PROTOTYPE_HPP
#include <uSignal/filterRepresentations/zerosPolesGain.hpp>

namespace USignal::FilterDesign::InfiniteImpulseResponse::AnalogPrototype
{

/// @brief Designs a Butterworth filter with the given order.
/// @param[in] order  The filter order.  This must be in the range [0, 25].
/// @result The Butterworth analog prototype filter.
USignal::FilterRepresentations::ZerosPolesGain<double> butterworth(int order);

/// @brief Designs a Bessel filter with the given order.
/// @param[in] order  The filter order.  This must be in the range [0, 24].
/// @result The Butterworth analog prototype filter.
USignal::FilterRepresentations::ZerosPolesGain<double> bessel(int order);

/// @brief Designs a Chebyshev Type I filter.  This has ripples in the
///        passband.
/// @param[in] order   The filter order.  The filter will have order + 1 poles.
/// @param[in] ripple  The ripple factor in decibels - i.e., the filter will
///                    have ripple deciebels of ripple in the passband.
/// @result The Chebyshev Type I analog prototype filter.
USignal::FilterRepresentations::ZerosPolesGain<double>
chebyshevTypeI(const int order, const double ripple);

/// @brief Designs a Chebyshev Type II filter.  This has ripples in the
///        stopband.
/// @param[in] order   The filter order.  The filter will have order + 1 poles.
/// @param[in] ripple  The ripple factor in decibels - i.e., the filter will
///                    have ripple deciebels of ripple in the stopband.
/// @result The Chebyshev Type II analog prototype filter.
USignal::FilterRepresentations::ZerosPolesGain<double>
chebyshevTypeII(const int order, const double ripple);



}

#endif
