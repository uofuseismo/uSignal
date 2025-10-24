#ifndef USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_ANALOG_PROTOTYPE_HPP
#define USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_ANALOG_PROTOTYPE_HPP
#include <uSignal/filterRepresentations/zerosPolesGain.hpp>

namespace USignal::FilterDesign::InfiniteImpulseResponse::AnalogPrototype
{

/// @brief Designs a Butterworth filter with the given order.
/// @param[in] order  The filter order.  This must be in the range [0, 25].
USignal::FilterRepresentations::ZerosPolesGain<double> butterworth(int order);


}

#endif
