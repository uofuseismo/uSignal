#ifndef USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_BILINEAR_TRANSFORM_HPP
#define USIGNAL_FILTER_DESIGN_INFINITE_IMPULSE_RESPONSE_BILINEAR_TRANSFORM_HPP
#include <uSignal/filterRepresentations/zerosPolesGain.hpp>

namespace USignal::FilterDesign::InfiniteImpulseResponse
{
template<typename T>
USignal::FilterRepresentations::ZerosPolesGain<T>
bilinearTransform(const USignal::FilterRepresentations::ZerosPolesGain<T> &analogZerosPolesGain, T samplingFrequency);
};

#endif
