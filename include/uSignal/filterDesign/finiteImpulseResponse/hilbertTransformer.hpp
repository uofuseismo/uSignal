#ifndef USIGNAL_FILTER_DESIGN_FINITE_IMPULSE_RESPONSE_HILBERT_TRANSFORMER_HPP
#define USIGNAL_FILTER_DESIGN_FINITE_IMPULSE_RESPONSE_HILBERT_TRANSFORMER_HPP
#include <complex>
#include <uSignal/filterRepresentations/finiteImpulseResponse.hpp>

namespace USignal::FilterDesign::FiniteImpulseResponse
{

///@brief Designs an FIR Hilbert transformer filter.
template<typename T>
USignal::FilterRepresentations::FiniteImpulseResponse<std::complex<T>>
hilbertTransformer(const int order, const double beta);

}
#endif
