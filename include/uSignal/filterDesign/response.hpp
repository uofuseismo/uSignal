#ifndef USIGNAL_FILTER_DESIGN_RESPONSE_HPP
#define USIGNAL_FILTER_DESIGN_RESPONSE_HPP
#include <uSignal/vector.hpp>
#include <complex>
namespace USignal::FilterRepresentations
{
  template<class T> class InfiniteImpulseResponse;
  template<class T> class FiniteImpulseResponse;
}
namespace USignal::FilterDesign::Response
{

/// @brief Comptutes the frequency response of an analog filter.
/// @param[in] filter       The analog filter infinite impulse response filter.
/// @param[in] frequencies  The angular frequencies, in rad/s, at which to
///                         compute the response.
/// @result The complex valued spectrum that defines the analog filter response
///         at each frequency.  This has the same length as frequencies.
template<typename T>
USignal::Vector<std::complex<T>> 
    computeAnalog(const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filter,
                  const USignal::Vector<T> &frequencies);
/// @brief Computes the frequency response of a digital filter.
/// @param[in] filter       The digital filter infinite impulse response filter.
/// @param[in] frequencies  The normalized angular frequencies at which to 
///                         compute the response.  Note, values should be in
///                         the range \f$ [0, \pi] \f$ where 0 is the zero
///                         frequency and \f$ \pi \f$ is the Nyquist freuqency.
/// @result The complex valued spectrum that defines the digital filter response
///         at each frequency.  This has the same length as frequencies.
template<typename T>
USignal::Vector<std::complex<T>>
    computeDigital(const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filter,
                   const USignal::Vector<T> &frequencies);

}
#endif
