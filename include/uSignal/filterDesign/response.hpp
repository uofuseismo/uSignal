#ifndef USIGNAL_FILTER_DESIGN_RESPONSE_HPP
#define USIGNAL_FILTER_DESIGN_RESPONSE_HPP
#include <uSignal/vector.hpp>
#include <complex>
#include <cmath>
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
/// @brief Computes the frequency response of a digital filter.
/// @param[in] filter       The digital filter finite impulse response filter.
/// @param[in] frequencies  The normalized angular frequencies at which to 
///                         compute the response.  Note, values should be in
///                         the range \f$ [0, \pi] \f$ where 0 is the zero
///                         frequency and \f$ \pi \f$ is the Nyquist freuqency.
/// @result The complex valued spectrum that defines the digital filter response
///         at each frequency.  This has the same length as frequencies.
template<typename T>
USignal::Vector<std::complex<T>>
    computeDigital(const USignal::FilterRepresentations::FiniteImpulseResponse<T> &filter,
                   const USignal::Vector<T> &frequencies);
/// @brief Computes the magnitude of the frequency response of a digital filter.
/// @param[in] filter       The digital filter infinite impulse response filter.
/// @param[in] frequencies  The normalized angular frequencies at which to 
///                         compute the response.  Note, values should be in
///                         the range \f$ [0, \pi] \f$ where 0 is the zero
///                         frequency and \f$ \pi \f$ is the Nyquist freuqency.
/// @result The magnitude of the of the digital filter response at each
///         frequency.  This has the same length as frequencies.
template<typename T>
USignal::Vector<T>
    computeDigitalAmplitudeSpectrum(const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filter,
                                    const USignal::Vector<T> &frequencies);
/// @brief Computes the magnitude of the frequency response of a digital filter.
/// @param[in] filter       The digital filter finite impulse response filter.
/// @param[in] frequencies  The normalized angular frequencies at which to 
///                         compute the response.  Note, values should be in
///                         the range \f$ [0, \pi] \f$ where 0 is the zero
///                         frequency and \f$ \pi \f$ is the Nyquist freuqency.
/// @result The magnitude of the of the digital filter response at each
///         frequency.  This has the same length as frequencies.
template<typename T>
USignal::Vector<T>
    computeDigitalAmplitudeSpectrum(const USignal::FilterRepresentations::FiniteImpulseResponse<T> &filter,
                                    const USignal::Vector<T> &frequencies);
/// @brief Computes the phase spectra of a digital filter.
/// @param[in] filter       The digital filter infinite impulse response filter.
/// @param[in] frequencies  The normalized angular frequencies at which to 
///                         compute the response.  Note, values should be in
///                         the range \f$ [0, \pi] \f$ where 0 is the zero
///                         frequency and \f$ \pi \f$ is the Nyquist freuqency.
/// @param[in] unwrap       This 
/// @result The unwrapped phase, in radians, of the of the digital filter response
///         at each frequency.  This has the same length as frequencies.
template<typename T>
USignal::Vector<T>
    computeDigitalPhaseSpectrum(const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filter,
                                const USignal::Vector<T> &frequencies);
/// @brief Computes the phase spectra of a digital filter.
/// @param[in] filter       The digital filter finite impulse response filter.
/// @param[in] frequencies  The normalized angular frequencies at which to 
///                         compute the response.  Note, values should be in
///                         the range \f$ [0, \pi] \f$ where 0 is the zero
///                         frequency and \f$ \pi \f$ is the Nyquist freuqency.
/// @param[in] unwrap       This 
/// @result The unwrapped phase, in radians, of the of the digital filter response
///         at each frequency.  This has the same length as frequencies.
template<typename T>
USignal::Vector<T>
    computeDigitalPhaseSpectrum(const USignal::FilterRepresentations::FiniteImpulseResponse<T> &filter,
                                const USignal::Vector<T> &frequencies);

/// @brief Unwraps the phase.
/// @param[in] phases     The phase angles, in radians, tabulated at each
///                       frequency.  Note, this assumes the frequencies were
///                       in order when this was tabulated.
/// @param[in] tolerance  When the difference between consecutive phase angles
///                       exceeds this tolerance, it's phase will be unwrapped
///                       by adding 2\pi to the second phase angle.
/// @result The unwrapped phases.
template<typename T>
USignal::Vector<T> unwrapPhase(const USignal::Vector<T> &phases, const T tolerance = M_PI);

template<typename T>
USignal::Vector<T>
    computeDigitalUnwrappedPhaseSpectrum(const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filter,
                                         const USignal::Vector<T> &frequencies,
                                         const double unwrap);
template<typename T>
USignal::Vector<T>
    computeDigitalUnwrappedPhaseSpectrum(const USignal::FilterRepresentations::FiniteImpulseResponse<T> &filter,
                                         const USignal::Vector<T> &frequencies,
                                         const double unwrap);




}
#endif
