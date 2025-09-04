#include <complex>
#include <cmath>
#include "uSignal/filterDesign/infiniteImpulseResponse/bilinearTransform.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"

namespace
{

template<typename T>
USignal::Vector<std::complex<T>> 
    bilinearTransform(const USignal::Vector<std::complex<T>> &x,
                      const std::complex<T> fs2)
{
    USignal::Vector<std::complex<T>> xBilinear(x.size());
    for (int i = 0; i < static_cast<int> (x.size()); ++i)
    {
        auto numerator = static_cast<std::complex<double>> (fs2 + x[i]);
        auto denominator = static_cast<std::complex<double>> (fs2 - x[i]);
        xBilinear[i] = static_cast<std::complex<T>> (numerator/denominator);
    }
    return xBilinear;
}

template<typename T>
std::complex<double> computeGain(const USignal::Vector<std::complex<T>> &x,
                                 const std::complex<T> fs2)
{
    // Always do this in double precision to mitigate overflow
    std::complex<double> gain{1 + 0i};
    for (const auto &xi : x)
    {
        gain = gain*( static_cast<std::complex<double>> (fs2)
                    - static_cast<std::complex<double>> (xi));
    }
    return static_cast<std::complex<T>> (gain);
}

}

template<typename T>
USignal::FilterRepresentations::ZerosPolesGain<T>
USignal::FilterDesign::InfiniteImpulseResponse::bilinearTransform(
    const USignal::FilterRepresentations::ZerosPolesGain<T> &analogZerosPolesGain,
    const T samplingFrequency)
{
    if (samplingFrequency <= 0)
    {
        throw std::invalid_argument("Sampling frequency must be positive");
    }
    const auto analogZeros = analogZerosPolesGain.getZerosReference();
    const auto analogPoles = analogZerosPolesGain.getPolesReference();
    auto analogGain = static_cast<double> (analogZerosPolesGain.getGain());
    if (analogZeros.size() > analogPoles.size())
    {
        throw std::invalid_argument("Cannot have more zeros than poles");
    }
    // Bilinear transform
    std::complex<T> fs2 = 2*samplingFrequency + 0i;
    auto bilinearZeros = ::bilinearTransform(analogZeros, fs2);
    auto bilinearPoles = ::bilinearTransform(analogPoles, fs2);
    // Zeros at infinite get moved to Nyquist frequency
    for (size_t i = bilinearZeros.size(); i < bilinearPoles.size(); ++i)
    {
        bilinearZeros.push_back(-1 + 0i);
    }
    // Compensate for gain change
    auto analogNumeratorGain = ::computeGain(analogZeros, fs2);
    auto analogDenominatorGain = ::computeGain(analogPoles, fs2);
    auto bilinearGain = analogGain
                       *std::real(analogNumeratorGain/analogDenominatorGain);
    // Pack up result and leave
    USignal::FilterRepresentations::ZerosPolesGain<T>
        digitalZerosPolesGain{bilinearZeros,
                              bilinearPoles,
                              static_cast<T> (bilinearGain)};
    return digitalZerosPolesGain;
}

///--------------------------------------------------------------------------///
///                              Instantiation                               ///
///--------------------------------------------------------------------------///
template
USignal::FilterRepresentations::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::bilinearTransform<double>(
    const USignal::FilterRepresentations::ZerosPolesGain<double> &,
    const double );
template
USignal::FilterRepresentations::ZerosPolesGain<float>
USignal::FilterDesign::InfiniteImpulseResponse::bilinearTransform<float>(
    const USignal::FilterRepresentations::ZerosPolesGain<float> &,
    const float );
