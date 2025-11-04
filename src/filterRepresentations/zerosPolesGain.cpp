#include <iostream>
#include <complex>
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include "src/utilities/math/polynomial.hpp"

using namespace USignal::FilterRepresentations;

template<class T>
class ZerosPolesGain<T>::ZerosPolesGainImpl
{
public:
    USignal::Vector<std::complex<T>> mZeros;
    USignal::Vector<std::complex<T>> mPoles;
    T mGain{1};
    bool mInitialized{false};
};

/// Constructor
template<class T>
ZerosPolesGain<T>::ZerosPolesGain(
    const USignal::Vector<std::complex<T>> &zeros,
    const USignal::Vector<std::complex<T>> &poles,
    const T gain) :
    pImpl(std::make_unique<ZerosPolesGainImpl> ())
{
    if (gain == 0){throw std::invalid_argument("Gain is zero");}
    pImpl->mZeros = zeros;
    pImpl->mPoles = poles;
    pImpl->mGain = gain;
    pImpl->mInitialized = true;
}

/// Construct from a transfer function IIR representations
template<class T>
ZerosPolesGain<T>::ZerosPolesGain(
    const InfiniteImpulseResponse<T> &ba)
{
    auto bs = ba.getNumeratorFilterCoefficients();
    if (bs.empty())
    {
        throw std::invalid_argument("No numerator coefficients");
    }
    if (bs.at(0) == 0){throw std::invalid_argument("b[0] == 0");}
    auto as = ba.getDenominatorFilterCoefficients();
    if (as.empty())
    {
        throw std::invalid_argument("No denominator coefficients");
    }
    const T a0{as.at(0)};
    if (as.at(0) == 0){throw std::invalid_argument("a[0] == 0");}
    // Normalize
    bs = bs/a0;
    as = as/a0;
    // Compute gain
    const T gain{bs.at(0)};
    bs = bs/gain;
    // Get the numerator/denominator roots
    USignal::Vector<std::complex<T>> 
        zeros( std::vector<std::complex<T>> (1, 0) );
    if (bs.size() > 1)
    {
        zeros = USignal::Utilities::Math::Polynomial::computeRoots(bs);
    }
    USignal::Vector<std::complex<T>>
        poles( std::vector<std::complex<T>> (1, 0) );
    if (as.size() > 1)
    {
        poles = USignal::Utilities::Math::Polynomial::computeRoots(as);
    }
    ZerosPolesGain<T> zpk(zeros, poles, gain);
    *this = std::move(zpk);
}

/// Copy constructor
template<class T>
ZerosPolesGain<T>::ZerosPolesGain(const ZerosPolesGain<T> &zpk)
{
    *this = zpk;
}

/// Move constructor
template<class T>
ZerosPolesGain<T>::ZerosPolesGain(ZerosPolesGain<T> &&zpk) noexcept
{   
    *this = std::move(zpk);
}   

/// Copy assignment
template<class T>
ZerosPolesGain<T> &
ZerosPolesGain<T>::operator=(const ZerosPolesGain &zpk)
{
    if (&zpk == this){return *this;}
    pImpl = std::make_unique<ZerosPolesGainImpl> (*zpk.pImpl);
    return *this;
}

/// Move assignment
template<class T> 
ZerosPolesGain<T> &
ZerosPolesGain<T>::operator=(ZerosPolesGain<T> &&zpk) noexcept
{   
    if (&zpk == this){return *this;}
    pImpl = std::move(zpk.pImpl);
    return *this;
}

/// The poles
template<class T>
USignal::Vector<std::complex<T>> ZerosPolesGain<T>::getPoles() const
{
    if (!pImpl->mInitialized)
    {
        throw std::runtime_error("System not initialized - check gain");
    }
    return pImpl->mPoles;
}

template<class T>
const USignal::Vector<std::complex<T>> 
&ZerosPolesGain<T>::getPolesReference() const
{
    if (!pImpl->mInitialized)
    {
        throw std::runtime_error("System not initialized - check gain");
    }
    return *&pImpl->mPoles;
}

/// The zeros
template<class T>
USignal::Vector<std::complex<T>> ZerosPolesGain<T>::getZeros() const
{
    if (!pImpl->mInitialized)
    {
        throw std::runtime_error("System not initialized - check gain");
    }
    return pImpl->mZeros;
}

template<class T>
const USignal::Vector<std::complex<T>> 
    &ZerosPolesGain<T>::getZerosReference() const
{
    if (!pImpl->mInitialized)
    {
        throw std::runtime_error("System not initialized - check gain");
    }
    return *&pImpl->mZeros;
}

/// The gain
template<class T> T ZerosPolesGain<T>::getGain() const
{
    if (!pImpl->mInitialized)
    {
        throw std::runtime_error("System not initialized - check gain");
    }
    return pImpl->mGain;
}

/// Destructor
template<class T>
ZerosPolesGain<T>::~ZerosPolesGain() = default;

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class USignal::FilterRepresentations::ZerosPolesGain<double>;
template class USignal::FilterRepresentations::ZerosPolesGain<float>;

