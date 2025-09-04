#include <complex>
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/vector.hpp"

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

