#include <complex>
#include "uSignal/system/system.hpp"
#include "uSignal/vector.hpp"

using namespace USignal::System;

template<class U, class T>
class ISystem<U, T>::ISystemImpl
{
public:
    USignal::Vector<U> mX;
    USignal::Vector<T> mY;
};

/// Constructor
template<class U, class T>
ISystem<U, T>::ISystem() :
    pImpl(std::make_unique<ISystemImpl> ())
{
}

/// Copy constructor
template<class U, class T>
ISystem<U, T>::ISystem(const ISystem<U, T> &system)
{
    *this = system;
}

/// Move constructor
template<class U, class T>
ISystem<U, T>::ISystem(ISystem<U, T> &&system) noexcept
{
    *this = std::move(system);
}

/// Copy assignment
template<class U, class T>
ISystem<U, T>& ISystem<U, T>::operator=(const ISystem<U, T> &system)
{
    if (&system == this){return *this;}
    pImpl = std::make_unique<ISystemImpl> (*system.pImpl);
    return *this;
}

/// Move assignment
template<class U, class T>
ISystem<U, T>& ISystem<U, T>::operator=(ISystem<U, T> &&system) noexcept
{
    if (&system == this){return *this;}
    pImpl = std::move(system.pImpl);
    return *this;
}

/// Set the input signal
template<class U, class T>
void ISystem<U, T>::setInput(const USignal::Vector<U> &x)
{
    pImpl->mX = x;
}

template<class U, class T>
void ISystem<U, T>::setInput(USignal::Vector<U> &&x) noexcept
{
    pImpl->mX = std::move(x);
}

/// Get the input signal
template<class U, class T>
USignal::Vector<U> ISystem<U, T>::getInput() const
{
    return pImpl->mX;
}

template<class U, class T>
const USignal::Vector<U>& ISystem<U, T>::getInputReference() const noexcept
{
    return pImpl->mX;
}

/// Reset class
template<class U, class T> void ISystem<U, T>::clear() noexcept
{
    pImpl = std::make_unique<ISystemImpl> ();
}

/// Destructor
template<class U, class T> ISystem<U, T>::~ISystem() = default;

/// Set the output signal
template<class U, class T>
void ISystem<U, T>::setOutput(const USignal::Vector<T> &y)
{
    pImpl->mY = y;
}

template<class U, class T>
void ISystem<U, T>::setOutput(USignal::Vector<T> &&y) noexcept
{
    pImpl->mY = std::move(y);
}

/// Get the output signal
template<class U, class T>
USignal::Vector<T> ISystem<U, T>::getOutput() const
{
    return pImpl->mY;
}

template<class U, class T>
const USignal::Vector<T>& ISystem<U, T>::getOutputReference() const noexcept
{
    return pImpl->mY;
} 

///--------------------------------------------------------------------------///
///                             Template Instantiation                       ///
///--------------------------------------------------------------------------///
template class USignal::System::ISystem<double, double>;
template class USignal::System::ISystem<float,  float>;
template class USignal::System::ISystem<float,  double>;
template class USignal::System::ISystem<double, float>;
template class USignal::System::ISystem<double, std::complex<double>>;
template class USignal::System::ISystem<std::complex<double>, double>;
template class USignal::System::ISystem<std::complex<float>, float>;
