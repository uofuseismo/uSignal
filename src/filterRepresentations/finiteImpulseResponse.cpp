#include <iostream>
#include <complex>
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/vector.hpp"

using namespace USignal::FilterRepresentations;

template<class T>
class FiniteImpulseResponse<T>::FiniteImpulseResponseImpl
{
public:
    class USignal::Vector<T> mFilterCoefficients;
};

/// Constructor
template<class T>
FiniteImpulseResponse<T>::FiniteImpulseResponse(
    const USignal::Vector<T> &filterCoefficients)
{
    constexpr T zero{0};
    if (filterCoefficients.empty())
    {
        throw std::invalid_argument("The filter coefficients are empty");
    }
    bool allZero{true};
    for (const auto &b : filterCoefficients)
    {
        if (b != zero)
        {
            allZero = false;
            break;
        }
    }
    if (allZero)
    {   
        pImpl = nullptr;
        throw std::runtime_error("The filter coefficients are all zero");
    }
    pImpl = std::make_unique<FiniteImpulseResponseImpl> ();
    pImpl->mFilterCoefficients = filterCoefficients;
/*
    int nCoefficients = static_cast<int> (pImpl->mFilterCoefficients.size());
    for (int i = 0; i < nCoefficients; ++i)
    {
        if (pImpl->mFilterCoefficients.at(i) != zero){break;}
        pImpl->mFilterCoefficients.pop_front();
    }
    for (int i = nCoefficients - 1; i >= 0; --i)
    {
        if (pImpl->mFilterCoefficients.at(i) != zero){break;}
        pImpl->mFilterCoefficients.pop_back(); 
    }
    if (pImpl->mFilterCoefficients.empty())
    {
        pImpl = nullptr;
        throw std::runtime_error("The filter coefficients are all zero");
    }
*/
}

/// Copy constructor
template<class T>
FiniteImpulseResponse<T>::FiniteImpulseResponse(
    const FiniteImpulseResponse<T> &firFilter)
{
    *this = firFilter;
}

/// Move constructor
template<class T>
FiniteImpulseResponse<T>::FiniteImpulseResponse(
    FiniteImpulseResponse<T> &&firFilter) noexcept
{
    *this = std::move(firFilter);
}

/// Copy assignment
template<class T>
FiniteImpulseResponse<T> &
FiniteImpulseResponse<T>::operator=(const FiniteImpulseResponse<T> &firFilter)
{
    if (&firFilter == this){return *this;}
    pImpl = std::make_unique<FiniteImpulseResponseImpl> (*firFilter.pImpl);
    return *this;
}

/// Move assignment
template<class T>
FiniteImpulseResponse<T> &
FiniteImpulseResponse<T>::operator=(
    FiniteImpulseResponse<T> &&firFilter) noexcept
{
    if (&firFilter == this){return *this;}
    pImpl = std::move(firFilter.pImpl);
    return *this;
}

/// Filter coefficients
template<class T>
const USignal::Vector<T> &
FiniteImpulseResponse<T>::getFilterCoefficientsReference() const noexcept
{
    return *&pImpl->mFilterCoefficients;
}

/// Filter coefficients
template<class T>
USignal::Vector<T> 
FiniteImpulseResponse<T>::getFilterCoefficients() const noexcept
{
    return pImpl->mFilterCoefficients;
}

/// Filter order
template<class T>
int FiniteImpulseResponse<T>::getOrder() const noexcept
{
    return static_cast<int> (pImpl->mFilterCoefficients.size() - 1);
}

/// Destructor
template<class T>
FiniteImpulseResponse<T>::~FiniteImpulseResponse() = default;

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class USignal::FilterRepresentations::FiniteImpulseResponse<double>;
template class USignal::FilterRepresentations::FiniteImpulseResponse<float>;
template class USignal::FilterRepresentations::FiniteImpulseResponse<int>;
template class USignal::FilterRepresentations::FiniteImpulseResponse<std::complex<double>>;
template class USignal::FilterRepresentations::FiniteImpulseResponse<std::complex<float>>;

