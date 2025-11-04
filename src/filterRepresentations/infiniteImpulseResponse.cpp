#include <iostream>
#include <complex>
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/vector.hpp"
#include "src/utilities/math/polynomial.hpp"

using namespace USignal::FilterRepresentations;

template<typename T>
class InfiniteImpulseResponse<T>::InfiniteImpulseResponseImpl
{
public:
    USignal::Vector<T> mNumeratorCoefficients;
    USignal::Vector<T> mDenominatorCoefficients;
    int mOrder{0};
};

/// Destructor
template<typename T>
InfiniteImpulseResponse<T>::~InfiniteImpulseResponse() = default;

/// Constructor
template<typename T>
InfiniteImpulseResponse<T>::InfiniteImpulseResponse(
    const USignal::Vector<T> &numeratorCoefficients,
    const USignal::Vector<T> &denominatorCoefficients) :
    pImpl(std::make_unique<InfiniteImpulseResponseImpl> ())
{
    if (numeratorCoefficients.empty())
    {
        throw std::invalid_argument("No numerator coefficients");
    }
    if (denominatorCoefficients.empty())
    {
        throw std::invalid_argument("No denominator coefficients");
    }
    if (denominatorCoefficients.at(0) == 0)
    {
        throw std::invalid_argument(
           "First denominator coefficient cannot be zero");
    }
    pImpl->mNumeratorCoefficients = numeratorCoefficients;
    pImpl->mDenominatorCoefficients = denominatorCoefficients; 
    pImpl->mOrder
        = std::max( static_cast<int> (numeratorCoefficients.size()),
                    static_cast<int> (denominatorCoefficients.size()) )
        - 1;
}

/// Constructor
template<typename T>
InfiniteImpulseResponse<T>::InfiniteImpulseResponse(
    const ZerosPolesGain<T> &zpk) :
    pImpl(std::make_unique<InfiniteImpulseResponseImpl> ())
{
    const auto zeros = zpk.getZerosReference();
    const auto poles = zpk.getPolesReference();
    auto gain = zpk.getGain();
    // Compute polynomial representation of zeros by expanding:
    // (z - z_1)*(z - z_2)*...*(z - z_n)
    auto b = Utilities::Math::Polynomial::expandToRealCoefficients(zeros);
    b = b*gain; // Introduce gain into zeros
    // Compute polynomial representation of poles by expanding:
    // (p - p_1)*(p - p_2)*...*(p - p_n)
    auto a = Utilities::Math::Polynomial::expandToRealCoefficients(poles);
    pImpl->mOrder
        = std::max( static_cast<int> (b.size()),
                    static_cast<int> (a.size()) )
        - 1;
    pImpl->mNumeratorCoefficients = std::move(b);
    pImpl->mDenominatorCoefficients = std::move(a);
}

/*
template<>
InfiniteImpulseResponse<std::complex<double>>::InfiniteImpulseResponse(
    const ZerosPolesGain<std::complex<double>> &zpk) :
    pImpl(std::make_unique<InfiniteImpulseResponseImpl> ()) 
{
    const auto zeros = zpk.getZerosReference();
    const auto poles = zpk.getPolesReference();
    auto gain = zpk.getGain();
    // Compute polynomial representation of zeros by expanding:
    // (z - z_1)*(z - z_2)*...*(z - z_n)
    auto b = Utilities::Math::Polynomial::expand(zeros);
    //b = b*gain; // Introduce gain into zeros
    // Compute polynomial representation of poles by expanding:
    // (p - p_1)*(p - p_2)*...*(p - p_n)
    auto a = Utilities::Math::Polynomial::expand(poles);
    pImpl->mNumeratorCoefficients = b;
    pImpl->mDenominatorCoefficients = a;
    pImpl->mOrder
        = std::max( static_cast<int> (b.size()),
                    static_cast<int> (a.size()) )
        - 1;
}
*/

/// Copy constructor
template<typename T>
InfiniteImpulseResponse<T>::InfiniteImpulseResponse(
    const InfiniteImpulseResponse<T> &iir)
{
    *this = iir;
}

/// Move constructor
template<typename T>
InfiniteImpulseResponse<T>::InfiniteImpulseResponse(
    InfiniteImpulseResponse<T> &&iir) noexcept
{
    *this = std::move(iir);
}

/// Copy assignment
template<typename T>
InfiniteImpulseResponse<T>& InfiniteImpulseResponse<T>::operator=(
    const InfiniteImpulseResponse &sections)
{
    if (&sections == this){return *this;}
    pImpl = std::make_unique<InfiniteImpulseResponseImpl> (*sections.pImpl);
    return *this;
}

/// Move assignment
template<typename T>
InfiniteImpulseResponse<T>& InfiniteImpulseResponse<T>::operator=(
    InfiniteImpulseResponse &&sections) noexcept
{
    if (&sections == this){return *this;}
    pImpl = std::move(sections.pImpl);
    return *this;
}

/// Numerator coefficients
template<class T>
const USignal::Vector<T> &
InfiniteImpulseResponse<T>::getNumeratorFilterCoefficientsReference() const noexcept
{
    return *&pImpl->mNumeratorCoefficients;
}

template<class T>
USignal::Vector<T> 
InfiniteImpulseResponse<T>::getNumeratorFilterCoefficients() const noexcept
{
    return pImpl->mNumeratorCoefficients;
}

/// Numerator coefficients
template<class T>
const USignal::Vector<T> &
InfiniteImpulseResponse<T>::getDenominatorFilterCoefficientsReference() const noexcept
{
    return *&pImpl->mDenominatorCoefficients;
}

template<class T>
USignal::Vector<T> 
InfiniteImpulseResponse<T>::getDenominatorFilterCoefficients() const noexcept
{
    return pImpl->mDenominatorCoefficients;
}

/// Order
template<class T>
int InfiniteImpulseResponse<T>::getOrder() const noexcept
{
    return pImpl->mOrder;
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class USignal::FilterRepresentations::InfiniteImpulseResponse<double>;
template class USignal::FilterRepresentations::InfiniteImpulseResponse<float>;
//template class USignal::FilterRepresentations::InfiniteImpulseResponse<int>;

