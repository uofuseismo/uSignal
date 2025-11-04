#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include "uSignal/filterRepresentations/secondOrderSections.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include "zpk2sos.hpp"

using namespace USignal::FilterRepresentations;

namespace
{
template<typename T>
USignal::Vector<T> reshape(const std::vector<std::array<T, 3>> &x,
                           const bool isNumerator)
{
    if (x.empty())
    {
        if (isNumerator)
        {
            throw std::invalid_argument("Numerator coefficient array empty");
        }
        else
        {
            throw std::invalid_argument("Denominator coefficient array empty");
        }
    }
    USignal::Vector<T> result;
    result.reserve(3*x.size());
    for (int section = 0; section < static_cast<int> (x.size()); ++section)
    {
        auto coefficient0 = x[section][0];
        if (coefficient0 == 0)
        {
            if (isNumerator)
            {
                throw std::invalid_argument(
                   "Numerator Leading coefficient zero for section "
                  + std::to_string(section + 1));
            }
            else
            {
                throw std::invalid_argument(
                   "Denominator leading coefficient zero for section "
                  + std::to_string(section + 1));
            }
        }
        result.push_back(coefficient0);
        result.push_back(x[section][1]);
        result.push_back(x[section][2]);
    }
    return result;
}

SecondOrderSections<double> 
zpk2sos(const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
        SecondOrderSections<double>::PairingStrategy strategy)
{
    auto isNearest = (strategy == SecondOrderSections<double>::PairingStrategy::Nearest) ? true : false;
    auto gain = static_cast<double> (zpk.getGain());
    auto zerosRef = zpk.getZerosReference(); 
    auto polesRef = zpk.getPolesReference();
    std::vector<std::complex<double>> zeros; zeros.reserve(zerosRef.size());
    for (const auto &z : zerosRef){zeros.push_back(z);}
    std::vector<std::complex<double>> poles; poles.reserve(polesRef.size());
    for (const auto &p : polesRef){poles.push_back(p);}
    auto [bsAll, asAll] = ::zpk2sos(zeros, poles, gain, isNearest);
    SecondOrderSections<double> sos{bsAll, asAll};
    return sos;
}

}

template<typename T>
class SecondOrderSections<T>::SecondOrderSectionsImpl
{
public:
    USignal::Vector<T> mNumeratorCoefficients;
    USignal::Vector<T> mDenominatorCoefficients;
    int mSections{0};
};

/// Destructor
template<typename T>
SecondOrderSections<T>::~SecondOrderSections() = default;

/// Constructor
template<typename T>
SecondOrderSections<T>::SecondOrderSections(
    const std::vector<std::array<T, 3>> &numeratorCoefficients,
    const std::vector<std::array<T, 3>> &denominatorCoefficients)
{
    auto b = ::reshape(numeratorCoefficients, true);
    auto a = ::reshape(denominatorCoefficients, false);
    SecondOrderSections<T> sos{b, a};
    *this = std::move(sos);
}

template<>
SecondOrderSections<double>::SecondOrderSections(
    const InfiniteImpulseResponse<double> &ba,
    const SecondOrderSections::PairingStrategy strategy)
{
    ZerosPolesGain<double> zpk{ba};
    *this = ::zpk2sos(zpk, strategy);
}

template<>
SecondOrderSections<float>::SecondOrderSections(
    const InfiniteImpulseResponse<float> &ba,
    const SecondOrderSections::PairingStrategy strategy)
{
    ZerosPolesGain<float> zpk{ba};
    const auto gain = static_cast<double> (zpk.getGain());
    auto zerosFloat = zpk.getZeros();
    auto polesFloat = zpk.getPoles();
    USignal::Vector<std::complex<double>> zeros(zerosFloat.size());
    for (int i = 0; i < static_cast<int> (zerosFloat.size()); ++i)
    {   
        zeros[i] = static_cast<std::complex<double>> (zerosFloat[i]);
    }   
    USignal::Vector<std::complex<double>> poles(polesFloat.size());
    for (int i = 0; i < static_cast<int> (polesFloat.size()); ++i)
    {   
        poles[i] = static_cast<std::complex<double>> (polesFloat[i]);
    }
    ZerosPolesGain<double> zpkDouble(zeros, poles, gain);
    auto strategyDouble
        = static_cast<SecondOrderSections<double>::PairingStrategy> (strategy);
    auto sos = ::zpk2sos(zpkDouble, strategyDouble); 
    auto bs = sos.getNumeratorFilterCoefficients();
    auto as = sos.getDenominatorFilterCoefficients();
    USignal::Vector<float> bsFloat(bs.size());
    for (int i = 0; i < static_cast<int> (bs.size()); ++i)
    {   
        bsFloat[i] = static_cast<float> (bs[i]);
    }   
    USignal::Vector<float> asFloat(as.size());
    for (int i = 0; i < static_cast<int> (as.size()); ++i)
    {   
        asFloat[i] = static_cast<float> (as[i]);
    }   
    *this = SecondOrderSections<float> (bsFloat, asFloat);
}

template<>
SecondOrderSections<double>::SecondOrderSections(
    const ZerosPolesGain<double> &zpk,
    const SecondOrderSections::PairingStrategy strategy)
{
    *this = ::zpk2sos(zpk, strategy); 
}

template<>
SecondOrderSections<float>::SecondOrderSections(
    const ZerosPolesGain<float> &zpk,
    const SecondOrderSections::PairingStrategy strategy)
{
    const auto gain = static_cast<double> (zpk.getGain());
    auto zerosFloat = zpk.getZeros();
    auto polesFloat = zpk.getPoles();
    USignal::Vector<std::complex<double>> zeros(zerosFloat.size());
    for (int i = 0; i < static_cast<int> (zerosFloat.size()); ++i)
    {
        zeros[i] = static_cast<std::complex<double>> (zerosFloat[i]);
    }
    USignal::Vector<std::complex<double>> poles(polesFloat.size());
    for (int i = 0; i < static_cast<int> (polesFloat.size()); ++i)
    {
        poles[i] = static_cast<std::complex<double>> (polesFloat[i]);
    }
    ZerosPolesGain<double> zpkDouble(zeros, poles, gain);
    auto strategyDouble
        = static_cast<SecondOrderSections<double>::PairingStrategy> (strategy);
    auto sos = ::zpk2sos(zpkDouble, strategyDouble); 
    auto bs = sos.getNumeratorFilterCoefficients();
    auto as = sos.getDenominatorFilterCoefficients();
    USignal::Vector<float> bsFloat(bs.size());
    for (int i = 0; i < static_cast<int> (bs.size()); ++i)
    {
        bsFloat[i] = static_cast<float> (bs[i]);
    }
    USignal::Vector<float> asFloat(as.size());
    for (int i = 0; i < static_cast<int> (as.size()); ++i)
    {   
        asFloat[i] = static_cast<float> (as[i]);
    }   
    *this = SecondOrderSections<float> (bsFloat, asFloat);
}


/// Constructor
template<typename T>
SecondOrderSections<T>::SecondOrderSections(
    const USignal::Vector<T> &numeratorCoefficients,
    const USignal::Vector<T> &denominatorCoefficients) :
    pImpl(std::make_unique<SecondOrderSectionsImpl> ())
{
    if (numeratorCoefficients.size() !=
        denominatorCoefficients.size())
    {
        throw std::invalid_argument(
           "numeratorCoefficients.size() != denominatorCoefficients.size()");
    }
    if (numeratorCoefficients.size() < 3)
    {
        throw std::invalid_argument(
           "At least 3 numerator and denominator coefficients required");
    }
    if (numeratorCoefficients.size()%3 != 0)
    {
        throw std::invalid_argument( 
           "Number of coefficients must be multiple of 3");
    }
    int nSections = static_cast<int> (numeratorCoefficients.size())/3;
    for (int section = 0; section < nSections; ++section)
    {
        auto numeratorCoefficient = numeratorCoefficients[3*section];
        if (numeratorCoefficient == 0)
        {
            throw std::invalid_argument(
                "Leading numerator coefficient for "
                + std::to_string(section + 1)
                + " section is zero");
        }
        auto denominatorCoefficient = denominatorCoefficients[3*section];
        if (denominatorCoefficient == 0)
        {
            throw std::invalid_argument(
                "Leading denominator coefficient for "
                + std::to_string(section + 1)
                + " section is zero");
        }
    }
    pImpl->mSections = nSections;
    pImpl->mNumeratorCoefficients = numeratorCoefficients;
    pImpl->mDenominatorCoefficients = denominatorCoefficients; 
}

/// Copy constructor
template<typename T>
SecondOrderSections<T>::SecondOrderSections(const SecondOrderSections<T> &sos)
{
    *this = sos;
}

/// Move constructor
template<typename T>
SecondOrderSections<T>::SecondOrderSections(
    SecondOrderSections<T> &&sos) noexcept
{
    *this = std::move(sos);
}

/// Copy assignment
template<typename T>
SecondOrderSections<T>& SecondOrderSections<T>::operator=(
    const SecondOrderSections &sections)
{
    if (&sections == this){return *this;}
    pImpl = std::make_unique<SecondOrderSectionsImpl> (*sections.pImpl);
    return *this;
}

/// Move assignment
template<typename T>
SecondOrderSections<T>& SecondOrderSections<T>::operator=(
    SecondOrderSections &&sections) noexcept
{
    if (&sections == this){return *this;}
    pImpl = std::move(sections.pImpl);
    return *this;
}

/// Numerator coefficients
template<class T>
const USignal::Vector<T> &
SecondOrderSections<T>::getNumeratorFilterCoefficientsReference() const noexcept
{
    return *&pImpl->mNumeratorCoefficients;
}

template<class T>
USignal::Vector<T> 
SecondOrderSections<T>::getNumeratorFilterCoefficients() const noexcept
{
    return pImpl->mNumeratorCoefficients;
}

/// Numerator coefficients
template<class T>
const USignal::Vector<T> &
SecondOrderSections<T>::getDenominatorFilterCoefficientsReference() const noexcept
{
    return *&pImpl->mDenominatorCoefficients;
}

template<class T>
USignal::Vector<T> 
SecondOrderSections<T>::getDenominatorFilterCoefficients() const noexcept
{
    return pImpl->mDenominatorCoefficients;
}

/// Number of sections
template<class T>
int SecondOrderSections<T>::getNumberOfSections() const noexcept
{
    return pImpl->mSections;
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class USignal::FilterRepresentations::SecondOrderSections<double>;
template class USignal::FilterRepresentations::SecondOrderSections<float>;
template class USignal::FilterRepresentations::SecondOrderSections<int>;

