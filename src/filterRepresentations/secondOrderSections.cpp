#include <array>
#include <vector>
#include "uSignal/filterRepresentations/secondOrderSections.hpp"
#include "uSignal/vector.hpp"

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

