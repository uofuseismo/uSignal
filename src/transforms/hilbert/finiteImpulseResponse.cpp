#include "uSignal/transforms/hilbert/finiteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/filterDesign/finiteImpulseResponse/hilbertTransformer.hpp"

using namespace USignal::Transforms::Hilbert;

#ifdef WITH_IPP
#include <ipp/ipps.h>
template<class T>
class FiniteImpulseResponse<T>::FiniteImpulseResponseImpl
{
public:
    FiniteImpulseResponseImpl(const FiniteImpulseResponseOptions &options)
    {
        namespace UFD = USignal::FilterDesign::FiniteImpulseResponse;
        auto filter
            = UFD::hilbertTransformer<double>
              (options.getOrder(), options.getBeta());
    }
    void apply( )
    {
        if (mIsTypeIII)
        {

        }
    }
//private:
    Ipp64f *pSpecReal64f{nullptr};
    bool mIsTypeIII{false};

};
#else
static_assert(false, "Only IPP FIR filter implemented");
#endif

///--------------------------------------------------------------------------///
///                                    Options                               ///
///--------------------------------------------------------------------------///
class FiniteImpulseResponseOptions::FiniteImpulseResponseOptionsImpl
{
public:
    double mBeta{8};
    int mOrder{300};
};

/// Constructor
FiniteImpulseResponseOptions::FiniteImpulseResponseOptions() :
    pImpl(std::make_unique<FiniteImpulseResponseOptionsImpl> ())
{
}

/// Copy constructor
FiniteImpulseResponseOptions::FiniteImpulseResponseOptions(
    const FiniteImpulseResponseOptions &options)
{
    *this = options;
}

/// Move constructor
FiniteImpulseResponseOptions::FiniteImpulseResponseOptions(
    FiniteImpulseResponseOptions &&options) noexcept
{
    *this = options;
}

/// Copy assignment
FiniteImpulseResponseOptions& 
FiniteImpulseResponseOptions::operator=(
    const FiniteImpulseResponseOptions &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<FiniteImpulseResponseOptionsImpl> (*options.pImpl);
    return *this;
}

/// Move assignment
FiniteImpulseResponseOptions&
FiniteImpulseResponseOptions::operator=(
    FiniteImpulseResponseOptions &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}


/// Destructor
FiniteImpulseResponseOptions::~FiniteImpulseResponseOptions() = default;

/// Beta
void FiniteImpulseResponseOptions::setBeta(const double beta)
{
    if (beta <= 0)
    {
        throw std::invalid_argument("beta " + std::to_string(beta)
                                  + " must be positive");
    }
    pImpl->mBeta = beta;
}

double FiniteImpulseResponseOptions::getBeta() const noexcept
{
    return pImpl->mBeta;
}

/// Order
void FiniteImpulseResponseOptions::setOrder(const uint16_t order) noexcept
{
    pImpl->mOrder = static_cast<int> (order);
}

int FiniteImpulseResponseOptions::getOrder() const noexcept
{
    return pImpl->mOrder;
}

///--------------------------------------------------------------------------///
///                             Filter                                       ///
///--------------------------------------------------------------------------///
template<class T>
FiniteImpulseResponse<T>::~FiniteImpulseResponse() = default;
