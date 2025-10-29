#define NDEBUG 1
#include <iostream>
#include <string>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uSignal/transforms/hilbert/finiteImpulseResponse.hpp"
#include "uSignal/filterImplementations/finiteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/filterDesign/finiteImpulseResponse/hilbertTransformer.hpp"
#include "src/alignment.hpp"

using namespace USignal::Transforms::Hilbert;

template<class T>
class FiniteImpulseResponse<T>::FiniteImpulseResponseImpl
{
public:
    FiniteImpulseResponseImpl(const FiniteImpulseResponseOptions &options,
                              const bool isRealTime) :
        mIsRealTime(isRealTime)
    {
        namespace UFD = USignal::FilterDesign::FiniteImpulseResponse;
        namespace UFR = USignal::FilterRepresentations;
        namespace UFI = USignal::FilterImplementations;
        constexpr auto implementation
        {
            UFI::FiniteImpulseResponse<T>::Implementation::Direct
        };
        auto filter
            = UFD::hilbertTransformer<double>
              (options.getOrder(), options.getBeta());
        mIsTypeIII = filter.getOrder()%2 == 0 ? true : false;
        if (mIsTypeIII)
        {
            auto filterTaps = filter.getFilterCoefficients();
            auto nFilterTaps = static_cast<int> (filterTaps.size());
            mTypeIIIShift = static_cast<int> ((nFilterTaps - 1)/2);
            USignal::Vector<T> realTaps(filterTaps.size()); // TODO 
            USignal::Vector<T> imaginaryTaps(filterTaps.size());
            for (int i = 0; i < nFilterTaps; ++i)
            {
                realTaps[i] = std::real(filterTaps[i]);
                imaginaryTaps[i] = std::imag(filterTaps[i]);
#ifndef NDEBUG
                if (i == (nFilterTaps - 1)/2)
                {
                    assert(std::abs(std::real(filterTaps[i] - 1) <
                           std::numeric_limits<T>::epsilon()*100);
                }
                else
                {
                    assert(std::abs(std::real(filterTaps[i]) <
                           std::numeric_limits<T>::epsilon()*100);
                }
#endif
            }
            // TODO suss out delay lines and just implement this signal
            UFR::FiniteImpulseResponse<T> realFilter{realTaps};
            mRealFIRFilter
                = std::make_unique<UFI::FiniteImpulseResponse<T>>
                  (realFilter, implementation, mIsRealTime);

            UFR::FiniteImpulseResponse<T> imaginaryFilter{imaginaryTaps};
            mImaginaryFIRFilter
                = std::make_unique<UFI::FiniteImpulseResponse<T>>
                  (imaginaryFilter, implementation, mIsRealTime);
/*
            if (nFilterTaps > 0)
            {
                constexpr T zero{0};
                mInitialConditionsTypeIII.resize(nFilterTaps - 1, zero);
                mFinalConditionsTypeIII.resize(nFilterTaps - 1, zero);
            }
*/
        }
        else
        {
            auto filterTaps = filter.getFilterCoefficients();
            USignal::Vector<T> realTaps(filterTaps.size());
            USignal::Vector<T> imaginaryTaps(filterTaps.size());
            for (int i = 0; i < static_cast<int> (filterTaps.size()); ++i)
            {
                realTaps[i] = std::real(filterTaps[i]);
                imaginaryTaps[i] = std::imag(filterTaps[i]);
            }
            UFR::FiniteImpulseResponse<T> realFilter{realTaps};
            UFR::FiniteImpulseResponse<T> imaginaryFilter{imaginaryTaps};
            mRealFIRFilter
                = std::make_unique<UFI::FiniteImpulseResponse<T>>
                  (realFilter, implementation, mIsRealTime);
            mImaginaryFIRFilter
                = std::make_unique<UFI::FiniteImpulseResponse<T>>
                  (imaginaryFilter, implementation, mIsRealTime);
        }
        mOrder = options.getOrder();
        mInitialized = true;
    }
    void apply(const USignal::Vector<T> &x, USignal::Vector<std::complex<T>> *y)
    {
#ifndef NDEBUG
        assert(y != nullptr);
#endif
        if (x.empty())
        {
            y->resize(0);
            return;
        }
        if (mIsTypeIII)
        {
            // TODO suss out delay lines for 
            mRealFIRFilter->setInput(x);
            mRealFIRFilter->apply();

            mImaginaryFIRFilter->setInput(x);
            mImaginaryFIRFilter->apply();

            const auto yRealRef
                = mRealFIRFilter->getOutputReference();
            const auto yImaginaryRef
                = mImaginaryFIRFilter->getOutputReference();
            if (y->size() != yImaginaryRef.size())
            {
                constexpr std::complex<T> zero{0 + 0i};
                y->resize(yImaginaryRef.size(), zero);
            }
            auto yPtr = std::assume_aligned<ALIGNMENT> (y->data()); 
            //const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
            const auto yRealPtr
                = std::assume_aligned<ALIGNMENT> (yRealRef.data());
            const auto yImaginaryPtr
                = std::assume_aligned<ALIGNMENT> (yImaginaryRef.data());
            for (int i = 0; i < static_cast<int> (yImaginaryRef.size()); ++i)
            {
                yPtr[i] = std::complex<T> (yRealPtr[i], yImaginaryPtr[i]);
            }
            if (!mIsRealTime)
            {
                mRealFIRFilter->resetInitialConditions();
                mImaginaryFIRFilter->resetInitialConditions();
            }
        }
        else
        {
            mRealFIRFilter->setInput(x);
            mRealFIRFilter->apply();
 
            mImaginaryFIRFilter->setInput(x);
            mImaginaryFIRFilter->apply();

            const auto yRealRef
                = mRealFIRFilter->getOutputReference();
            const auto yImaginaryRef
                = mImaginaryFIRFilter->getOutputReference();
#ifndef NDEBUG
            assert(yRealRef.size() == yImaginaryRef.size());
#endif
            if (y->size() != yRealRef.size())
            {
                constexpr std::complex<T> zero{0 + 0i};
                y->resize(yRealRef.size(), zero);
            }
            auto yPtr = std::assume_aligned<ALIGNMENT> (y->data());
            const auto yRealPtr
                = std::assume_aligned<ALIGNMENT> (yRealRef.data());
            const auto yImaginaryPtr
                = std::assume_aligned<ALIGNMENT> (yImaginaryRef.data());
            for (int i = 0; i < static_cast<int> (yRealRef.size()); ++i)
            {
                yPtr[i] = std::complex<T> (yRealPtr[i], yImaginaryPtr[i]);
            }
            if (!mIsRealTime)
            {
                mRealFIRFilter->resetInitialConditions();
                mImaginaryFIRFilter->resetInitialConditions();
            }
        }
    }
    void resetInitialConditions()
    {
        mRealFIRFilter->resetInitialConditions();
        mImaginaryFIRFilter->resetInitialConditions();
    } 
    void setInitialConditions(const USignal::Vector<T> &initialConditions)
    {
        mRealFIRFilter->setInitialConditions(initialConditions);
        mImaginaryFIRFilter->setInitialConditions(initialConditions);
    }
//private:
    std::unique_ptr<USignal::FilterImplementations::FiniteImpulseResponse<T>>
        mRealFIRFilter{nullptr};
    std::unique_ptr<USignal::FilterImplementations::FiniteImpulseResponse<T>>
        mImaginaryFIRFilter{nullptr};
    USignal::Vector<T> mInitialConditionsTypeIII;
    USignal::Vector<T> mFinalConditionsTypeIII;
    int64_t mTypeIIIShift{0};
    int mOrder{0};
    bool mIsRealTime{false};
    bool mIsTypeIII{false};
    bool mInitialized{false};
};

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
FiniteImpulseResponse<T>::FiniteImpulseResponse(
    const FiniteImpulseResponseOptions &options,
    const bool isRealTime) :
    pImpl(std::make_unique<FiniteImpulseResponseImpl> (options, isRealTime))
{
}

template<class T>
bool FiniteImpulseResponse<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

template<class T>
void FiniteImpulseResponse<T>::apply()
{
    if (!isInitialized())
    {
        throw std::invalid_argument("FIR Hilbert transform not initialized");
    }   
    const auto &xRef = this->getInputReference();
    if (xRef.empty()){throw std::invalid_argument("Input signal is empty");}
    // Transform
    USignal::Vector<std::complex<T>> y;
    pImpl->apply(xRef, &y);
    this->setOutput(std::move(y)); 
}

template<class T>
void FiniteImpulseResponse<T>::setInitialConditions(
    const USignal::Vector<T> &initialConditions)
{
    if (!isInitialized())
    {
        throw std::invalid_argument("FIR Hilbert transform not initialized");
    }
    auto order = pImpl->mOrder;
    if (static_cast<int> (initialConditions.size()) != order)
    {
        throw std::invalid_argument("Initial conditions length = " 
                                  + std::to_string(initialConditions.size())
                                  + " should equal " 
                                  + std::to_string(order));
    }
    pImpl->setInitialConditions(initialConditions);
}


template<class T>
void FiniteImpulseResponse<T>::resetInitialConditions()
{
    if (!isInitialized())
    {
        throw std::invalid_argument("FIR Hilbert transform not initialized");
    }
    pImpl->resetInitialConditions();
}


template<class T>
FiniteImpulseResponse<T>::~FiniteImpulseResponse() = default;

///--------------------------------------------------------------------------///
///                           Template Instantiation                         ///
///--------------------------------------------------------------------------///
template class USignal::Transforms::Hilbert::FiniteImpulseResponse<double>;
template class USignal::Transforms::Hilbert::FiniteImpulseResponse<float>;

