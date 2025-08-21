#include <iostream>
#include <algorithm>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uSignal/filterImplementations/infiniteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "src/alignment.hpp"

namespace UFR = USignal::FilterRepresentations;
using namespace USignal::FilterImplementations;

namespace
{

enum class Implementation
{
    DirectForm2Fast,
    DirectForm2Slow
};

template<typename T>
void iirDF2Transpose(const int order,
                     const USignal::Vector<T> bIn,
                     const USignal::Vector<T> aIn,
                     const int nSamples,
                     const USignal::Vector<T> xIn,
                     USignal::Vector<T> *yOut,
                     USignal::Vector<T> &vInputOutputDelay)
{
    const auto b = std::assume_aligned<ALIGNMENT> (bIn.data());
    const auto a = std::assume_aligned<ALIGNMENT> (aIn.data());
    const auto x = std::assume_aligned<ALIGNMENT> (xIn.data());
    auto y = std::assume_aligned<ALIGNMENT> (yOut->data());
    auto state
        = std::assume_aligned<ALIGNMENT> (vInputOutputDelay.data());
#ifndef NDEBUG
    assert(std::abs(a[0] - 1) < 10*std::numeric_limits<T>::epsilon());
#endif
    for (int i = 0; i < nSamples; ++i)
    {
        y[i] = state[0] + b[0]*x[i];
        for (int k = 1; k < order; ++k)
        {
            state[k - 1] = state[k] + b[k]*x[i] - a[k]*y[i];
        }
        state[order - 1] = b[order]*x[i] - a[order]*y[i];
    }
}


template<typename T>
void iirDF2TransposeWrong(const int order,
                     const USignal::Vector<T> bIn,
                     const USignal::Vector<T> aIn,
                     const int nSamples,
                     const USignal::Vector<T> xIn,
                     USignal::Vector<T> *yOut,
                     USignal::Vector<T> &vInputOutputDelay)
{
    if (nSamples == 0){return;}
    constexpr T zero{0};
    const auto b = std::assume_aligned<ALIGNMENT> (bIn.data());
    const auto a = std::assume_aligned<ALIGNMENT> (aIn.data());
    const auto x = std::assume_aligned<ALIGNMENT> (xIn.data());
    auto y = std::assume_aligned<ALIGNMENT> (yOut->data());
    auto vInputOutput
         = std::assume_aligned<ALIGNMENT> (vInputOutputDelay.data());
    for (int i = 0; i < nSamples; ++i)
    {
        // N.B. v[0] not used in sum loop
        T Xi{zero};
        T Yi{zero};
        for (int j = 1; j <= order; ++j)
        {
            Xi = Xi - a[j]*vInputOutput[j - 1];
            Yi = Yi + b[j]*vInputOutput[j - 1];
        }
        Xi   = Xi + x[i];
        T v0 = Xi;
        y[i] = Yi + b[0]*v0;
        std::shift_right(vInputOutput + 0, vInputOutput + order + 1, 1);
        vInputOutput[0] = v0;
    }
}

}

#ifdef WITH_IPP
#include <ipp/ipps.h>
#include <ipp/ippcore.h>
#include <ipp/ippvm.h>

namespace
{

::Implementation getImplementation(const int order)
{
return ::Implementation::DirectForm2Slow;
    if (order == 0){return ::Implementation::DirectForm2Slow;}
    if (order > 8)
    {
        Ipp64u featureMask;
        auto status = ippGetCpuFeatures(&featureMask, nullptr);
        if (status == ippStsNoErr)
        {
            auto enabledMask = ippGetEnabledCpuFeatures();
            if ((featureMask & ippCPUID_SSE42) &&
                (enabledMask & ippCPUID_SSE42))
            {
                return ::Implementation::DirectForm2Slow;
            }
            if ((featureMask & ippCPUID_AVX2) &&
                (enabledMask & ippCPUID_AVX2))
            {
                return ::Implementation::DirectForm2Slow;
            }
            if ((featureMask & ippCPUID_AVX512F) &&
                (enabledMask & ippCPUID_AVX512F))
            {
                return ::Implementation::DirectForm2Slow;
            }
        }
        else
        {
            return ::Implementation::DirectForm2Slow;
        }
    }
    return ::Implementation::DirectForm2Fast;
}
}

template<>
class InfiniteImpulseResponse<double>::InfiniteImpulseResponseImpl
{
public:
    explicit InfiniteImpulseResponseImpl(
        const UFR::InfiniteImpulseResponse<double> &filterCoefficients)
    {
        auto bs = filterCoefficients.getNumeratorFilterCoefficients();
        auto as = filterCoefficients.getDenominatorFilterCoefficients();        
        mOrder = filterCoefficients.getOrder();
        mImplementation = ::getImplementation(mOrder);
        // Normalize the filter coefficients
        auto a0 = as.at(0);
#ifndef NDEBUG
        assert(a0 != 0);
#endif
        bs = bs*(1./a0);
        as = as*(1./a0);
        as[0] = 1;
        mDelayInitialConditions.resize(bs.size(), 0);
        mDelayFinalConditions.resize(bs.size(), 0);
        if (mImplementation == ::Implementation::DirectForm2Fast)
        {
            auto status = ippsIIRGetStateSize_64f(mOrder, &mBufferSize);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error("Failed to get state size");
            } 
            mTaps.resize(2*(mOrder + 1), 0);
            mBuffer = ippsMalloc_8u(mBufferSize); 
            std::copy(bs.begin(), bs.end(), mTaps.begin() + 0);
            std::copy(as.begin(), as.end(), mTaps.begin() + mOrder + 1);
            mDelayFinalConditionsPtr
                = static_cast<Ipp64f *> (mDelayFinalConditions.data());
            mTapsPtr = static_cast<Ipp64f *> (mTaps.data());
            status = ippsIIRInit_64f(&mIIRState, mTapsPtr, mOrder,
                                     mDelayFinalConditionsPtr, mBuffer);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error("Failed to initialize IIR state");
            }
            mDelayInitialConditionsPtr
                = static_cast<Ipp64f *> (mDelayInitialConditions.data());
        }
        else
        {
            // Need these to be the same size for DF2Transpose
            mB.resize(mOrder + 1, 0);
            mA.resize(mOrder + 1, 0);
            std::copy(bs.begin(), bs.end(), mB.begin());
            std::copy(as.begin(), as.end(), mA.begin());
            mDelay.resize(mOrder + 1, 0);
        } 
        mInitialized = true;
    }
    ~InfiniteImpulseResponseImpl()
    {
        clear();
    }
    void clear() noexcept
    {
        mB.clear();
        mA.clear();
        if (mBuffer)
        {
            ippFree(mBuffer);
            mBuffer = nullptr;
        }
        mTapsPtr = nullptr;
        mDelayInitialConditionsPtr = nullptr;
        mDelayFinalConditionsPtr = nullptr;
        mTaps.clear();
        mDelayInitialConditions.clear();
        mDelayFinalConditions.clear();
        mIIRState = nullptr;
        mDelay.clear();
        mInitialConditions.clear();
    }
    void apply(const USignal::Vector<double> &x, USignal::Vector<double> *y)
    {
        auto nSamples = static_cast<int> (x.size());
        if (y->size() != nSamples){y->resize(nSamples, 0);}
        if (mImplementation == ::Implementation::DirectForm2Fast)
        {
            auto xPtr = static_cast<const Ipp64f *> (x.data());
            auto yPtr = static_cast<Ipp64f *> (y->data());
            auto status = ippsIIR_64f(xPtr, yPtr, nSamples, mIIRState);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error("Failed to apply filter");
            }
            if (!mIsPostProcessing)
            {
                //ippsIIRGetDlyLine_64f(mIIRState, pBufIPP64f_);
            }
        }
        else
        {
            iirDF2Transpose(mOrder,
                            mB, mA,
                            nSamples, x, y,
                            mDelay);
            if (mIsPostProcessing)
            {
/*
                std::copy(mInitialConditions.begin(),
                          mInitialConditions.end(),
                          mDelay.begin());
*/
            }
        }
    } 
//private:
    IppsIIRState_64f *mIIRState{nullptr};
    Ipp64f *mTapsPtr{nullptr}; // Filter taps (b and a) 
    Ipp64f *mDelayInitialConditionsPtr{nullptr};
    Ipp64f *mDelayFinalConditionsPtr{nullptr};
    Ipp8u *mBuffer{nullptr};
    USignal::Vector<double> mB;
    USignal::Vector<double> mA;
    USignal::Vector<double> mTaps;
    USignal::Vector<double> mDelayInitialConditions;
    USignal::Vector<double> mDelayFinalConditions;
    USignal::Vector<double> mDelay;
    USignal::Vector<double> mInitialConditions;
    ::Implementation mImplementation{::Implementation::DirectForm2Fast};
    int mBufferSize{0};
    int mOrder{0};
    bool mIsPostProcessing{false};
    bool mInitialized{false};
};

template<>
class InfiniteImpulseResponse<float>::InfiniteImpulseResponseImpl
{
public:
    explicit InfiniteImpulseResponseImpl(
        const UFR::InfiniteImpulseResponse<float> &filterCoefficients)
    {

    }
    bool mInitialized{false};
};
#else
static_assert(false, "Only IPP IIR filter implemented");
#endif

/// Constructor
template<class T>
InfiniteImpulseResponse<T>::InfiniteImpulseResponse(
    const UFR::InfiniteImpulseResponse<T> &filterCoefficients) :
    pImpl(std::make_unique<InfiniteImpulseResponseImpl> (filterCoefficients))
{

}

/// Destructor
template<class T>
InfiniteImpulseResponse<T>::~InfiniteImpulseResponse() = default;

/// Initialized?
template<class T>
bool InfiniteImpulseResponse<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Apply
template<class T>
void InfiniteImpulseResponse<T>::apply()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    const auto x = this->getInputReference();
    if (x.empty())
    {   
        USignal::Vector<T> y;
        this->setOutput(std::move(y));
    }
    USignal::Vector<T> y;
    pImpl->apply(x, &y);
    this->setOutput(std::move(y));
}

///--------------------------------------------------------------------------///
template class USignal::FilterImplementations::InfiniteImpulseResponse<double>;
//template class USignal::FilterImplementations::InfiniteImpulseResponse<float>;

