#include <iostream>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uSignal/filterImplementations/finiteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include "src/alignment.hpp"

using namespace USignal::FilterImplementations;
namespace UFR = USignal::FilterRepresentations;

#ifdef WITH_IPP
#include <ipp/ipps.h>
template<>
class FiniteImpulseResponse<double>::FiniteImpulseResponseImpl
{
public:
    FiniteImpulseResponseImpl(
        const UFR::FiniteImpulseResponse<double> &filterCoefficients,
        Implementation implementation)
    {
        mOrder = filterCoefficients.getOrder();
        if (mOrder < 0)
        {
            throw std::invalid_argument("Order must be at least zero");
        }
        auto bs = filterCoefficients.getFilterCoefficients();
#ifndef NDEBUG
        assert(!bs.empty());
        assert(static_cast<int> (bs.size()) == mOrder + 1);
#endif
        // State
        int stateSize{0};
        int bufferSize{0};
        auto status = ippsFIRSRGetSize(mOrder + 1,
                                       ipp64f,
                                       &stateSize,
                                       &bufferSize);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error("Failed to get FIR state size");
        }
        pState = reinterpret_cast<IppsFIRSpec_64f *> (ippsMalloc_8u(stateSize));
        pBuffer = ippsMalloc_8u(bufferSize);

        // Taps
        pTaps = ippsMalloc_64f(bs.size());
        std::copy(bs.begin(), bs.end(), pTaps);

        // Delay lines
        auto workSpace = std::max(mOrder + 1, 8);
        pDelaySource = ippsMalloc_64f(workSpace);
        pDelayDest = ippsMalloc_64f(workSpace);
        constexpr double zero{0};
        std::fill(pDelaySource, pDelaySource + workSpace, zero);
        std::fill(pDelayDest ,  pDelayDest + workSpace,   zero);
        if (mOrder > 0){mInitialConditions.resize(mOrder, zero);}

        // Initialize FIR filter
        IppAlgType algorithmType = ippAlgDirect;
        if (implementation == Implementation::FFT)
        {
            algorithmType = ippAlgFFT;
        }
        if (implementation == Implementation::Automatic)
        {
            algorithmType = ippAlgAuto;
        }
        status = ippsFIRSRInit_64f(pTaps,
                                   static_cast<int> (bs.size()),
                                   algorithmType,
                                   pState);
        if (status != ippStsNoErr)
        {
            releaseMemory();
            throw std::runtime_error("Failed to initialize FIR state");
        }
        mInitialized = true;
    }
    /// Release memory
    void releaseMemory() noexcept
    {
        mInitialConditions.clear();
        if (pTaps)
        {
            ippsFree(pTaps);
            pTaps = nullptr;
        }
        if (pBuffer)
        {
            ippsFree(pBuffer);
            pBuffer = nullptr;
        }
        if (pDelaySource)
        {
            ippsFree(pDelaySource);
            pDelaySource = nullptr;
        }
        if (pDelayDest)
        {
            ippsFree(pDelayDest);
            pDelayDest = nullptr;
        }
        if (pState)
        {
            ippsFree(pState);
            pState = nullptr;
        }
        mInitialized = false;
    }
    /// Reset initial conditions
    void resetInitialConditions() noexcept
    {
        if (mOrder > 0)
        {
            std::copy(mInitialConditions.begin(), mInitialConditions.end(),
                      pDelaySource);
        }
    }
    /// Apply the filter
    void apply(const Vector<double> &x, Vector<double> *y)
    {
        if (y == nullptr){throw std::runtime_error("y is null");}
        if (x.empty()){return;}
        if (y->size() != x.size()){y->resize(x.size(), 0);}
        auto n = static_cast<int> (x.size());
        const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
        auto yPtr = std::assume_aligned<ALIGNMENT> (y->data());
        auto status = ippsFIRSR_64f(xPtr, yPtr, n,
                                    pState,
                                    pDelaySource, pDelayDest, pBuffer);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error("Failed to apply FIR single rate filter");
        }
        if (mIsRealTime && mOrder > 0)
        {
            std::copy(pDelayDest, pDelayDest + mOrder, pDelaySource);
        }
    }

//private:
    Vector<double> mInitialConditions;
    IppsFIRSpec_64f *pState{nullptr};
    Ipp64f *pTaps{nullptr};
    Ipp64f *pDelaySource{nullptr};
    Ipp64f *pDelayDest{nullptr};
    Ipp8u *pBuffer{nullptr};
    int mOrder{0}; 
    bool mIsRealTime{false};
    bool mInitialized{false};
};

#else
static_assert(false, "Only IPP FIR filter implemented");
#endif

/// Constructor
template<class T>
FiniteImpulseResponse<T>::FiniteImpulseResponse(
    const UFR::FiniteImpulseResponse<T> &filterCoefficients,
    const Implementation implementation) :
    pImpl(std::make_unique<FiniteImpulseResponseImpl>
             (filterCoefficients, implementation)
          )
{
}

/// Initialized?
template<class T>
bool FiniteImpulseResponse<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Apply
template<class T>
void FiniteImpulseResponse<T>::apply()
{
    if (!isInitialized())
    {
        throw std::runtime_error("FIR filter not initialized");
    }
    const auto x = this->getInputReference();
    Vector<T> y;
    if (x.empty())
    {
        y.resize(0);
        this->setOutput(std::move(y));
    }
    y.resize(x.size());
    pImpl->apply(x, &y);
    this->setOutput(std::move(y));
}

/// Destructor
template<typename T>
FiniteImpulseResponse<T>::~FiniteImpulseResponse() = default;

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class USignal::FilterImplementations::FiniteImpulseResponse<double>;
//template class USignal::FilterImplementations::FiniteImpulseResponse<float>;

