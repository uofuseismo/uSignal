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
template<typename T>
class FiniteImpulseResponse<T>::FiniteImpulseResponseImpl
{
public:
    FiniteImpulseResponseImpl(
        const UFR::FiniteImpulseResponse<T> &filterCoefficients,
        Implementation implementation,
        const bool isRealTime) :
        mIsRealTime(isRealTime)
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
        // FIR implementation
        IppAlgType algorithmType = ippAlgDirect;
        if (implementation == Implementation::FFT)
        {
            algorithmType = ippAlgFFT;
        }
        if (implementation == Implementation::Automatic)
        {
            algorithmType = ippAlgAuto;
        }

        // State
        int stateSize{0};
        int bufferSize{0};
        if constexpr (std::is_same<T, double>::value)
        {
            auto status = ippsFIRSRGetSize(mOrder + 1,
                                           ipp64f,
                                           &stateSize,
                                           &bufferSize);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error("Failed to get FIR state size");
            }
            pState64f
                = reinterpret_cast<IppsFIRSpec_64f *>
                  (ippsMalloc_8u(stateSize));
            pBuffer = ippsMalloc_8u(bufferSize);

            // Taps
            pTaps64f = ippsMalloc_64f(bs.size());
            std::copy(bs.begin(), bs.end(), pTaps64f);

            // Delay lines
            auto workSpace = std::max(mOrder + 1, 8);
            pDelaySource64f = ippsMalloc_64f(workSpace);
            pDelayDest64f = ippsMalloc_64f(workSpace);
            constexpr double zero{0};
            std::fill(pDelaySource64f, pDelaySource64f + workSpace, zero);
            std::fill(pDelayDest64f ,  pDelayDest64f + workSpace,   zero);
            if (mOrder > 0){mInitialConditions.resize(mOrder, zero);}

            // Initialize FIR filter
            status = ippsFIRSRInit_64f(pTaps64f,
                                       static_cast<int> (bs.size()),
                                       algorithmType,
                                       pState64f);
            if (status != ippStsNoErr)
            {
                releaseMemory();
                throw std::runtime_error(
                    "Failed to initialize double FIR state");
            }
        }
        else if constexpr (std::is_same<T, float>::value)
        {
            auto status = ippsFIRSRGetSize(mOrder + 1,
                                           ipp32f,
                                           &stateSize,
                                           &bufferSize);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error("Failed to get FIR state size");
            }
            pState32f
                = reinterpret_cast<IppsFIRSpec_32f *>
                  (ippsMalloc_8u(stateSize));
            pBuffer = ippsMalloc_8u(bufferSize);

            // Taps
            pTaps32f = ippsMalloc_32f(bs.size());
            std::copy(bs.begin(), bs.end(), pTaps32f);

            // Delay lines
            auto workSpace = std::max(mOrder + 1, 8); 
            pDelaySource32f = ippsMalloc_32f(workSpace);
            pDelayDest32f = ippsMalloc_32f(workSpace);
            constexpr float zero{0};
            std::fill(pDelaySource32f, pDelaySource32f + workSpace, zero);
            std::fill(pDelayDest32f ,  pDelayDest32f + workSpace,   zero);
            if (mOrder > 0){mInitialConditions.resize(mOrder, zero);}

            // Initialize FIR filter
            status = ippsFIRSRInit_32f(pTaps32f,
                                       static_cast<int> (bs.size()),
                                       algorithmType,
                                       pState32f);
            if (status != ippStsNoErr)
            {
                releaseMemory();
                throw std::runtime_error(
                    "Failed to initialize float FIR state");
            }
        }
        else
        {
            static_assert(false, "Unhandled precision");
        }
        mInitialized = true;
    }
    /// Sets the initial conditions
    void setInitialConditions(const USignal::Vector<T> &initialConditions)
    {
        if (static_cast<int> (initialConditions.size()) != mOrder)
        {
            throw std::invalid_argument("Initial conditions length "
                                      + std::to_string(initialConditions.size())
                                      + " must equal filter order " 
                                      + std::to_string(mOrder));
        } 
        if (mOrder > 0)
        {
            std::copy(initialConditions.begin(),
                      initialConditions.end(),
                      mInitialConditions.begin());
        }
    }
    /// Release memory
    void releaseMemory() noexcept
    {
        mInitialConditions.clear();
        if (pBuffer)
        {
            ippsFree(pBuffer);
            pBuffer = nullptr;
        }
        if (pTaps64f)
        {
            ippsFree(pTaps64f);
            pTaps64f = nullptr;
        }
        if (pTaps32f)
        {
            ippsFree(pTaps32f);
            pTaps32f = nullptr;
        }
        if (pDelaySource64f)
        {
            ippsFree(pDelaySource64f);
            pDelaySource64f = nullptr;
        }
        if (pDelaySource32f)
        {
            ippsFree(pDelaySource32f);
            pDelaySource32f = nullptr;
        }
        if (pDelayDest64f)
        {
            ippsFree(pDelayDest64f);
            pDelayDest64f = nullptr;
        }
        if (pDelayDest32f)
        {
            ippsFree(pDelayDest32f);
            pDelayDest32f = nullptr;
        }
        if (pState64f)
        {
            ippsFree(pState64f);
            pState64f = nullptr;
        }
        if (pState32f)
        {
            ippsFree(pState32f);
            pState32f = nullptr;
        }
        mInitialized = false;
    }
    /// Reset initial conditions
    void resetInitialConditions() noexcept
    {
        if (mOrder > 0)
        {
            if constexpr (std::is_same<T, double>::value)
            {
                std::copy(mInitialConditions.begin(), mInitialConditions.end(),
                          pDelaySource64f);
            }
            else if constexpr (std::is_same<T, float>::value)
            {
                std::copy(mInitialConditions.begin(), mInitialConditions.end(),
                          pDelaySource32f);
            }
            else
            {
                static_assert(false, "Unhandled precision");
            }
        }
    }
    /// Apply the filter
    void apply(const USignal::Vector<T> &x,
               USignal::Vector<T> *y)
    {
        if (y == nullptr){throw std::runtime_error("y is null");}
        if (x.empty()){return;}
        if (y->size() != x.size()){y->resize(x.size(), 0);}
        auto n = static_cast<int> (x.size());
        if constexpr (std::is_same<T, double>::value)
        {
            const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
            auto yPtr = std::assume_aligned<ALIGNMENT> (y->data());
            auto status = ippsFIRSR_64f(xPtr, yPtr, n,
                                        pState64f,
                                        pDelaySource64f, pDelayDest64f,
                                        pBuffer);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error(
                    "Failed to apply double FIR single rate filter");
            }
            if (mIsRealTime && mOrder > 0)
            {
                std::copy(pDelayDest64f, pDelayDest64f + mOrder,
                          pDelaySource64f);
            }
        }
        else if constexpr (std::is_same<T, float>::value)
        {
            const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
            auto yPtr = std::assume_aligned<ALIGNMENT> (y->data());
            auto status = ippsFIRSR_32f(xPtr, yPtr, n,
                                        pState32f,
                                        pDelaySource32f, pDelayDest32f, 
                                        pBuffer);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error(
                    "Failed to apply double FIR single rate filter");
            }
            if (mIsRealTime && mOrder > 0)
            {
                std::copy(pDelayDest32f, pDelayDest32f + mOrder, 
                          pDelaySource32f);
            }
        }
        else
        {
            static_assert(false, "Unhandled precision");
        }
    }

//private:
    USignal::Vector<T> mInitialConditions;
    IppsFIRSpec_64f *pState64f{nullptr};
    Ipp64f *pTaps64f{nullptr};
    Ipp64f *pDelaySource64f{nullptr};
    Ipp64f *pDelayDest64f{nullptr};
    IppsFIRSpec_32f *pState32f{nullptr};
    Ipp32f *pTaps32f{nullptr};
    Ipp32f *pDelaySource32f{nullptr};
    Ipp32f *pDelayDest32f{nullptr};
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
    const Implementation implementation,
    const bool isRealTime) :
    pImpl(std::make_unique<FiniteImpulseResponseImpl>
             (filterCoefficients, implementation, isRealTime)
          )
{
}

/// Initialized?
template<class T>
bool FiniteImpulseResponse<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/*
/// Initial conditions
template<class T>
void FiniteImpulseResponse<T>::setInitialConditions(
    const USignal::Vector<T> &initialConditions)
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class not initialized");
    }
    pImpl->setInitialConditions(initialConditions);
}
*/

/// Apply
template<class T>
void FiniteImpulseResponse<T>::apply()
{
    if (!isInitialized())
    {
        throw std::runtime_error("FIR filter not initialized");
    }
    const auto x = this->getInputReference();
    USignal::Vector<T> y;
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
template class USignal::FilterImplementations::FiniteImpulseResponse<float>;

