#include <iostream>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uSignal/filterImplementations/secondOrderSections.hpp"
#include "uSignal/filterRepresentations/secondOrderSections.hpp"
#include "uSignal/vector.hpp"
#include "src/alignment.hpp"

using namespace USignal::FilterImplementations;
namespace UFR = USignal::FilterRepresentations;

#ifdef WITH_IPP
#include <ipp/ipps.h>
template<>
class SecondOrderSections<double>::SecondOrderSectionsImpl
{
public:
    SecondOrderSectionsImpl(
        const UFR::SecondOrderSections<double> &filterCoefficients)
    {
        mSections = filterCoefficients.getNumberOfSections();
        if (mSections < 1)
        {
            throw std::invalid_argument("No sections");
        }
        auto bs = filterCoefficients.getNumeratorFilterCoefficients();
        auto as = filterCoefficients.getDenominatorFilterCoefficients();
#ifndef NDEBUG
        assert(bs.size() == 3*mSections);
        assert(as.size() == 3*mSections);
#endif
        // Buffer
        int bufferSize{0};
        auto status = ippsIIRGetStateSize_BiQuad_64f(mSections,
                                                     &bufferSize);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error("Failed to get state size");
        }
        pBuffer = ippsMalloc_8u(bufferSize);
        // Taps
        pTaps = ippsMalloc_64f(2*3*mSections);
        for (int i = 0; i < mSections; ++i)
        {
            pTaps[6*i + 0] = bs.at(3*i + 0);
            pTaps[6*i + 1] = bs.at(3*i + 1);
            pTaps[6*i + 2] = bs.at(3*i + 2);
            pTaps[6*i + 3] = as.at(3*i + 0);
            pTaps[6*i + 4] = as.at(3*i + 1);
            pTaps[6*i + 5] = as.at(3*i + 2);
        }
        // Source/delay lines
        mInitialConditions.resize(2*mSections, 0);
        auto workSpace = std::max(8, 2*mSections);
        pDelaySource = ippsMalloc_64f(workSpace);
        pDelayDest = ippsMalloc_64f(workSpace);
        constexpr double zero{0};
        std::fill(pDelaySource, pDelaySource + workSpace, zero);
        std::fill(pDelayDest,   pDelayDest + workSpace,   zero);
        // State
        status = ippsIIRInit_BiQuad_64f(&pState,
                                        pTaps,
                                        mSections,
                                        pDelaySource,
                                        pBuffer);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error(
                "Failed to initialize IPP SOS filter state");
            releaseMemory();
        }
        mInitialized = true;
    }
    /// Reset the initial conditions
    void resetInitialConditions() noexcept
    {
        std::copy(mInitialConditions.begin(), mInitialConditions.end(),
                  pDelaySource);
    }
    /// Destructor
    ~SecondOrderSectionsImpl()
    {
        releaseMemory();
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
        pState = nullptr;    
        mInitialized = false;
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
        auto status = ippsIIRSetDlyLine_64f(pState, pDelaySource);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error("Failed to set delay line");
        }
        status = ippsIIR_64f(xPtr, yPtr, n, pState);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error("Failed to apply SOS filter");
        }
        if (mIsRealTime)
        {
            status = ippsIIRGetDlyLine_64f(pState, pDelayDest);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error(
                    "Failed to get filter final conditions");
            }
            std::copy(pDelayDest, pDelayDest + 2*mSections, pDelaySource);
        }
    }
//private:
    USignal::Vector<double> mInitialConditions;
    IppsIIRState_64f *pState{nullptr};
    Ipp64f *pTaps{nullptr};
    Ipp64f *pDelaySource{nullptr};
    Ipp64f *pDelayDest{nullptr};
    Ipp8u *pBuffer{nullptr};
    int mSections{0};
    bool mIsRealTime{false};
    bool mInitialized{false};
};

///--------------------------------------------------------------------------///
///                                   Float                                  ///
///--------------------------------------------------------------------------///

template<>
class SecondOrderSections<float>::SecondOrderSectionsImpl
{
public:
    SecondOrderSectionsImpl(
        const UFR::SecondOrderSections<float> &filterCoefficients)
    {
        mSections = filterCoefficients.getNumberOfSections();
        if (mSections < 1)
        {
            throw std::invalid_argument("No sections");
        }
        auto bs = filterCoefficients.getNumeratorFilterCoefficients();
        auto as = filterCoefficients.getDenominatorFilterCoefficients();
#ifndef NDEBUG
        assert(bs.size() == 3*mSections);
        assert(as.size() == 3*mSections);
#endif
        // Buffer
        int bufferSize{0};
        auto status = ippsIIRGetStateSize_BiQuad_32f(mSections,
                                                     &bufferSize);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error("Failed to get state size");
        }
        pBuffer = ippsMalloc_8u(bufferSize);
        // Taps
        pTaps = ippsMalloc_32f(2*3*mSections);
        for (int i = 0; i < mSections; ++i)
        {
            pTaps[6*i + 0] = bs.at(3*i + 0);
            pTaps[6*i + 1] = bs.at(3*i + 1);
            pTaps[6*i + 2] = bs.at(3*i + 2);
            pTaps[6*i + 3] = as.at(3*i + 0);
            pTaps[6*i + 4] = as.at(3*i + 1);
            pTaps[6*i + 5] = as.at(3*i + 2);
        }
        // Source/delay lines
        mInitialConditions.resize(2*mSections, 0);
        auto workSpace = std::max(8, 2*mSections);
        pDelaySource = ippsMalloc_32f(workSpace);
        pDelayDest = ippsMalloc_32f(workSpace);
        constexpr float zero{0};
        std::fill(pDelaySource, pDelaySource + workSpace, zero);
        std::fill(pDelayDest,   pDelayDest + workSpace,   zero);
        // State
        status = ippsIIRInit_BiQuad_32f(&pState,
                                        pTaps,
                                        mSections,
                                        pDelaySource,
                                        pBuffer);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error(
                "Failed to initialize IPP SOS filter state");
            releaseMemory();
        }
        mInitialized = true;
    }
    /// Destructor
    ~SecondOrderSectionsImpl()
    {   
        releaseMemory();
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
        pState = nullptr;    
        mInitialized = false;
    }   
    /// Apply the filter
    void apply(const Vector<float> &x, Vector<float> *y)
    {
        if (y == nullptr){throw std::runtime_error("y is null");}
        if (x.empty()){return;}
        if (y->size() != x.size()){y->resize(x.size(), 0);}
        auto n = static_cast<int> (x.size());
        const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
        auto yPtr = std::assume_aligned<ALIGNMENT> (y->data());
        auto status = ippsIIRSetDlyLine_32f(pState, pDelaySource);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error("Failed to set delay line");
        }
        status = ippsIIR_32f(xPtr, yPtr, n, pState);
        if (status != ippStsNoErr)
        {
            throw std::runtime_error("Failed to apply SOS filter");
        }
        if (mIsRealTime)
        {
            status = ippsIIRGetDlyLine_32f(pState, pDelayDest);
            if (status != ippStsNoErr)
            {
                throw std::runtime_error(
                    "Failed to get filter final conditions");
            }
            std::copy(pDelayDest, pDelayDest + 2*mSections, pDelaySource);
        }
    }
//private:
    USignal::Vector<float> mInitialConditions;
    IppsIIRState_32f *pState{nullptr};
    Ipp32f *pTaps{nullptr};         
    Ipp32f *pDelaySource{nullptr};
    Ipp32f *pDelayDest{nullptr};
    Ipp8u *pBuffer{nullptr};
    int mSections{0};
    bool mIsRealTime{false};
    bool mInitialized{false};
};

#else
static_assert(false, "Only IPP SOS filter implemented");
#endif

/// Constructor
template<class T>
SecondOrderSections<T>::SecondOrderSections(
    const UFR::SecondOrderSections<T> &filterCoefficients) :
    pImpl(std::make_unique<SecondOrderSectionsImpl> (filterCoefficients))
{
}

/// Initialized?
template<class T>
bool SecondOrderSections<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Apply
template<class T>
void SecondOrderSections<T>::apply()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
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
template<class T>
SecondOrderSections<T>::~SecondOrderSections() = default;

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class USignal::FilterImplementations::SecondOrderSections<double>;
template class USignal::FilterImplementations::SecondOrderSections<float>;

