#include <iostream>
#include <string>
#include <cmath>
#include <type_traits>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uSignal/transforms/fourier/forward.hpp"
#include "uSignal/vector.hpp"

#ifdef WITH_IPP
#include <ipp/ipps.h>
#else
static_assert(false, "Only IPP DFT implemented");
#endif

using namespace USignal::Transforms::Fourier;

namespace
{

///@param[out] specSize        FFT specification structure size value
///@param[out] specBufferSize  Buffer size value for FFT initialization function
///@param[out] bufferSize      Size of the FFT external work buffer
///@param[in] order            FFT order.  The input signal length N = 2^{order}
template<typename T>
void getFFTSize(int *specSize,
                   int *specBufferSize, 
                   int *bufferSize,
                   const int order,
                   const IppHintAlgorithm hint = ippAlgHintNone)
{
    if constexpr (std::is_same<T, double>::value)
    { 
        auto status = ippsFFTGetSize_R_64f(order,
                                           IPP_FFT_DIV_INV_BY_N,
                                           hint,
                                           specSize,
                                           specBufferSize,
                                           bufferSize);
        if (status != ippStsNoErr)
        {
            *specSize = -1;
            *specBufferSize =-1;
            *bufferSize =-1;
            throw std::runtime_error("Failed to get double fft sizes");
        }
    }
    else if constexpr (std::is_same<T, float>::value)
    {
        auto status = ippsFFTGetSize_R_32f(order,
                                           IPP_FFT_DIV_INV_BY_N,
                                           hint,
                                           specSize,
                                           specBufferSize,
                                           bufferSize);
        if (status != ippStsNoErr)
        {
            *specSize = -1; 
            *specBufferSize =-1;
            *bufferSize =-1;
            throw std::runtime_error("Failed to get float fft sizes");
        } 
    } 
    else
    {
        throw std::runtime_error("Unhandled template value in fft size");
    }
}

template<typename T>
void getDFTSize(int *specSize,
                int *specBufferSize,
                int *bufferSize,
                const int signalLength,
                const IppHintAlgorithm hint = ippAlgHintNone)
{
    if constexpr (std::is_same<T, double>::value)
    { 
        auto status = ippsDFTGetSize_R_64f(signalLength,
                                           IPP_FFT_DIV_INV_BY_N,
                                           hint,
                                           specSize,
                                           specBufferSize,
                                           bufferSize);
        if (status != ippStsNoErr)
        {
            *specSize = -1; 
            *specBufferSize =-1;
            *bufferSize =-1;
            throw std::runtime_error("Failed to get double dft sizes");
        }   
    }
    else if constexpr (std::is_same<T, float>::value)
    {                                      
        auto status = ippsDFTGetSize_R_32f(signalLength,
                                           IPP_FFT_DIV_INV_BY_N,
                                           hint,
                                           specSize,
                                           specBufferSize,
                                           bufferSize);
        if (status != ippStsNoErr)
        {
            *specSize = -1;
            *specBufferSize =-1;
            *bufferSize =-1;
            throw std::runtime_error("Failed to get float dft sizes");
        }
    }
    else
    {
        throw std::runtime_error("Unhandled template value in dft size");
    }

}

}

template<class T>
class Forward<T>::ForwardImpl
{
public:
    explicit ForwardImpl(const ForwardParameters &parameters) :
        mParameters(parameters)
    {
        mInitialized = true;
    }
    ~ForwardImpl()
    {
        clear();
    }
    void clear()
    {
        if (mFFTSpec64f){ippsFree(mFFTSpec64f);}
        if (mDFTSpec64f){ippsFree(mDFTSpec64f);}
        if (mFFTSpec32f){ippsFree(mFFTSpec32f);}
        if (mDFTSpec32f){ippsFree(mDFTSpec32f);}
        if (mBuffer){ippsFree(mBuffer);}
        mFFTSpec64f = nullptr;
        mDFTSpec64f = nullptr;
        mFFTSpec32f = nullptr;
        mDFTSpec32f = nullptr;
        mBuffer = nullptr;
        mDoFFT = false;
        mAllocatedBuffers = false;
    }
    void initialize(const int signalLength)
    {
        auto doFFT = mParameters.getImplementation()
                     == ForwardParameters::Implementation::Fast ?
                        true : false;
        int order{-1};
        auto orderWork
            = static_cast<int>
              (std::floor(std::log2(static_cast<double> (signalLength))));
        auto signalLengthNextPowerOf2
            = static_cast<int> (std::pow(2, orderWork));
        if (signalLengthNextPowerOf2 < signalLength)
        {
            signalLengthNextPowerOf2 = signalLengthNextPowerOf2*2;
#ifndef NDEBUG
            assert(signalLengthNextPowerOf2 >= signalLength);
#endif 
        }
        mSignalLengthNextPowerOf2 = signalLengthNextPowerOf2;
        if (doFFT || signalLengthNextPowerOf2 == signalLength)
        {
            // An FFT is desired but the input signal isn't the right length
            if (signalLengthNextPowerOf2 != signalLength)
            {
                orderWork = orderWork + 1; 
                signalLengthNextPowerOf2
                    = static_cast<int> (std::pow(2, orderWork));
            }
            if (signalLengthNextPowerOf2 < signalLength)
            {
                throw std::invalid_argument(
                    "Failed to compute pad2 length for signal length "
                   + std::to_string(signalLength));
            }
            doFFT = true;
            order = orderWork;
        }
        int transformLength
            = (signalLength%2 == 0) ? signalLength/2 + 1 :
              (signalLength - 1)/2 + 1;
        if (doFFT)
        {
            transformLength = signalLengthNextPowerOf2/2 + 1;
        } 
        // There's a chance we can leave early
        if (static_cast<int> (doFFT) == static_cast<int> (mDoFFT) &&
            signalLength == mSignalLength &&
            transformLength == mTransformLength &&
            mAllocatedBuffers)
        {
            return;
        }
        // Ensure everything is cleared prior to reallocation
        clear();
        // Figure out sizes 
        int specSize, specBufferSize, bufferSize;
        if (doFFT)
        {
            ::getFFTSize<T>(&specSize,
                         &specBufferSize, // Init only
                          &bufferSize,
                          order,
                          ippAlgHintNone);
            mFFTSpec64f = nullptr;
            Ipp8u *pSpecBuffer = ippsMalloc_8u(specBufferSize);
            Ipp8u *pSpec = ippsMalloc_8u(specSize);
            mBuffer = ippsMalloc_8u(bufferSize);
            if constexpr (std::is_same<T, double>::value)
            {
                auto status = ippsFFTInit_R_64f(&mFFTSpec64f,
                                                order,
                                                IPP_FFT_DIV_INV_BY_N,
                                                ippAlgHintNone,
                                                pSpec,
                                                pSpecBuffer);
                if (pSpecBuffer){ippsFree(pSpecBuffer);}
                if (status != ippStsNoErr)
                {
                    ippsFree(pSpec);
                    clear();
                    throw std::runtime_error("Failed to initialize double FFT");
                }
            }
            else if constexpr (std::is_same<T, float>::value)
            {
                auto status = ippsFFTInit_R_32f(&mFFTSpec32f,
                                                order,
                                                IPP_FFT_DIV_INV_BY_N,
                                                ippAlgHintNone,
                                                pSpec,
                                                pSpecBuffer);
                if (pSpecBuffer){ippsFree(pSpecBuffer);}
                if (status != ippStsNoErr)
                {
                    ippsFree(pSpec);
                    clear();
                    throw std::runtime_error("Failed to initialize float FFT");
                }
            }
            else
            {
                static_assert(false, "Unhandled FFT precision");
            }
        }
        else
        {
            getDFTSize<T>(&specSize,
                          &specBufferSize, // Init only
                          &bufferSize,
                          signalLength,
                          ippAlgHintNone);
            mDFTSpec64f = nullptr;
            Ipp8u *pSpecBuffer = ippsMalloc_8u(specBufferSize);
            mBuffer = ippsMalloc_8u(bufferSize);
            if constexpr (std::is_same<T, double>::value)
            {
                mDFTSpec64f
                    = reinterpret_cast<IppsDFTSpec_R_64f *>
                      (ippsMalloc_8u(specSize));
                auto status = ippsDFTInit_R_64f(signalLength,
                                                IPP_FFT_DIV_INV_BY_N,
                                                ippAlgHintNone,
                                                mDFTSpec64f, 
                                                pSpecBuffer);
                if (pSpecBuffer){ippsFree(pSpecBuffer);}
                if (status != ippStsNoErr)
                {   
                    clear();
                    throw std::runtime_error("Failed to initialize double DFT");
                }   
            }
            else if constexpr (std::is_same<T, float>::value)
            {
                mDFTSpec32f
                    = reinterpret_cast<IppsDFTSpec_R_32f *>
                      (ippsMalloc_8u(specSize));
                auto status = ippsDFTInit_R_32f(signalLength,
                                                IPP_FFT_DIV_INV_BY_N,
                                                ippAlgHintNone,
                                                mDFTSpec32f,
                                                pSpecBuffer);
                if (pSpecBuffer){ippsFree(pSpecBuffer);}
                if (status != ippStsNoErr)
                {
                    clear();
                    throw std::runtime_error("Failed to initialize float DFT");
                }
            }
            else
            {   
                static_assert(false, "Unhandled DFT precision");
            }   
        }
        // Copy some results
        mSignalLength = signalLength;
        mTransformLength = transformLength;
        mDoFFT = doFFT;
        mAllocatedBuffers = true;
    }
    void apply(const USignal::Vector<T> &input,
               USignal::Vector<std::complex<T>> *output)
    {
        if (!mAllocatedBuffers)
        {
            throw std::runtime_error("Buffers not allocated");
        }
        if (static_cast<int> (output->size()) != getTransformLength())
        {
            constexpr std::complex<T> zero{0 + 2i};
            output->resize(getTransformLength(), zero);
        }
        if (mDoFFT)
        {
            if constexpr (std::is_same<T, double>::value)
            {
                const double *inputPtr = input.data();
                USignal::Vector<double> workSpace;
                // Pad?
                if (mSignalLength != mSignalLengthNextPowerOf2)
                {
                    workSpace.resize(mSignalLengthNextPowerOf2, 0);
                    std::copy(input.begin(), input.end(), workSpace.begin());
                    inputPtr = workSpace.data();
                }
                auto pSrc = reinterpret_cast<const Ipp64f *> (inputPtr);
                auto pDst = reinterpret_cast<Ipp64f *> (output->data());
                // Result for even length signals is:
                //   {R_0, 0}, {R_1, I_1}, ..., {R_{{N-1}/2}, I_{{N-1}/2}}
                // and for odd length signals:
                //   {R_0, 0}, {R_1, I_1}, ..., {R_{N/2-1}, I_{N/2-1}}
                auto status
                    = ippsFFTFwd_RToCCS_64f(pSrc, pDst, mFFTSpec64f, mBuffer);
                if (status != ippStsNoErr)
                {
                    throw std::runtime_error("Failed to compute double fft");
                }
            }
            else if constexpr (std::is_same<T, float>::value)
            {
                const float *inputPtr = input.data();
                USignal::Vector<float> workSpace;
                // Pad?
                if (mSignalLength != mSignalLengthNextPowerOf2)
                {
                    workSpace.resize(mSignalLengthNextPowerOf2, 0);
                    std::copy(input.begin(), input.end(), workSpace.begin());
                    inputPtr = workSpace.data();
                }
                auto pSrc = reinterpret_cast<const Ipp32f *> (inputPtr);
                auto pDst = reinterpret_cast<Ipp32f *> (output->data());
                auto status
                    = ippsFFTFwd_RToCCS_32f(pSrc, pDst, mFFTSpec32f, mBuffer);
                if (status != ippStsNoErr)
                {
                    throw std::runtime_error("Failed to compute float fft");
                }
            }
            else
            {
                throw std::runtime_error("Unhandled FFT implementation");
            }
        }
        else
        {
            if constexpr (std::is_same<T, double>::value)
            {
                const auto pSrc
                    = reinterpret_cast<const Ipp64f *> (input.data());
                auto pDst = reinterpret_cast<Ipp64f *> (output->data());
                auto status
                    = ippsDFTFwd_RToCCS_64f(pSrc, pDst, mDFTSpec64f, mBuffer);
                if (status != ippStsNoErr) 
                {
                    throw std::runtime_error("Failed to compute double dft");
                }
            }
            else if constexpr (std::is_same<T, float>::value)
            {
                const auto pSrc
                    = reinterpret_cast<const Ipp32f *> (input.data());
                auto pDst = reinterpret_cast<Ipp32f *> (output->data());
                auto status
                    = ippsDFTFwd_RToCCS_32f(pSrc, pDst, mDFTSpec32f, mBuffer);
                if (status != ippStsNoErr)
                {
                    throw std::runtime_error("Failed to compute float dft");
                }
            }
            else
            {
                throw std::runtime_error("Unhandled DFT implementation");
            }
        } 
    }

    [[nodiscard]] int getTransformLength() const
    {
        return mTransformLength;
    }
//private:
    ForwardParameters mParameters;
    IppsFFTSpec_R_64f *mFFTSpec64f{nullptr};
    IppsDFTSpec_R_64f *mDFTSpec64f{nullptr};
    IppsFFTSpec_R_32f *mFFTSpec32f{nullptr};
    IppsDFTSpec_R_32f *mDFTSpec32f{nullptr};
    Ipp8u *mBuffer{nullptr};
    int mSignalLength{-1};
    int mSignalLengthNextPowerOf2{-1};
    int mTransformLength{-1};
    bool mAllocatedBuffers{false};
    bool mDoFFT{false};
    bool mInitialized{false};
};

/// Constructor
template<typename T>
Forward<T>::Forward(const ForwardParameters &parameters) :
    pImpl(std::make_unique<ForwardImpl> (parameters))
{
}

template<typename T>
void Forward<T>::apply()
{
    if (!isInitialized())
    {
        throw std::invalid_argument("Forward transformer not initialized");
    }
    const auto &xRef = this->getInputReference();
    if (xRef.empty()){throw std::invalid_argument("Input signal is empty");}
    pImpl->initialize(static_cast<int> (xRef.size())); // throws
    // Allocate the output
    USignal::Vector<std::complex<T>> y;
    pImpl->apply(xRef, &y);
    this->setOutput(std::move(y)); 
}

/// Initialized?
template<typename T>
bool Forward<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Destructor
template<typename T>
Forward<T>::~Forward() = default;

///--------------------------------------------------------------------------///
///                                    Parameters                            ///
///--------------------------------------------------------------------------///

class ForwardParameters::ForwardParametersImpl
{
public:
    ForwardParameters::Implementation
       mImplementation{ForwardParameters::Implementation::Discrete};
};

/// Constructor
ForwardParameters::ForwardParameters() :
    pImpl(std::make_unique<ForwardParametersImpl> ())
{
}

/// Copy constructor
ForwardParameters::ForwardParameters(const ForwardParameters &parameters)
{
    *this = parameters;
}

/// Move constructor
ForwardParameters::ForwardParameters(ForwardParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Constructor
ForwardParameters::ForwardParameters(const Implementation implementation) :
    pImpl(std::make_unique<ForwardParametersImpl> ())
{
    setImplementation(implementation);
}

/// Copy assignment
ForwardParameters&
ForwardParameters::operator=(const ForwardParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<ForwardParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move assignment
ForwardParameters&
ForwardParameters::operator=(ForwardParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Implementation
void ForwardParameters::setImplementation(
    const Implementation implementation) noexcept
{
    pImpl->mImplementation = implementation;
}

ForwardParameters::Implementation
ForwardParameters::getImplementation() const noexcept
{
    return pImpl->mImplementation;
}

/// Destructor
ForwardParameters::~ForwardParameters() = default;

///--------------------------------------------------------------------------///
///                           Template Instantiation                         ///
///--------------------------------------------------------------------------///
template class USignal::Transforms::Fourier::Forward<double>;
template class USignal::Transforms::Fourier::Forward<float>;
