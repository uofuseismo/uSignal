#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <mkl.h>
/*
#include "rtseis/enums.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"
*/
#include "uSignal/utilities/interpolation/weightedAverageSlopes.hpp"
#include "src/alignment.hpp"

#define MKL_SPLINE_ORDER 4

using namespace USignal::Utilities::Interpolation;

namespace
{

#pragma omp declare simd
template<class T>
inline void computeWiWiMi(const T mi, T *wi, T *wimi)
{
    constexpr T huge = std::numeric_limits<T>::max();
    constexpr T eps = std::numeric_limits<T>::epsilon();
    // The goal is to compute w_i and w_i m_i = 1/max(|m_i|, epsilon)*m_i
    // The trick is that to let w_i m_i safely go to 0 and let the
    // denominator become really big when the slope is 0.
    *wi = huge;
    *wimi = 0;
    auto ami = std::abs(mi);
    if (ami > eps)
    {
        *wi = 1/ami;
        // "If m_i and m_{i+1} have opposite signs then s_i is zero".
        // Realize, the copysign in the ensuing numerator will ensure
        // operations like 1 - 1, 1 + 1, -1 + 1, or -1 - 1.
        *wimi = std::copysign(1, mi);
    }
}

/// Implements first equation in Wiggins' Interpolation fo Digitized Curves
/// pg. 2077
template<typename T>
void computeUniformSlopes(const T dx,
                          const USignal::Vector<T> &y,
                          USignal::Vector<T> *splineCoeffs)
{
    auto n = static_cast<int> (y.size());
#ifndef NDEBUG
    assert(n > 1);
#endif
    USignal::Vector<T> slopes(n, 0); 
    const auto yPtr = std::assume_aligned<ALIGNMENT> (y.data());
    auto slopesPtr = std::assume_aligned<ALIGNMENT> (slopes.data());
    constexpr T one{1};
    T dxi = one/dx;
    // Handle the initial conditions
    slopesPtr[0] = (yPtr[1] - yPtr[0])/dx;
    for (int i = 1; i < n - 1; ++i)
    {
        T mi  = (yPtr[i] - yPtr[i-1])*dxi;
        T mi1 = (yPtr[i+1] - yPtr[i])*dxi;
        // w_i and w_i m_i = 1/max(|m_i|, epsilon)*m_i
        T wi, wimi;
        computeWiWiMi(mi, &wi, &wimi);
        // w_{i+1} and w_{i+1} m_{i+1} = 1/(max(|m_{i+1}|, epsilon)*m_{i+1}
        T wi1, wi1mi1;
        computeWiWiMi(mi1, &wi1, &wi1mi1);
        // s_i = (w_i*m_i + w_{i+1}*m_{i+1})/(w_{i} + w_{i+1})
        slopesPtr[i] = (wimi + wi1mi1)/(wi + wi1);
    }
    // Handle the final conditions
    slopes[n - 1] = (yPtr[n - 1] - yPtr[n - 2])/dx;
    // Compute the spline coefficients: Eqn 4 from
    // Monotone Piecewise Cubic Interpolation - Fritsch and Carlson 1980
    auto splineCoeffsPtr
         = std::assume_aligned<ALIGNMENT> (splineCoeffs->data());
    auto dxi2 = one/(dx*dx);
    for (int i = 0; i < n - 1; ++i)
    {
        auto di = slopesPtr[i];
        auto di1 = slopesPtr[i + 1];
        auto delta = (yPtr[i + 1] - yPtr[i])*dxi;
        auto c1 = yPtr[i];
        auto c2 = di;
        auto c3 = (-2*di - di1 + 3*delta)*dxi;
        auto c4 = (di + di1 - 2*delta)*dxi2;
        splineCoeffsPtr[4*i + 0] = c1; 
        splineCoeffsPtr[4*i + 1] = c2;
        splineCoeffsPtr[4*i + 2] = c3;
        splineCoeffsPtr[4*i + 3] = c4;
    }
}

template<typename T>
void computeUniformSlopes(const USignal::Vector<T> &x,
                          const USignal::Vector<T> &y,
                          USignal::Vector<T> *splineCoeffs)
{
    auto n = static_cast<int> (x.size());
#ifndef NDEBUG
    assert(x.size() == y.size());
    assert(n > 1);
#endif
    const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
    const auto yPtr = std::assume_aligned<ALIGNMENT> (y.data());
    USignal::Vector<T> slopes(n, 0);
    auto slopesPtr = std::assume_aligned<ALIGNMENT> (slopes.data());
    // Handle the initial conditions
    slopesPtr[0] = (yPtr[1] - yPtr[0])/(xPtr[1] - xPtr[0]);
    for (int i = 1; i < n - 1; ++i)
    {
        T mi  = (yPtr[i] - yPtr[i - 1])/(xPtr[i] - xPtr[i - 1]);
        T mi1 = (yPtr[i + 1] - yPtr[i])/(xPtr[i + 1] - xPtr[i]);
        // w_i and w_i m_i = 1/max(|m_i|, epsilon)*m_i
        T wi, wimi;
        computeWiWiMi(mi, &wi, &wimi);
        // w_{i+1} and w_{i+1} m_{i+1} = 1/(max(|m_{i+1}|, epsilon)*m_{i+1}
        T wi1, wi1mi1;
        computeWiWiMi(mi1, &wi1, &wi1mi1);
        // s_i = (w_i*m_i + w_{i+1}*m_{i+1})/(w_{i} + w_{i+1})
        slopesPtr[i] = (wimi + wi1mi1)/(wi + wi1);
    }
    // Handle the final conditions
    slopesPtr[n - 1] = (yPtr[n - 1] - yPtr[n - 2])/(xPtr[n - 1] - xPtr[n - 2]);
    // Compute the spline coefficients: Eqn 4 from
    // Monotone Piecewise Cubic Interpolation - Fritsch and Carlson 1980
    auto splineCoeffsPtr
         = std::assume_aligned<ALIGNMENT> (splineCoeffs->data());
    constexpr T one{1};
    for (int i = 0; i < n - 1; ++i)
    {
        T dx = xPtr[i+1] - xPtr[i]; 
        auto dxi = one/dx; 
        auto dxi2 = dxi*dxi;
        auto di = slopesPtr[i];
        auto di1 = slopesPtr[i + 1];
        auto delta = (yPtr[i + 1] - yPtr[i])*dxi;
        auto c1 = yPtr[i];
        auto c2 = di;
        auto c3 = (-2*di - di1 + 3*delta)*dxi;
        auto c4 = (di + di1 - 2*delta)*dxi2;
        splineCoeffsPtr[4*i + 0] = c1;
        splineCoeffsPtr[4*i + 1] = c2;
        splineCoeffsPtr[4*i + 2] = c3;
        splineCoeffsPtr[4*i + 3] = c4;
    }
}

}

template<class T>
class WeightedAverageSlopes<T>::WeightedAverageSlopesImpl
{
public:
    /// Copy assignment operator
    WeightedAverageSlopesImpl&
    operator=(const WeightedAverageSlopesImpl &slopes)
    {
        if (&slopes == this){return *this;}
        clear();
        mXiEqual[0] = slopes.mXiEqual[0];
        mXiEqual[1] = slopes.mXiEqual[1];
        mRange = slopes.mRange;
        mCoeffs = slopes.mCoeffs;
        mSites = slopes.mSites;
        mUniformPartition = slopes.mUniformPartition;
        mInitialized = slopes.mInitialized; 
        if (mCoeffs > 0)
        {
            mSplineCoeffs = slopes.mSplineCoeffs;
        }
        if (!mUniformPartition && mSites > 0)
        {
            mXi = slopes.mXi; 
        }
        return *this;
    }
    /// Destructor
    ~WeightedAverageSlopesImpl()
    {
        clear(); 
    }
    /// Releases memory on the class
    void clear() noexcept
    {
        if (mHaveTask64f){dfDeleteTask(&mTask64f);}
        if (mHaveTask32f){dfDeleteTask(&mTask32f);}
        mSplineCoeffs.clear();
        mXi.clear();
        mTask64f = nullptr;
        mTask32f = nullptr;
        mXiEqual[0] = 0;
        mXiEqual[1] = 0;
        mRange = std::make_pair(0, 0);
        mCoeffs = 0;
        mSites = 0;
        mUniformPartition = true;
        mHaveTask64f = false;
        mHaveTask32f = false;
        mInitialized = false;        
    }
//private: 
    /// The task pointers for the data fitting 
    DFTaskPtr mTask64f{nullptr};
    DFTaskPtr mTask32f{nullptr};
    /// The spline coefficients.  This has dimension [mCoeffs].
    USignal::Vector<T> mSplineCoeffs;
    /// The abscissas for non-uniform inteprolation. This has dimension [mSites]
    USignal::Vector<T> mXi;
    /// The range for an equal interpolation
    std::array<T, 2> mXiEqual{0, 0};
    /// The min/max x for interpolation
    std::pair<double, double> mRange{0, 0};
    /// The number of spline coefficients.  This is splineOrder*(mSites - 1)
    int mCoeffs = 0;
    /// The number of interpolation sites
    int mSites = 0;
    /// The custom spline order
    //const MKL_INT splineOrder = 4;
    /// Notes the precision of the module
    //RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
    /// Using a uniform partition 
    bool mUniformPartition{true};
    /// Have the double MKL task.
    bool mHaveTask64f{false};
    /// Have the float MKL task.
    bool mHaveTask32f{false};
    /// Class is initialized.
    bool mInitialized{false};
};
/*
/// Constructors
template<class T>
WeightedAverageSlopes<T>::WeightedAverageSlopes() :
    pImpl(std::make_unique<WeightedAverageSlopesImpl> ())
{
}

/// Copy constructor
template<class T>
WeightedAverageSlopes<T>::WeightedAverageSlopes(
    const WeightedAverageSlopes &slopes)
{
    *this = slopes;
}

/// Move constructor
template<class T>
WeightedAverageSlopes<T>::WeightedAverageSlopes(
    WeightedAverageSlopes &&slopes) noexcept
{
    *this = std::move(slopes);
}

/// Copy assignment operator
template<class T>
WeightedAverageSlopes<T>& WeightedAverageSlopes<T>::operator=(
    const WeightedAverageSlopes &slopes)
{
    if (&slopes == this){return *this;}
    pImpl = std::make_unique<WeightedAverageSlopesImpl> (*slopes.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
WeightedAverageSlopes<T>& WeightedAverageSlopes<T>::operator=(
    WeightedAverageSlopes &&slopes) noexcept
{
    if (&slopes == this){return *this;}
    pImpl = std::move(slopes.pImpl);
    return *this;
}
*/

/// Destructors
template<class T>
WeightedAverageSlopes<T>::~WeightedAverageSlopes() = default;

/*
/// Releases memory on the module
template<class T>
void WeightedAverageSlopes<T>::clear() noexcept
{
    pImpl->clear();
}
*/

/// Determines if the class is initialized
template<class T>
bool WeightedAverageSlopes<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Initialize the class
template<class T>
WeightedAverageSlopes<T>::WeightedAverageSlopes(
    const std::pair<T, T> &endPoints,
    const USignal::Vector<T> &y)
{
    if (y.size() < 2)
    {
        throw std::invalid_argument("At least two points required");
    }
    if (endPoints.first >= endPoints.second)
    {
        throw std::invalid_argument("endPoints.first = "
                                  + std::to_string(endPoints.first)
                                  + " must be less than endPoints.second = "
                                  + std::to_string(endPoints.second));
    }
    //clear();
    // Compute the spline coefficients
    auto n = static_cast<int> (y.size());
    auto dx = static_cast<T> ((endPoints.second - endPoints.first)/static_cast<double> (n - 1));
    pImpl->mSites = n;
    pImpl->mCoeffs = MKL_SPLINE_ORDER*(pImpl->mSites - 1);
    pImpl->mSplineCoeffs.resize(pImpl->mCoeffs, 0);
    ::computeUniformSlopes(dx, y, &pImpl->mSplineCoeffs);
    // Create a custom piecewise 4th order spline
    pImpl->mTask64f = nullptr;
    pImpl->mRange.first = endPoints.first;
    pImpl->mRange.second = endPoints.second;
    pImpl->mXiEqual[0] = endPoints.first;
    pImpl->mXiEqual[1] = endPoints.second;
    if (std::is_same<T, double>::value)
    {
        auto *xiEqual
             = reinterpret_cast<const double *> (pImpl->mXiEqual.data());
        auto *yPtr = reinterpret_cast<const double *> (y.data());
        auto status = dfdNewTask1D(&pImpl->mTask64f, n, xiEqual,
                                   DF_UNIFORM_PARTITION, 1, yPtr, DF_NO_HINT); 
        if (status != DF_STATUS_OK)
        {
            dfDeleteTask(&pImpl->mTask64f);
            throw std::runtime_error("Failed to create task");
        }
        auto *splineCoeffsPtr
            = reinterpret_cast<const double *> (pImpl->mSplineCoeffs.data());
        status = dfdEditPPSpline1D(
            pImpl->mTask64f, MKL_SPLINE_ORDER,
            DF_PP_DEFAULT, DF_NO_BC, NULL, DF_NO_IC, NULL,
            splineCoeffsPtr, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            dfDeleteTask(&pImpl->mTask64f);
            throw std::runtime_error("Failed to edit spline pipeline");
        }
        pImpl->mHaveTask64f = true;
    }
    else if (std::is_same<T, float>::value)
    {
        auto *xiEqual 
             = reinterpret_cast<const float *> (pImpl->mXiEqual.data());
        auto *yPtr = reinterpret_cast<const float *> (y.data());
        auto status = dfsNewTask1D(&pImpl->mTask32f, n, xiEqual,
                                   DF_UNIFORM_PARTITION, 1, yPtr, DF_NO_HINT); 
        if (status != DF_STATUS_OK)
        {   
            dfDeleteTask(&pImpl->mTask32f);
            throw std::runtime_error("Failed to create task");
        }   
        auto *splineCoeffsPtr 
            = reinterpret_cast<const float *> (pImpl->mSplineCoeffs.data());
        status = dfsEditPPSpline1D(
            pImpl->mTask32f, MKL_SPLINE_ORDER,
            DF_PP_DEFAULT, DF_NO_BC, NULL, DF_NO_IC, NULL,
            splineCoeffsPtr, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            dfDeleteTask(&pImpl->mTask32f);
            throw std::runtime_error("Failed to edit spline pipeline");
        }
        pImpl->mHaveTask32f = true;
    }
    else
    {
        throw std::runtime_error("Unhandled precision");
    }
    pImpl->mInitialized = true;
}

/// Initialize the class
template<class T>
WeightedAverageSlopes<T>::WeightedAverageSlopes(
    const USignal::Vector<T> &x,
    const USignal::Vector<T> &y)
{
    if (x.size() < 2)
    {
        throw std::invalid_argument("There must be at least 2 abscissas");
    }
    if (x.size() != y.size())
    {
        throw std::invalid_argument("x.size() = "
                                  + std::to_string(x.size())
                                  + " must equal y.size() = "
                                  + std::to_string(y.size()));
    }
    if (!std::is_sorted(x.begin(), x.end()))
    {
        throw std::invalid_argument("x abscissas are not sorted");
    }
    //clear();
    auto n = static_cast<int> (y.size()); 
    pImpl->mSites = n;
    pImpl->mCoeffs = MKL_SPLINE_ORDER*(pImpl->mSites - 1);
    pImpl->mSplineCoeffs.resize(pImpl->mCoeffs, 0);
    computeUniformSlopes(x, y, &pImpl->mSplineCoeffs);
    // Create a custom piecewise 4th order spline
    pImpl->mTask64f = nullptr;
    pImpl->mRange.first = x.at(0);
    pImpl->mRange.second = x.at(n - 1);
    pImpl->mXiEqual[0] = x.at(0);
    pImpl->mXiEqual[1] = x.at(n - 1);
    pImpl->mXi = x;
    if (std::is_same<T, double>::value)
    {
        auto *xiPtr
             = reinterpret_cast<const double *> (pImpl->mXi.data());
        auto *yPtr = reinterpret_cast<const double *> (y.data());
        auto status = dfdNewTask1D(&pImpl->mTask64f, n, xiPtr,
                                   DF_NON_UNIFORM_PARTITION,
                                   1, yPtr, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            dfDeleteTask(&pImpl->mTask64f);
            throw std::runtime_error("Failed to create task");
        }
        auto *splineCoeffsPtr
            = reinterpret_cast<const double *> (pImpl->mSplineCoeffs.data());
        status = dfdEditPPSpline1D(pImpl->mTask64f, MKL_SPLINE_ORDER,
                                   DF_PP_DEFAULT, DF_NO_BC, NULL, DF_NO_IC, NULL,
                                   splineCoeffsPtr, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            dfDeleteTask(&pImpl->mTask64f);
            throw std::runtime_error("Failed to edit spline pipeline");
        }
        pImpl->mHaveTask64f = true;
    }
    else if (std::is_same<T, float>::value)
    {
        auto *xiPtr
             = reinterpret_cast<const float *> (pImpl->mXi.data());
        auto *yPtr = reinterpret_cast<const float *> (y.data());
        auto status = dfsNewTask1D(&pImpl->mTask32f, n, xiPtr,
                                   DF_NON_UNIFORM_PARTITION,
                                   1, yPtr, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            dfDeleteTask(&pImpl->mTask32f);
            throw std::runtime_error("Failed to create task");
        }
        auto *splineCoeffsPtr
            = reinterpret_cast<const float *> (pImpl->mSplineCoeffs.data());
        status = dfsEditPPSpline1D(pImpl->mTask32f, MKL_SPLINE_ORDER,
                                   DF_PP_DEFAULT, DF_NO_BC, NULL, DF_NO_IC, NULL,
                                   splineCoeffsPtr, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            dfDeleteTask(&pImpl->mTask32f);
            throw std::runtime_error("Failed to edit spline pipeline");
        }
        pImpl->mHaveTask32f = true;
    }
    else
    {
        throw std::runtime_error("Unhandled precision");
    }
    pImpl->mInitialized = true;
}

// Interpolate
template<class T>
USignal::Vector<T> WeightedAverageSlopes<T>::interpolate(
    const USignal::Vector<T> &xq) const
{
    USignal::Vector<T> yq;
    auto nq = static_cast<int> (xq.size());
    // Checks
    if (nq < 1){return yq;} // Nothing to do
    double xMin = getMinimumX(); // Throws on initialization
    double xMax = getMaximumX(); // Throws on initialization
    auto [xqMin, xqMax] = std::minmax_element(xq.begin(), xq.end());
    if (*xqMin < xMin || *xqMax > xMax)
    {
        throw std::invalid_argument("Min/max of xq = ("
                                  + std::to_string(*xqMin) + ","
                                  + std::to_string(*xqMax) 
                                  + ") Must be in range ["
                                  + std::to_string(xMin) + ","
                                  + std::to_string(xMax) + "]");
    }
    auto isSorted = std::is_sorted(xq.begin(), xq.end());
    // Interpolate
    const MKL_INT nsite{nq};
    MKL_INT sortedHint = isSorted ? DF_SORTED_DATA : DF_NO_HINT;
    constexpr MKL_INT nOrder{1};   // Length of dorder
    const MKL_INT dOrder[1] = {0}; // Order of derivatives
    if (std::is_same<T, double>::value)
    {   
        const auto xqPtr = reinterpret_cast<const double *> (xq.data());
        auto yqPtr = reinterpret_cast<double *> (yq.data());
        auto status = dfdInterpolate1D(pImpl->mTask64f, DF_INTERP, DF_METHOD_PP,
                                       nsite, xqPtr,
                                       sortedHint, nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yqPtr,
                                       DF_MATRIX_STORAGE_ROWS, NULL);
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Interpolation failed");
        }
    }
    else
    {
        const auto xqPtr = reinterpret_cast<const float *> (xq.data());
        auto yqPtr = reinterpret_cast<float *> (yq.data());
        auto status = dfsInterpolate1D(pImpl->mTask64f, DF_INTERP, DF_METHOD_PP,
                                       nsite, xqPtr,
                                       sortedHint, nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yqPtr,
                                       DF_MATRIX_STORAGE_ROWS, NULL);
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Interpolation failed");
        }
    }
    return yq;
}

// Uniform interpolation
template<class T>
USignal::Vector<T> WeightedAverageSlopes<T>::interpolate(
    const int nq, const std::pair<T, T> &xInterval) const
{
    USignal::Vector<T> yq;
    // Checks
    if (nq < 1){return yq;} // Nothing to do
    auto xMin = getMinimumX(); // Throws on initialization
    auto xMax = getMaximumX(); // Throws on initialization
    yq.resize(nq, 0);
    double xqMin = xInterval.first;
    double xqMax = xInterval.second;
    if (xqMin > xqMax)
    {
        throw std::invalid_argument("xInterval.first = "
                                  + std::to_string(xqMin)
                                  + " > xInterval.second = "
                                  + std::to_string(xqMax));
    }
    if (xqMin < xMin || xqMax > xMax)
    {
        throw std::invalid_argument("Min/max of xq = ("
                                  + std::to_string(xqMin) + ","
                                  + std::to_string(xqMax)
                                  + ") must be in range ["    
                                  + std::to_string(xMin) + ","
                                  + std::to_string(xMax) + "]");
    }
    // Interpolate
    const MKL_INT nsite{nq};
    constexpr MKL_INT sortedHint{DF_UNIFORM_PARTITION};
    constexpr MKL_INT nOrder{1};  // Length of dorder
    const MKL_INT dOrder[1] = {0}; // Order of derivatives
    std::array<T, 2> xq;
    xq[0] = xqMin;
    xq[1] = xqMax;
    if (std::is_same<T, double>::value)
    {
        const auto xqPtr = reinterpret_cast<double *> (xq.data());
        auto yqPtr = reinterpret_cast<double *> (yq.data());
        auto status = dfdInterpolate1D(pImpl->mTask64f, DF_INTERP, DF_METHOD_PP,
                                       nsite, xqPtr,
                                       sortedHint, nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yqPtr,
                                       DF_MATRIX_STORAGE_ROWS, NULL);
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Interpolation failed");
        }
    }
    else if (std::is_same<T, float>::value)
    {
        const auto xqPtr = reinterpret_cast<float *> (xq.data());
        auto yqPtr = reinterpret_cast<float *> (yq.data());
        auto status = dfsInterpolate1D(pImpl->mTask32f, DF_INTERP, DF_METHOD_PP,
                                       nsite, xqPtr,
                                       sortedHint, nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yqPtr,
                                       DF_MATRIX_STORAGE_ROWS, NULL);
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Interpolation failed");
        }
    }
    else
    {
        throw std::invalid_argument("Unhandled precision");
    }
    return yq;
}
/// Get minimum x
template<class T>
T WeightedAverageSlopes<T>::getMinimumX() const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    return pImpl->mRange.first;
}

/// Get maximum x
template<class T>
T WeightedAverageSlopes<T>::getMaximumX() const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    return pImpl->mRange.second;
}

/// Template class instantiation
template class USignal::Utilities::Interpolation::WeightedAverageSlopes<double>;
template class USignal::Utilities::Interpolation::WeightedAverageSlopes<float>;
