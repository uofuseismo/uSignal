#include <iostream>
#include <algorithm>
#include <complex>
#include <cmath>
#include <vector>
#include <boost/align.hpp>
#include "uSignal/vector.hpp"
#include "src/alignment.hpp"

using namespace USignal;

template<class T>
class Vector<T>::VectorImpl
{
public:
    VectorImpl() = default;
    explicit VectorImpl(size_t n) :
        mX(n)
    {
    }
    VectorImpl(size_t n, T value) :
        mX(n, value)
    {
    }
    std::vector<T, boost::alignment::aligned_allocator<T, ALIGNMENT>> mX;
    //std::vector<double, boost::alignment::aligned_allocator<double, ALIGNMENT> > vector;
};

/// Constructor
template<class T>
Vector<T>::Vector() :
    pImpl(std::make_unique<VectorImpl> ())
{
}

/// Copy constructor
template<class T>
Vector<T>::Vector(const Vector<T> &v)
{
    *this = v;
}

/// Move constructor
template<class T>
Vector<T>::Vector(Vector<T> &&v) noexcept
{
    *this = std::move(v);
}

/// Construct from vector
template<class T>
Vector<T>::Vector(const std::vector<T> &v) :
    pImpl(std::make_unique<VectorImpl> ())
{
    pImpl->mX.resize(v.size());
    std::copy(v.begin(), v.end(), pImpl->mX.begin()); 
}

/// Construct vector of a given size
template<class T>
Vector<T>::Vector(const size_t n) :
    pImpl(std::make_unique<VectorImpl> (n))
{
}

/// Construct vector of a given size with default values
template<class T>
Vector<T>::Vector(const size_t n, const T value) :
    pImpl(std::make_unique<VectorImpl> (n, value))
{
}

/// Copy assignment
template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T> &v)
{
    if (&v == this){return *this;}
    pImpl = std::make_unique<VectorImpl> (*v.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Vector<T>& Vector<T>::operator=(Vector<T> &&v) noexcept
{
    if (&v == this){return *this;}
    pImpl = std::move(v.pImpl);
    return *this;
}

/// Resize
template<class T>
void Vector<T>::resize(const size_t n)
{
    pImpl->mX.resize(n);
}

template<class T>
void Vector<T>::resize(const size_t n, const T value)
{
    pImpl->mX.resize(n, value);
}

/// Reserve
template<class T>
void Vector<T>::reserve(const size_t n)
{
    pImpl->mX.reserve(n);
}

/// Data
template<class T>
T* Vector<T>::data() noexcept
{
    return pImpl->mX.data();
}

template<class T>
const T* Vector<T>::data() const noexcept
{
    return pImpl->mX.data();
}

/// Size
template<class T>
size_t Vector<T>::size() const noexcept
{
    return pImpl->mX.size();
}

/// Release memory
template<class T>
void Vector<T>::clear() noexcept
{
    pImpl->mX.clear();
}

/// Empty?
template<class T>
bool Vector<T>::empty() const noexcept
{
    return pImpl->mX.empty();
}

/// Pop front
template<class T>
void Vector<T>::pop_front()
{
    // Awkward but I want to keep that memory alignment
    std::shift_left(pImpl->mX.begin(), pImpl->mX.end(), 1);
    pImpl->mX.pop_back();
}

/// Pop back
template<class T>
void Vector<T>::pop_back()
{
    pImpl->mX.pop_back();
}

/// Push back
template<class T>
void Vector<T>::push_back(const T value)
{
    pImpl->mX.push_back(value);
}

/// Destructor
template<class T>
Vector<T>::~Vector() = default;

/// Iterators
template<class T>
Vector<T>::iterator Vector<T>::begin()
{
    return pImpl->mX.begin();
}

template<class T>
Vector<T>::iterator Vector<T>::end()
{
    return pImpl->mX.end();
}

template<class T>
Vector<T>::const_iterator Vector<T>::begin() const
{
    return pImpl->mX.begin();
}

template<class T>
Vector<T>::const_iterator Vector<T>::end() const
{
    return pImpl->mX.end();
}

template<class T>
Vector<T>::const_iterator Vector<T>::cbegin() const
{
    return pImpl->mX.cbegin();
}

template<class T>
Vector<T>::const_iterator Vector<T>::cend() const
{
    return pImpl->mX.cend();
}

template<class T>
T& Vector<T>::operator[](const size_t index)
{
    return pImpl->mX[index];
}

template<class T>
T& Vector<T>::operator[](const size_t index) const
{
    return pImpl->mX[index];
}

template<class T>
T& Vector<T>::at(size_t index)
{
    return pImpl->mX.at(index);
}

template<class T>
T& Vector<T>::at(const size_t index) const
{
    return pImpl->mX.at(index);
}

template<typename T>
USignal::Vector<T> 
USignal::operator+(const T a, const Vector<T> &x)
{
    auto n = static_cast<int> (x.size());
    USignal::Vector<T> y;
    if (n == 0){return y;}
    y.resize(x.size());
    const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
    auto yPtr = std::assume_aligned<ALIGNMENT> (y.data());
    std::transform(xPtr, xPtr + n, yPtr,
                   [=](const auto xi)
                   {
                       return a + xi;
                   });
    return y;
}

template<typename T>
USignal::Vector<T> 
USignal::operator+(const Vector<T> &x, const T a) 
{
    return a + x;
}

template<typename T>
USignal::Vector<T> 
USignal::operator*(const T a, const Vector<T> &x) 
{
    auto n = static_cast<int> (x.size());
    USignal::Vector<T> y;
    if (n == 0){return y;} 
    y.resize(x.size());
    const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
    auto yPtr = std::assume_aligned<ALIGNMENT> (y.data());
    std::transform(xPtr, xPtr + n, yPtr,
                   [=](const auto xi)
                   {
                       return a*xi;
                   });
    return y;
}

template<typename T>
USignal::Vector<std::complex<T>> 
USignal::operator*(const std::complex<T> a, const Vector<std::complex<T>> &x) 
{
    auto n = static_cast<int> (x.size());
    USignal::Vector<std::complex<T>> y;
    if (n == 0){return y;} 
    y.resize(x.size());
    const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
    auto yPtr = std::assume_aligned<ALIGNMENT> (y.data());
    std::transform(xPtr, xPtr + n, yPtr,
                   [=](const auto xi) 
                   {   
                       return a*xi;
                   }); 
    return y;
}

template<typename T>
USignal::Vector<T> 
USignal::operator*(const Vector<T> &x, const T a)
{   
    return a*x;
}   

template<typename T>
USignal::Vector<T>
USignal::operator/(const Vector<T> &x, const T a)
{
    constexpr T zero{0};
    if (a == zero){throw std::invalid_argument("Division by zero");}
/*
    auto n = static_cast<int> (x.size());
    USignal::Vector<T> y;
    if (n == 0){return y;} 
    y.resize(x.size());
    const auto xPtr = std::assume_aligned<ALIGNMENT> (x.data());
    auto yPtr = std::assume_aligned<ALIGNMENT> (y.data());
    std::transform(xPtr, xPtr + n, yPtr,
                   [=](const auto xi) 
                   {
                       return xi/a;
                   });
    return y;
*/
    // Typically accurate enough
    auto aInverse = static_cast<T> (1.0/a);
    return aInverse*x;
}

template class USignal::Vector<double>;
template class USignal::Vector<float>;
template class USignal::Vector<int>;
template class USignal::Vector<std::complex<double>>;
template class USignal::Vector<std::complex<float>>;
 
template USignal::Vector<double> USignal::operator+(const USignal::Vector<double> &x, double a);
template USignal::Vector<float> USignal::operator+(const USignal::Vector<float> &x, float a);
template USignal::Vector<double> USignal::operator*(const USignal::Vector<double> &x, double a); 
template USignal::Vector<float> USignal::operator*(const USignal::Vector<float> &x, float a); 
template USignal::Vector<std::complex<double>> USignal::operator*(std::complex<double> a, const USignal::Vector<std::complex<double>> &x);
template USignal::Vector<std::complex<float>> USignal::operator*(std::complex<float> a, const USignal::Vector<std::complex<float>> &x);
template USignal::Vector<std::complex<double>> USignal::operator*(const USignal::Vector<std::complex<double>> &x, std::complex<double> a);
template USignal::Vector<std::complex<float>> USignal::operator*(const USignal::Vector<std::complex<float>> &x, std::complex<float> a); 
template USignal::Vector<double> USignal::operator/(const USignal::Vector<double> &x, double a); 
template USignal::Vector<float> USignal::operator/(const USignal::Vector<float> &x, float a); 

