#ifndef PRIVATE_INFINITY_NORM_HPP
#define PRIVATE_INFINITY_NORM_HPP
#include <cmath>
#include <uSignal/vector.hpp>
namespace
{

template<typename T>
T computeInfinityNorm(const USignal::Vector<T> &yTrue,
                      const USignal::Vector<T> &yEstimate)
{
    if (yTrue.size() != yEstimate.size())
    {
        throw std::invalid_argument("Inconsistent input sizes");
    }
    T error{0};
    for (int i = 0; i < static_cast<int> (yTrue.size()); ++i)
    {
        error = std::max(error, std::abs(yTrue[i] - yEstimate[i]));
    }
    return error;
}

}
#endif
