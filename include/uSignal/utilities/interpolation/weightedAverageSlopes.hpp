#ifndef USIGNAL_UTILITIES_INTERPOLATION_WEIGHTED_AVERAGE_SLOPES_HPP
#define USIGNAL_UTILITIES_INTERPOLATION_WEIGHTED_AVERAGE_SLOPES_HPP
#include <memory>
#include <uSignal/vector.hpp>
namespace USignal::Utilities::Interpolation
{

class WeightedAverageSlopesOptions
{
public:

};

template<class T>
class WeightedAverageSlopes //final : public USignal::System::ISystem<T, T>
{
public:
    WeightedAverageSlopes(const std::pair<T, T> &interval, const USignal::Vector<T> &y);
    WeightedAverageSlopes(const USignal::Vector<T> &abscissas, const USignal::Vector<T> &y);

    [[nodiscard]] USignal::Vector<T> interpolate(const USignal::Vector<T> &xq) const;
    [[nodiscard]] USignal::Vector<T> interpolate(int nq, const std::pair<T, T> &xInterval) const;

    /// @result True indicates the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @result The minimum x at which you can interoplate.
    [[nodiscard]] T getMinimumX() const;
    /// @result The maximum x at which you can interpolate.
    [[nodiscard]] T getMaximumX() const; 
    ~WeightedAverageSlopes();
private:
    class WeightedAverageSlopesImpl;
    std::unique_ptr<WeightedAverageSlopesImpl> pImpl;
};

}
#endif
