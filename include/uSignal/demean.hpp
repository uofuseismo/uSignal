#ifndef USIGNAL_DEMEAN_HPP
#define USIGNAL_DEMEAN_HPP
#include <memory>
#include <uSignal/system/system.hpp>
namespace USignal
{
template<class T = double>
/// @class Demean
/// @brief Removes the mean from a signal.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Demean final : public USignal::System::ISystem<T, T>
{
public:
    /// @brief Constructor.
    Demean();
    /// @brief Copy constuctor.
    Demean(const Demean &demean);
    /// @brief Move constructor.
    Demean(Demean &&demean) noexcept;

    /// @brief Copy assignment.
    Demean& operator=(const Demean &demean);
    /// @brief Move assignment.
    Demean& operator=(Demean &&demean) noexcept;
    

    /// @result True indicates the class is initialized which, in this case,
    ///         is always true.
    [[nodiscard]] bool isInitialized() const noexcept override;

    /// @brief Demeans the signal.
    void apply() override;

    /// @result The mean of the input signal.
    [[nodiscard]] double getMean() const noexcept;
 
    void clear() noexcept override;
    /// @brief Constructor.
    ~Demean() override;
private:
    class DemeanImpl;
    std::unique_ptr<DemeanImpl> pImpl;
};
}
#endif
