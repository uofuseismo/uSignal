#ifndef USIGNAL_FILTER_IMPLEMENTATIONS_SECOND_ORDER_SECTIONS_HPP
#define USIGNAL_FILTER_IMPLEMENTATIONS_SECOND_ORDER_SECTIONS_HPP
#include <memory>
#include <uSignal/filterRepresentations/secondOrderSections.hpp>
#include <uSignal/system/system.hpp>
namespace USignal::FilterImplementations
{

/// @brief Implements the cascade of second order section (biquad) infinite
///        impulse response filters.  This algorithm should be preferred
//         for high-order IIR filters.
/// @copyright Ben Baker (UUSS) distributed under the MIT No AI license.
template<class T>
class SecondOrderSections final : public USignal::System::ISystem<T, T>
{
public:
    /// @brief Constructor.
    /// @param[in] filterCoefficients   The second order section filter
    ///                                 coefficients.
    explicit SecondOrderSections(const USignal::FilterRepresentations::SecondOrderSections<T> &filterCoefficients, const bool isRealTime = false);

    /// @result True indicates the 
    [[nodiscard]] bool isInitialized() const noexcept final;   

    /// @brief Filters the signal.
    void apply() final;

    /// @brief Resets the delay lines to the initial conditions.
    void resetInitialConditions();

    /// @brief Destructor.
    ~SecondOrderSections() override;

    SecondOrderSections() = delete;
private:
    class SecondOrderSectionsImpl;
    std::unique_ptr<SecondOrderSectionsImpl> pImpl;
};

}
#endif
