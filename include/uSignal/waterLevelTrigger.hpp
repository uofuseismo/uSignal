#ifndef USIGNAL_WATER_LEVEL_TRIGGER_HPP
#define USIGNAL_WATER_LEVEL_TRIGGER_HPP
#include <memory>
#include <uSignal/system/system.hpp>
namespace USignal
{
/// @brief A water-level trigger works as follows: When a characteristic
///        signal exceeds some on tolerance, the window trigger window
///        is activated.  Then, when the characteristic signal falls below
///        an off tolerance the trigger window is terminated.
class WaterLevelTriggerOptions
{
public:
    /// @param[in] triggerOnAndOffThreshold  The trigger window commences
    ///               when the characteristic signal exceeds 
    ///               triggerOnAndOffThreshold.first and the characteristic
    ///               signals falls below triggerOnAndOffThreshold.second. 
    explicit WaterLevelTriggerOptions(const std::pair<double, double> &triggerOnAndOffThreshold);
    /// @brief Copy constructor.
    WaterLevelTriggerOptions(const WaterLevelTriggerOptions &options);
    /// @brief Move constructor.
    WaterLevelTriggerOptions(WaterLevelTriggerOptions &&options) noexcept;
        
    /// @result When the signal exceeds result.first the detection
    ///         window commences.  When the signal then drops below
    ///         result.second the detection window ends.
    [[nodiscard]] std::pair<double, double> getOnAndOffThreshold() const noexcept;

    /// @brief Copy assignment operator.
    WaterLevelTriggerOptions& operator=(const WaterLevelTriggerOptions &options);
    /// @brief Move assignment operator.
    WaterLevelTriggerOptions& operator=(WaterLevelTriggerOptions &&options) noexcept;

    /// @brief Destructor.
    ~WaterLevelTriggerOptions();
private:
    class WaterLevelTriggerOptionsImpl;
    std::unique_ptr<WaterLevelTriggerOptionsImpl> pImpl;
};

template<class T>
class WaterLevelTrigger final : public USignal::System::ISystem<T, int>
{
public:
    explicit WaterLevelTrigger(const WaterLevelTriggerOptions &options,
                               const bool realTime = false);

    /// @result True indicates the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;

    /// @brief Applies the water level trigger.
    void apply() final;

    /// @brief Resets the trigger to a state where it is awaiting the
    ///        characteristic signal to exceed the on threshold.
    void resetInitialConditions();

    /// @brief Destructor.
    ~WaterLevelTrigger() final;
private:
    class WaterLevelTriggerImpl;
    std::unique_ptr<WaterLevelTriggerImpl> pImpl;
};

}
#endif
