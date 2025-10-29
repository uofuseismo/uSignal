#ifndef USIGNAL_FILTER_IMPLEMENTATIONS_FINITE_IMPULSE_RESPONSE_HPP
#define USIGNAL_FILTER_IMPLEMENTATIONS_FINITE_IMPULSE_RESPONSE_HPP
#include <memory>
#include <uSignal/filterRepresentations/finiteImpulseResponse.hpp>
#include <uSignal/system/system.hpp>
namespace USignal::FilterImplementations
{

/// @class FiniteImpulseResponse "finiteImpulseResponse.hpp"
/// @brief Implements the finite impulse response filter for post-processing.
/// @copyright Ben Baker (UUSS) distributed under the MIT No AI license.
template<class T>
class FiniteImpulseResponse final : public USignal::System::ISystem<T, T>
{
public:
    enum class Implementation
    {
        Direct,   /*!< Use a direct, time-domain implementation. */
        FFT,      /*!< Use a frequency domain, Fast Fourier transform implementation. */
        Automatic /*!< The the algorithm choose. */
    };
public:
    /// @brief Constructor.
    /// @param[in] filterCoefficients  The finite impulse response filter
    ///                                coefficients.
    /// @param[in] implementation      Influences the choice of filter implementation.
    /// @param[in] isRealTime          If false then this filter is for post-prcoessing.
    ///                                In this case, initial conditions will be restored
    ///                                after every application.
    explicit FiniteImpulseResponse(const USignal::FilterRepresentations::FiniteImpulseResponse<T> &filterCoefficients,
                                   Implementation implementation = Implementation::Direct,
                                   const bool isRealTime = false);

    /// @result True indicates the class is ready for application.
    [[nodiscard]] bool isInitialized() const noexcept final;

    /// @brief Sets the filter initial conditions.  By default the system
    ///        is assumed to be at rest - i.e., all 0.
    /// @throws std::invalid_argument if initialConditions.size() does not
    ///         equal the filter order.
    void setInitialConditions(const USignal::Vector<T> &initialConditions);

    /// @brief Filters the signal.
    void apply() final;

    /// @brief Resets the delay lines to the initial conditions.
    void resetInitialConditions();

    /// @brief Destructor.
    ~FiniteImpulseResponse() override;

    FiniteImpulseResponse() = delete;
    FiniteImpulseResponse(const FiniteImpulseResponse &) = delete;
    FiniteImpulseResponse(FiniteImpulseResponse &&) noexcept = delete;
private:
    class FiniteImpulseResponseImpl;
    std::unique_ptr<FiniteImpulseResponseImpl> pImpl;
};

}
#endif
