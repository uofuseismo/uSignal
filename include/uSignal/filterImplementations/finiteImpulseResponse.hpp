#ifndef USIGNAL_FILTER_IMPLEMENTATIONS_FINITE_IMPULSE_RESPONSE_HPP
#define USIGNAL_FILTER_IMPLEMENTATIONS_FINITE_IMPULSE_RESPONSE_HPP
#include <memory>
#include <uSignal/filterRepresentations/finiteImpulseResponse.hpp>
#include <uSignal/system/system.hpp>
namespace USignal::FilterImplementations
{

/// @brief Implements the finite impulse response filter.
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
    /// @param[in] filterCoefficients  The finite impulse response filter
    ///                                coefficients.
    explicit FiniteImpulseResponse(const USignal::FilterRepresentations::FiniteImpulseResponse<T> &filterCoefficients,
                                   Implementation implementation = Implementation::Direct);
    /// @result True indicates the 
    [[nodiscard]] bool isInitialized() const noexcept final;   

    void apply() final;
    /// @brief Destructor.
    ~FiniteImpulseResponse() override;

    FiniteImpulseResponse() = delete;
private:
    class FiniteImpulseResponseImpl;
    std::unique_ptr<FiniteImpulseResponseImpl> pImpl;
};

}
#endif
