#ifndef USIGNAL_TRANSFORMS_HILBERT_FINITE_IMPULSE_RESPONSE_HPP
#define USIGNAL_TRANSFORMS_HILBERT_FINITE_IMPULSE_RESPONSE_HPP
#include <complex>
#include <memory>
#include <uSignal/vector.hpp>
#include <uSignal/system/system.hpp>

namespace USignal::Transforms::Hilbert
{

class FiniteImpulseResponseOptions
{
public:
    /// @brief Constructor.
    FiniteImpulseResponseOptions();
    /// @brief Copy constructor.
    FiniteImpulseResponseOptions(const FiniteImpulseResponseOptions &options);
    /// @brief Move constructor.
    FiniteImpulseResponseOptions(FiniteImpulseResponseOptions &&options) noexcept;
    
    /// @brief Sets the FIR filter order. 
    /// @param[in] order  The filter order.
    ///                   If order is even then the real FIR filter will have
    ///                   one non-zero coefficient and every other imaginary
    ///                   coefficient will be non-zero.  The Type III filter
    ///                   is computationally advantageous.
    ///                   If order is odd then then neither the real nor
    ///                   imaginary FIR filters will be sparse.  However, the
    ///                   filter amplitude near the Nyquist frequency,
    ///                   \f$ \pi \f$, will be constant.  The Type IV filter
    ///                   produces a better approximation to Hilbert transform.
    /// @note The number of filter taps is the order + 1.
    void setOrder(uint16_t order) noexcept;
    /// @result The filter order.
    /// @result By default this is 300.
    [[nodiscard]] int getOrder() const noexcept;

    /// @brief Sets the beta in the FIR Kaiser window filter design.
    /// @param[in] beta  The beta in design.
    /// @throws std::invalid_argument if beta is not positive. 
    void setBeta(double beta);
    /// @result The beta in the FIR Kaiser window filter design.
    /// @note By default this is 8.
    [[nodiscard]] double getBeta() const noexcept;

    /// @brief Destructor.
    ~FiniteImpulseResponseOptions();

    /// @brief Copy assignment.
    FiniteImpulseResponseOptions& operator=(const FiniteImpulseResponseOptions &options);
    /// @brief Move assignment.
    FiniteImpulseResponseOptions& operator=(FiniteImpulseResponseOptions &&options) noexcept;
private:
    class FiniteImpulseResponseOptionsImpl;
    std::unique_ptr<FiniteImpulseResponseOptionsImpl> pImpl;
};

/// @class FiniteImpulseResopnse finiteImpulseResponse.hpp
/// @brief Computes the analytic signal via an FIR Hilbert transform.
/// @note The Hilbert transform is obtained from the analytic signal, 
///       \f$ x_a = x_r + i x_h \f$ where \f$ x_r \f$ is the input signal and
///       \f$ x_h \f$ is the Hilbert transform of the signal.
template<class T>
class FiniteImpulseResponse final : public USignal::System::ISystem<T, std::complex<T>>
{
public:
    explicit FiniteImpulseResponse(const FiniteImpulseResponseOptions &options);
 
    /// @brief Destructor.
    ~FiniteImpulseResponse() override;
    FiniteImpulseResponse() = delete;
private:
    class FiniteImpulseResponseImpl;
    std::unique_ptr<FiniteImpulseResponseImpl> pImpl;
};

}
#endif
