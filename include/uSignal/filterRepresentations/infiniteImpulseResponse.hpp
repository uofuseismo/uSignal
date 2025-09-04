#ifndef USIGNAL_FILTER_REPRESENTATIONS_INFINITE_IMPULSE_RESPONSE_HPP
#define USIGNAL_FILTER_REPRESENTATIONS_INFINITE_IMPULSE_RESPONSE_HPP
#include <memory>
#include <array>
#include <uSignal/vector.hpp>
#include <uSignal/filterRepresentations/zerosPolesGain.hpp>
namespace USignal::FilterRepresentations
{
/// @name InfiniteImpulseResponse "infiniteImpulseResponse.hpp"
/// @brief An IIR filter whose representation is given by a numerator
///        and denominator coefficients.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<class T>
class InfiniteImpulseResponse
{
public:
    /// @brief Constructs the IIR representation of the filter.
    /// @param[in] numeratorCoefficients    The feed-forward numerator coefficients.
    ///                                     These are commonly denoted by b.
    /// @param[in] denominatorCoefficients  The feed-backward denominator coefficients.
    ///                                     These are commonly denoted by a.
    InfiniteImpulseResponse(
        const USignal::Vector<T> &numeratorCoefficients,
        const USignal::Vector<T> &denominatorCoefficients);
    /// @brief Constructs the IIR from its zeros, poles, and gain.
    explicit InfiniteImpulseResponse(const ZerosPolesGain<T> &zpk);
    /// @brief Copy constructor.
    InfiniteImpulseResponse(const InfiniteImpulseResponse &secondOrderSections);
    /// @brief Move constructor.
    InfiniteImpulseResponse(InfiniteImpulseResponse &&secondOrderSections) noexcept;

    /// @result A vector holding the numerator coefficients.
    [[nodiscard]] USignal::Vector<T> getNumeratorFilterCoefficients() const noexcept;
    /// @result A vector holding the denominator coefficients.
    [[nodiscard]] USignal::Vector<T> getDenominatorFilterCoefficients() const noexcept;
    /// @result A reference to the vector holding the numerator coefficients.
    [[nodiscard]] const USignal::Vector<T> &getNumeratorFilterCoefficientsReference() const noexcept;
    /// @result A reference to the vector holding the denominator coefficients.
    [[nodiscard]] const USignal::Vector<T> &getDenominatorFilterCoefficientsReference() const noexcept;
    /// @result The filter order (this is the max(b.size(), a.size()) - 1).
    [[nodiscard]] int getOrder() const noexcept;

    /// @brief Destructor.
    ~InfiniteImpulseResponse();

    /// @brief Copy assignment.
    InfiniteImpulseResponse& operator=(const InfiniteImpulseResponse &iir);
    /// @brief Move assignment.
    InfiniteImpulseResponse& operator=(InfiniteImpulseResponse &&iir) noexcept;

    InfiniteImpulseResponse() = delete;
private:
    class InfiniteImpulseResponseImpl;
    std::unique_ptr<InfiniteImpulseResponseImpl> pImpl;
};
}
#endif
