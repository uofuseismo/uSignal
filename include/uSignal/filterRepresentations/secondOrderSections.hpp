#ifndef USIGNAL_FILTER_REPRESENTATIONS_SECOND_ORDER_SECTIONS_HPP
#define USIGNAL_FILTER_REPRESENTATIONS_SECOND_ORDER_SECTIONS_HPP
#include <memory>
#include <array>
#include <uSignal/vector.hpp>
namespace USignal::FilterRepresentations
{
  template<class T> class ZerosPolesGain;
  template<class T> class InfiniteImpulseResponse;
}
namespace USignal::FilterRepresentations
{
/// @name SecondOrderSections "secondOrderSections.hpp"
/// @brief An IIR filter whose representation is given by a cascade
///        of second order (biquad) IIR filters.
/// @param[in] 
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<class T>
class SecondOrderSections
{
public:
    enum class PairingStrategy
    {
        None,
        Nearest
    };
public:
    /// @brief Constructs the second order sections representation of the filter.
    /// @param[in] numeratorCoefficients    These are the numerator coefficients,
    ///                                     typically denoted by B.  The number of
    ///                                     sections is given by the length of the
    ///                                     vector.
    /// @param[in] denominatorCoefficients  These are the denominator coefficients,
    ///                                     typically denoted by A.  The number of
    ///                                     sections is given by the length of the
    ///                                     vector.  This must match the length of
    ///                                     the numerator coefficients vector.
    ///                                     Additionally, coefficient of each section
    ///                                     cannot be zero.
    SecondOrderSections(
        const std::vector<std::array<T, 3>> &numeratorCoefficients,
        const std::vector<std::array<T, 3>> &denominatorCoefficients);
    /// @brief Constructs the second order sections representation of the filter.
    /// @param[in] numeratorCoefficients    These are the numerator coefficients,
    ///                                     typically denoted by B.  This vector
    ///                                     represents a [ns x 3] matrix in row 
    ///                                     major format. 
    /// @param[in] denominatorCoefficients  These are the denominator coefficients,
    ///                                     typically denoted by A.  This vector
    ///                                     represents a [ns x 3] matrix in row
    ///                                     major format.
    SecondOrderSections(
        const USignal::Vector<T> &numeratorCoefficients,
        const USignal::Vector<T> &denominatorCoefficients);

    explicit SecondOrderSections(const InfiniteImpulseResponse<T> &ba,
                                 PairingStrategy strategy = PairingStrategy::Nearest);
    explicit SecondOrderSections(const ZerosPolesGain<T> &zpk,
                                 PairingStrategy strategy = PairingStrategy::Nearest);
    /// @brief Copy constructor.
    SecondOrderSections(const SecondOrderSections &secondOrderSections);
    /// @brief Move constructor.
    SecondOrderSections(SecondOrderSections &&secondOrderSections) noexcept;

    /// @result A vector holding the numerator coefficients.  This has
    ///         dimension [ns x 3] and represents row major matrix where
    ///         ns is given by \c getNumberOfSections().
    [[nodiscard]] USignal::Vector<T> getNumeratorFilterCoefficients() const noexcept;
    /// @result A vector holding the denominator coefficients.  This has
    ///         dimension [ns x 3] and represents a row major matrix where
    ///         ns is given by \c getNumberOfSections().
    [[nodiscard]] USignal::Vector<T> getDenominatorFilterCoefficients() const noexcept;
    /// @result A reference to the vector holding the numerator coefficients.
    [[nodiscard]] const USignal::Vector<T> &getNumeratorFilterCoefficientsReference() const noexcept;
    /// @result A reference to the vector holding the denominator coefficients.
    [[nodiscard]] const USignal::Vector<T> &getDenominatorFilterCoefficientsReference() const noexcept;
    /// @result The number of second order sections.
    [[nodiscard]] int getNumberOfSections() const noexcept;

    /// @brief Destructor.
    ~SecondOrderSections();

    /// @brief Copy assignment.
    SecondOrderSections& operator=(const SecondOrderSections &sections);
    /// @brief Move assignment.
    SecondOrderSections& operator=(SecondOrderSections &&sections) noexcept;

    SecondOrderSections() = delete;
private:
    class SecondOrderSectionsImpl;
    std::unique_ptr<SecondOrderSectionsImpl> pImpl;
};
}
#endif
