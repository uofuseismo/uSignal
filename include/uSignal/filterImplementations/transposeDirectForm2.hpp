#ifndef USIGNAL_FILTER_IMPLEMENTATIONS_TRANSPOSE_DIRECT_FORM_2_HPP
#define USIGNAL_FILTER_IMPLEMENTATIONS_TRANSPOSE_DIRECT_FORM_2_HPP
#include <memory>
#include <uSignal/filterRepresentations/infiniteImpulseResponse.hpp>
#include <uSignal/system/system.hpp>
namespace USignal::FilterImplementations
{

/// @brief Implements the Tranpose Direct Form 2 infinite impulse response
///        filter implementation.  This is a relatively slow algorithm but
///        has good numerical properties.  However, for high-order IIR filters
///        you should be looking at Second Order Section filters.
/// @copyright Ben Baker (UUSS) distributed under the MIT No AI license.
template<class T>
class TransposeDirectForm2 final : public USignal::System::ISystem<T, T>
{
public:
    /// @param[in] filterCoefficients  The infinite impulse repsonse filter
    ///                                coefficients.
    explicit TransposeDirectForm2(const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filterCoefficients);
    /// @result True indicates the 
    [[nodiscard]] bool isInitialized() const noexcept final;   
    void apply() final;
    /// @brief Destructor.
    ~TransposeDirectForm2() override;

    TransposeDirectForm2() = delete;
private:
    class TransposeDirectForm2Impl;
    std::unique_ptr<TransposeDirectForm2Impl> pImpl;
};

}
#endif
