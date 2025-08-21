#ifndef USIGNAL_FILTER_IMPLEMENTATIONS_INFINITE_IMPULSE_RESPONSE_HPP
#define USIGNAL_FILTER_IMPLEMENTATIONS_INFINITE_IMPULSE_RESPONSE_HPP
#include <memory>
#include <uSignal/filterRepresentations/infiniteImpulseResponse.hpp>
#include <uSignal/system/system.hpp>
namespace USignal::FilterImplementations
{

template<class T>
class InfiniteImpulseResponse final : public USignal::System::ISystem<T, T>
{
public:
    explicit InfiniteImpulseResponse(const USignal::FilterRepresentations::InfiniteImpulseResponse<T> &filterCoefficients);

    [[nodiscard]] bool isInitialized() const noexcept final;   
    void apply() final;

    ~InfiniteImpulseResponse() override;

    InfiniteImpulseResponse() = default;
private:
    class InfiniteImpulseResponseImpl;
    std::unique_ptr<InfiniteImpulseResponseImpl> pImpl;
};

}
#endif
