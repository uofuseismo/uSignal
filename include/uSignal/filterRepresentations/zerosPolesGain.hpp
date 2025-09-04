#ifndef USIGNAL_FILTER_REPRESENTATIONS_ZEROS_POLES_GAIN_HPP
#define USIGNAL_FILTER_REPRESENTATIONS_ZEROS_POLES_GAIN_HPP
#include <uSignal/vector.hpp>
#include <complex>
#include <memory>
namespace USignal::FilterRepresentations
{
template<class T = double>
class ZerosPolesGain
{
public:
    /// @brief Initializes the system in terms of its zeros, poles, and gain.
    /// @param[in] zeros  The zeros of the system.
    /// @param[in] poles  The poles of the system. 
    /// @param[in] gain   The system gain.
    ZerosPolesGain(const USignal::Vector<std::complex<T>> &zeros,
                   const USignal::Vector<std::complex<T>> &poles,
                   T gain);
    /// @brief Copy constructor.
    ZerosPolesGain(const ZerosPolesGain &zerosPolesGain);
    /// @brief Move constructor.
    ZerosPolesGain(ZerosPolesGain &&zerosPolesGain) noexcept;

    /// @result The poles of the system.
    [[nodiscard]] USignal::Vector<std::complex<T>> getPoles() const;
    /// @result The poles  of the system.
    [[nodiscard]] const USignal::Vector<std::complex<T>> &getPolesReference() const;
    /// @result The zeros of the system.
    [[nodiscard]] USignal::Vector<std::complex<T>> getZeros() const;
    /// @result The zeros of the system.
    [[nodiscard]] const USignal::Vector<std::complex<T>> &getZerosReference() const;
    /// @result The system gain.
    [[nodiscard]] T getGain() const;
      
    /// @brief destructor.
    ~ZerosPolesGain();

    /// @brief Copy assignment.
    ZerosPolesGain& operator=(const ZerosPolesGain &zerosPolesGain);
    /// @brief Move assignment.
    ZerosPolesGain& operator=(ZerosPolesGain &&zerosPolesGain) noexcept;

    ZerosPolesGain() = delete;
   
private:
    class ZerosPolesGainImpl;
    std::unique_ptr<ZerosPolesGainImpl> pImpl;
};
}
#endif
