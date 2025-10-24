#ifndef USIGNAL_TRANSFORMS_FOURIER_FORWARD_HPP
#define USIGNAL_TRANSFORMS_FOURIER_FORWARD_HPP
#include <complex>
#include <memory>
#include <uSignal/vector.hpp>
#include <uSignal/system/system.hpp>
namespace USignal::Transforms::Fourier
{

class ForwardOptions
{
public:
    /// @brief Defines the Fourier transform implementation.
    enum class Implementation
    {
        Discrete, /*!< Implements the discrete Fourier transform - i.e., no padding. */
        Fast      /*!< Implements the fast Fourier transform - i.e., pads signal to the next power of 2.
                       Note, for smaller signals it is typically faster to power through
                       a DFT calculation. */
    };
public:
    /// @brief Constructor.
    ForwardOptions();
    /// @brief Constructor with an implementation.
    explicit ForwardOptions(Implementation implementation);
    /// @brief Copy constructor.
    /// @param[in] options  The parameters from which to initialize this
    ///                        class.
    ForwardOptions(const ForwardOptions &options);
    /// @brief Move constructor.
    /// @param[in,out] options  The parameters from which to initialize this
    ///                            class.  On exit, options's behavior is
    ///                            undefined.
    ForwardOptions(ForwardOptions &&options) noexcept;

    /// @brief Defines the implementation.
    void setImplementation(Implementation implementation) noexcept;
    /// @result The DFT implementation.  By default this is Discrete.
    [[nodiscard]] Implementation getImplementation() const noexcept;

    /// @brief Copy assignment.
    /// @param[in] options  The forward parameters to copy to this.
    /// @result A deep copy of the forward transform options.
    ForwardOptions& operator=(const ForwardOptions &options);
    /// @brief Move assignment.
    /// @param[in,out] options  The forward parameters whose memory will be
    ///                            moved this.  On exit, options's behavior
    ///                            is undefined.
    /// @result The memory from options moved to this.
    ForwardOptions& operator=(ForwardOptions &&options) noexcept;

    /// @brief Destructor.
    ~ForwardOptions(); 
private:
    class ForwardOptionsImpl;
    std::unique_ptr<ForwardOptionsImpl> pImpl;
};

/// @class Forward "forward.hpp"
/// @brief Performs the forward Fourier transform, here, defined as
/// TODO
/// @copyright Ben Baker (University of Utah) distributed under the
///            MIT NO AI license.
template<class T>
class Forward final : public USignal::System::ISystem<T, std::complex<T>>
{
public:
    /// @brief Initializes the forward Fourier transformer.
    /// @param[in] options  The options influencing the transform calculation.
    explicit Forward(const ForwardOptions &options);
 
    /// @result True indicates the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept override;

    /// @brief Computes the Fourier tarnsform of the input signal.
    void apply() override;

    /// @brief Destructor.
    ~Forward() override;

    Forward(const Forward &) = delete;
    Forward& operator=(const Forward &) = delete;
    Forward() = delete; 
private:
    class ForwardImpl;
    std::unique_ptr<ForwardImpl> pImpl;
};

}
#endif
