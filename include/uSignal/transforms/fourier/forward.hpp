#ifndef USIGNAL_TRANSFORMS_FOURIER_FORWARD_HPP
#define USIGNAL_TRANSFORMS_FOURIER_FORWARD_HPP
#include <complex>
#include <memory>
#include <uSignal/vector.hpp>
#include <uSignal/system/system.hpp>
namespace USignal::Transforms::Fourier
{

class ForwardParameters
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
    ForwardParameters();
    /// @brief Constructor with an implementation.
    explicit ForwardParameters(Implementation implementation);
    /// @brief Copy constructor.
    /// @param[in] parameters  The parameters from which to initialize this
    ///                        class.
    ForwardParameters(const ForwardParameters &parameters);
    /// @brief Move constructor.
    /// @param[in,out] parameters  The parameters from which to initialize this
    ///                            class.  On exit, parameters's behavior is
    ///                            undefined.
    ForwardParameters(ForwardParameters &&parameters) noexcept;

    /// @brief Defines the implementation.
    void setImplementation(Implementation implementation) noexcept;
    /// @result The DFT implementation.  By default this is Discrete.
    [[nodiscard]] Implementation getImplementation() const noexcept;

    /// @brief Copy assignment.
    /// @param[in] parameters  The forward parameters to copy to this.
    /// @result A deep copy of the forward transform parameters.
    ForwardParameters& operator=(const ForwardParameters &parameters);
    /// @brief Move assignment.
    /// @param[in,out] parameters  The forward parameters whose memory will be
    ///                            moved this.  On exit, parameters's behavior
    ///                            is undefined.
    /// @result The memory from parameters moved to this.
    ForwardParameters& operator=(ForwardParameters &&parameters) noexcept;

    /// @brief Destructor.
    ~ForwardParameters(); 
private:
    class ForwardParametersImpl;
    std::unique_ptr<ForwardParametersImpl> pImpl;
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
    explicit Forward(const ForwardParameters &parameters);
 
    /// @result True indicates the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept override;

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
