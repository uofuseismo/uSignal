#ifndef USIGNAL_TAPER_HPP
#define USIGNAL_TAPER_HPP
#include <memory>
#include <uSignal/system/system.hpp>
namespace USignal
{
class TaperParameters
{
public:
    /// @brief Defines the window function used in the taper.
    enum class Window
    {
        Bartlett, /*!< Barlett (triangle) window */
        Blackman, /*!< Standard Blackman */
        Boxcar,   /*!< Boxcar */
        Hamming,  /*!< Hamming window */
        Hanning,  /*!< Hann(ing) window */
        Kaiser,   /*!< Kaiser window.  This requires specification of \f$ \beta \f$. */
        Sine      /*!< Sine window. */
    };  
public:
    /// @brief Constructor.
    /// @param[in] window      The window type.
    /// @param[in] percentage  The percentage of the signal to which the taper
    ///                        is applied.  For example, if this is 10 then the
    ///                        first 5 pct and last 5 pct of the signal will
    ///                        be tapered using the specified window type.
    /// @param[in] beta        If window is Kaiser then this is the
    ///                        \f$ \beta \f$ parameter. 
    TaperParameters(TaperParameters::Window window,
                    double percentage,
                    double beta = 0.5);
    /// @brief Copy constructor.
    TaperParameters(const TaperParameters &parameters);
    /// @brief Move constructor.
    TaperParameters(TaperParameters &&parameters) noexcept; 

    /// @result The window type.
    [[nodiscard]] TaperParameters::Window getWindow() const noexcept;
    /// @result The taper percentage.
    [[nodiscard]] double getPercentage() const noexcept;
    /// @result The \f$ \beta \f$ in the Kaiser window.
    /// @note By default this is 0.5.
    [[nodiscard]] double getBeta() const noexcept;

    /// @brief Copy assignment.
    TaperParameters& operator=(const TaperParameters &parameters);
    /// @brief Move assignment.
    TaperParameters& operator=(TaperParameters &&parameters) noexcept;

    /// @brief Destructor.
    ~TaperParameters();
    TaperParameters() = delete;
private:
    class TaperParametersImpl;
    std::unique_ptr<TaperParametersImpl> pImpl;
};

template<class T>
/// @class Taper "taper.hpp"
/// @brief Applies a window function to a signal.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Taper final : public USignal::System::ISystem<T, T>
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Taper(const TaperParameters &parameters);
    /// @brief Copy constuctor.
    Taper(const Taper &taper);
    /// @brief Move constructor.
    Taper(Taper &&taper) noexcept;
    /// @}

    /// @brief Initializes the filter.
    /// @param[in] window       The window type.
    /// @param[in] percentage   The percentage of the signal to which the taper
    ///                         is applied.  For example, if this is 10 then the
    ///                         first 5 pct and last 5 pct of the signal will
    ///                         be tapered using the specified window type.
    /// @throws std::invalid_argument if the percentage is out of the
    ///         range [0,100].
    //void initialize(Window window, double percentage, double beta = 0.5);
    /// @result True indicates the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept override;
    /// @result The window type.
    /// @throws std::runtime_error if \c isInitialized() is false.
    //[[nodiscard]] Window getWindow() const;
    /// @result The percentage to apply to the signal.
    /// @throws std:runtime_error if \c isInitialized() is false.
    //[[nodiscard]] double getPercentage() const;

    /// @brief The signal to window.
    /// @param[in] x  The signal to window.
    void setInput(const Vector<T> &x) override;

    /// @brief Tapers the signal.
    /// @throws std::runtime_error if \c isInitialized() is false.
    void apply() override;

    /// @name Destructors
    /// @{

    /// @brief Resets the class.
    void clear() noexcept override;
    /// @brief Constructor.
    ~Taper() override;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    Taper& operator=(const Taper &taper);
    /// @brief Move assignment.
    Taper& operator=(Taper &&taper) noexcept;
    /// @}
    Taper() = delete;
private:
    class TaperImpl;
    std::unique_ptr<TaperImpl> pImpl;
};
}
#endif
