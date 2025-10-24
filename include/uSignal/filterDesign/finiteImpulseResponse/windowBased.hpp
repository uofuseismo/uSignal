#ifndef USIGNAL_FILTER_DESIGN_FINITE_IMPULSE_RESPONSE_WINDOW_BASED_HPP
#define USIGNAL_FILTER_DESIGN_FINITE_IMPULSE_RESPONSE_WINDOW_BASED_HPP
#include <uSignal/filterRepresentations/finiteImpulseResponse.hpp>
namespace USignal::FilterDesign::FiniteImpulseResponse::WindowBased
{

enum class Window
{
    Hamming,
    Bartlett,
    Hanning,
    OptimalBlackman
};

template<typename T = double>
FilterRepresentations::FiniteImpulseResponse<T> 
    lowPass(int order, const T normalizedCutoffFrequency,
            Window window = Window::Hamming);

template<typename T = double>
FilterRepresentations::FiniteImpulseResponse<T> 
    highPass(int order, const T normalizedCutoffFrequency,
             Window window = Window::Hamming);

template<typename T = double>
FilterRepresentations::FiniteImpulseResponse<T> 
    bandPass(int order, 
             const std::pair<T, T> &lowAndHighNormalizedCutoffFrequencies,
             Window window = Window::Hamming);

template<typename T = double> 
FilterRepresentations::FiniteImpulseResponse<T>
    bandStop(int order, 
             const std::pair<T, T> &lowAndHighNormalizedCutoffFrequencies,
             Window window = Window::Hamming);
}
#endif
