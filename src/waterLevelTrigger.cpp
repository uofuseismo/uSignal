#include <limits>
#include <cmath>
#include "uSignal/waterLevelTrigger.hpp"
#include "uSignal/vector.hpp"

using namespace USignal;

class WaterLevelTriggerOptions::WaterLevelTriggerOptionsImpl
{
public:
    double mOnThreshold{0};
    double mOffThreshold{0};
};

/// Constructor
WaterLevelTriggerOptions::WaterLevelTriggerOptions(
    const std::pair<double, double> &onAndOffThreshold) :
    pImpl(std::make_unique<WaterLevelTriggerOptionsImpl> ())
{
    pImpl->mOnThreshold = onAndOffThreshold.first;
    pImpl->mOffThreshold = onAndOffThreshold.second;
}

/// Copy constructor
WaterLevelTriggerOptions::WaterLevelTriggerOptions(
    const WaterLevelTriggerOptions &options)
{
    *this = options;
}

/// Move constructor
WaterLevelTriggerOptions::WaterLevelTriggerOptions(
    WaterLevelTriggerOptions &&options) noexcept
{
    *this = std::move(options);
}

/// Copy assignment
WaterLevelTriggerOptions&
WaterLevelTriggerOptions::operator=(const WaterLevelTriggerOptions &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<WaterLevelTriggerOptionsImpl> (*options.pImpl);
    return *this;
}

/// Move assignment
WaterLevelTriggerOptions&
WaterLevelTriggerOptions::operator=(WaterLevelTriggerOptions &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}

/// On and off threshold
std::pair<double, double> 
WaterLevelTriggerOptions::getOnAndOffThreshold() const noexcept
{
    return std::pair {pImpl->mOnThreshold, pImpl->mOffThreshold};
}

/// Destructor
WaterLevelTriggerOptions::~WaterLevelTriggerOptions() = default;

///--------------------------------------------------------------------------///
///                              Trigger                                     ///
///--------------------------------------------------------------------------///

template<class T>
class WaterLevelTrigger<T>::WaterLevelTriggerImpl
{
public:
    std::unique_ptr<WaterLevelTriggerOptions> mOptions{nullptr};
    double mOnThreshold{0};
    double mOffThreshold{0};
    T mDelayLine{std::numeric_limits<T>::lowest()};
    bool mIsRealTime{false};
    bool mArmed{false};
    bool mInitialized{false};
};

/// Constructor
template<class T>
WaterLevelTrigger<T>::WaterLevelTrigger(
    const WaterLevelTriggerOptions &options,
    const bool isRealTime) :
    pImpl(std::make_unique<WaterLevelTriggerImpl> ())
{
    pImpl->mOptions = std::make_unique<WaterLevelTriggerOptions> (options);
    pImpl->mOnThreshold = pImpl->mOptions->getOnAndOffThreshold().first;
    pImpl->mOffThreshold = pImpl->mOptions->getOnAndOffThreshold().second;
    pImpl->mIsRealTime = isRealTime;
    pImpl->mInitialized = true;
}

/// Initialized?
template<class T>
bool WaterLevelTrigger<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Apply
template<class T>
void WaterLevelTrigger<T>::apply() 
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class not initialized");
    }
    const auto &x = this->getInputReference();
    auto nSamples = static_cast<int> (x.size());
    USignal::Vector<int> y;
    if (nSamples == 0)
    {
        this->setOutput(std::move(y));
        return;
    }
    y.resize(nSamples, 0);
    // Get to work
    auto onTolerance = pImpl->mOnThreshold;
    auto offTolerance = pImpl->mOffThreshold;
    bool armed = pImpl->mArmed;
    T xim1 = std::numeric_limits<T>::lowest();
    if (!pImpl->mIsRealTime){armed = false;}
    for (int i = 0; i < nSamples; ++i)
    {
        if (!armed)
        {
            // Searching for start of window
            if (xim1 < onTolerance && x[i] >= onTolerance)
            {
                armed = true;
            }
        }
        else
        {
            if (xim1 >= offTolerance && x[i] < offTolerance)
            {
                armed = false;
            }
        }
        xim1 = x[i];
        y[i] = static_cast<int> (armed);
    }
    if (pImpl->mIsRealTime)
    {
        pImpl->mDelayLine = xim1;
        pImpl->mArmed = armed;
    }
    this->setOutput(std::move(y));
}

/// Reset initial conditions
template<class T>
void WaterLevelTrigger<T>::resetInitialConditions()
{
    if (isInitialized()){throw std::runtime_error("Trigger not initialized");}
    pImpl->mDelayLine = std::numeric_limits<T>::lowest();
    pImpl->mArmed = false;
}

/// Destructor
template<class T>
WaterLevelTrigger<T>::~WaterLevelTrigger() = default;

///--------------------------------------------------------------------------///
///                        Template instantiation                            ///
///--------------------------------------------------------------------------///
template class USignal::WaterLevelTrigger<double>;
template class USignal::WaterLevelTrigger<float>;

