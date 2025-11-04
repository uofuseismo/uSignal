#include <string>
#include <complex>
#include <cmath>
#include "uSignal/filterDesign/infiniteImpulseResponse/convertBand.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/vector.hpp"

USignal::FilterRepresentations::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::
    convertLowpassAnalogPrototypeToLowpass(
        const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
        double omega0)
{
    if (omega0 <= 0)
    {
        throw std::invalid_argument("omega0 = " + std::to_string(omega0)
                                  + " must be positive");
    }
    const auto zeros = zpk.getZerosReference(); 
    const auto poles = zpk.getPolesReference();
    if (poles.size() < zeros.size())
    {
        throw std::runtime_error("BUG! npoles < nzeros not implemented");
    }
    // Scale all points radially from origin to shift cutoff frequency
    auto zerosLowpass = std::complex<double> (omega0, 0)*zeros;
    auto polesLowpass = std::complex<double> (omega0, 0)*poles;
    // Each shifted pole decreases gain by omega0 and each shifted zero
    // increases by omega0.  Cancel out the net change to keep overall
    // gain the same
    auto gain = static_cast<double> (zpk.getGain());
    auto degree = static_cast<int> (poles.size() - zeros.size());
    double gainLowpass = gain*std::pow(omega0, degree);
    USignal::FilterRepresentations::ZerosPolesGain<double>
        zpkLowpass{zerosLowpass, polesLowpass, gainLowpass};
    return zpkLowpass;
}

USignal::FilterRepresentations::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::
    convertLowpassAnalogPrototypeToHighpass(
        const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
        double omega0)
{
    if (omega0 <= 0)
    {   
        throw std::invalid_argument("omega0 = " + std::to_string(omega0)
                                  + " must be positive");
    }   
    const auto zeros = zpk.getZerosReference(); 
    const auto poles = zpk.getPolesReference();
    if (poles.size() < zeros.size())
    {   
        throw std::runtime_error("BUG! npoles < nzeros not implemented");
    }   
    // Invert positions radially about unit circle to convert LPF TO HPF
    // Scale all points radially from origin to shift cutoff frequency
    // Note, this initialization sets every element to zero which,
    // if the lowpass had zeros at infinity, this would invert them to
    // the origin.
    USignal::Vector<std::complex<double>>
        zerosHighpass(std::max(zeros.size(), poles.size()),
                      std::complex<double> (0, 0));
    std::complex<double> zerosProduct(1, 0);
    for (int i = 0; i < static_cast<int> (zeros.size()); ++i)
    {
        zerosHighpass[i] = omega0/zeros[i];
        zerosProduct = -zeros[i]*zerosProduct;
    }
    // Invert positions radially about unit circle to convert LPF TO HPF
    // Scale all points radially from origin to shift cutoff frequency
    USignal::Vector<std::complex<double>> polesHighpass(poles.size());
    std::complex<double> polesProduct(1, 0);
    for (int i = 0; i < static_cast<int> (poles.size()); ++i)
    {
        polesHighpass[i] = omega0/poles[i];
        polesProduct = -poles[i]*polesProduct;
    }
    // Compute scale factors and cancel out gain caused by inversion
    auto gain = static_cast<double> (zpk.getGain());
    double gainHighpass = gain*std::real(zerosProduct/polesProduct);
    USignal::FilterRepresentations::ZerosPolesGain<double>
         zpkHighpass{zerosHighpass, polesHighpass, gainHighpass};
    return zpkHighpass;
}

USignal::FilterRepresentations::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::
    convertLowpassAnalogPrototypeToBandpass(
        const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
        const std::pair<double, double> &passband)
{
    auto omega0 = passband.first;
    auto omega1 = passband.second;
    auto bandwidth = omega1 - omega0;
    if (omega0 <= 0)
    {
        throw std::invalid_argument("omega0 = " + std::to_string(omega0)
                                  + " must be positive");
    }
    if (omega1 <= omega0)
    {
        throw std::invalid_argument("omega1 = " + std::to_string(omega1)
                                   + " must be greater than omega0 = "
                                   + std::to_string(omega0));
    }
    const auto zeros = zpk.getZerosReference();
    const auto poles = zpk.getPolesReference();
    auto nZerosBandpass = zeros.size() + poles.size();
    auto nPolesBandpass = 2*poles.size();
    // Define some constants
    const std::complex<double> zw02(omega0*omega0, 0);
    const std::complex<double> zbw2(bandwidth/2.0, 0);
    // Duplicate and zeros and shift from baseband to +w0 and -w0
    constexpr std::complex<double> zero(0, 0);
    USignal::Vector<std::complex<double>> zerosBandpass(nZerosBandpass, zero);
    auto j = static_cast<int> (zeros.size());
    for (int i = 0; i < static_cast<int> (zeros.size()); ++i)
    {
        auto zlp = zbw2*zeros[i];     //Scale zeros to desired bw
        auto zlp2 = zlp*zlp;          //zlp^2
        auto zdif = zlp2 - zw02;      //zlp^2 - w0^2
        auto zsqrt = std::sqrt(zdif); //sqrt(zlp^2 - w0^2)
        zerosBandpass[i] = zlp + zsqrt; //zlp+sqrt(zlp^2-w0^2)
        zerosBandpass[j] = zlp - zsqrt; //zlp-sqrt(zlp^2-w0^2)
        j = j + 1;
    }
    // Move degree zeros to origin leaving degree zeros at infinity for BPF 
    /*
    for (int i = 2*zeros.size(); i < zeros_bp; i++)
    {
        zerosBandpass[i] = zero;
    }
    */
    // Duplicate poles and shift from baseband to +w0 and -w0
    USignal::Vector<std::complex<double>> polesBandpass(nPolesBandpass, zero);
    j = static_cast<int> (poles.size());
    for (int i = 0; i < static_cast<int> (poles.size()); ++i)
    {
        auto plp = zbw2*poles[i];     //Scale poles to desired bw
        auto plp2 = plp*plp;          //plp^2
        auto pdif = plp2 - zw02;      //plp^2 - w0^2
        auto psqrt = std::sqrt(pdif); //sqrt(plp^2 - w0^2)
        polesBandpass[i] = plp + psqrt;   //plp+sqrt(plp^2-w0^2)
        polesBandpass[j] = plp - psqrt;   //plp-sqrt(plp^2-w0^2)
        j = j + 1;
    }
    // Cancel out gain change from frequency scaling
    auto gain = static_cast<double> (zpk.getGain());
    auto degree = static_cast<int> (poles.size() - zeros.size());
    double gainBandpass = gain*std::pow(bandwidth, degree); //k*bw*bw^degree
    USignal::FilterRepresentations::ZerosPolesGain<double>
         zpkBandpass{zerosBandpass, polesBandpass, gainBandpass};
    return zpkBandpass;
}

USignal::FilterRepresentations::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::
    convertLowpassAnalogPrototypeToBandstop(
        const USignal::FilterRepresentations::ZerosPolesGain<double> &zpk,
        const std::pair<double, double> &passband)
{   
    auto omega0 = passband.first;
    auto omega1 = passband.second;
    auto bandwidth = omega1 - omega0;
    if (omega0 <= 0)
    {   
        throw std::invalid_argument("omega0 = " + std::to_string(omega0)
                                  + " must be positive");
    }
    if (omega1 <= omega0)
    {   
        throw std::invalid_argument("omega1 = " + std::to_string(omega1)
                                   + " must be greater than omega0 = "
                                   + std::to_string(omega0));
    }
    const auto zeros = zpk.getZerosReference();
    const auto poles = zpk.getPolesReference();
    auto nZerosBandstop = 2*poles.size();
    auto nPolesBandstop = 2*poles.size();
    // Define some constants
    const std::complex<double> zw02(omega0*omega0, 0);
    const std::complex<double> zbw2(bandwidth/2.0, 0);
    // Duplicate the zeros and shift from baseband to +w0 and -w0
    constexpr std::complex<double> zero(0, 0);
    USignal::Vector<std::complex<double>> zerosBandstop(nZerosBandstop, zero);
    std::complex<double> zProduct(1, 0); // Init product
    auto j = static_cast<int> (zeros.size());
    for (int i = 0; i < static_cast<int> (zeros.size()); ++i)
    {
        // Invert to a highpass filter with desired bandwidth
        auto zhp = zbw2/zeros[i];     // bw2/2/z
        auto zhp2 = zhp*zhp;          // zhp^2
        auto zdif = zhp2 - zw02;      // zhp^2 - w0^2
        auto zsqrt = std::sqrt(zdif); // sqrt(zhp^2 - w0^2)
        zerosBandstop[i] = zhp + zsqrt; // zhp+sqrt(zhp^2-w0^2)
        zerosBandstop[j] = zhp - zsqrt; // zhp-sqrt(zhp^2-w0^2)
        // Update product 
        zProduct = -zeros[i]*zProduct;
        j = j + 1;
    }
    // Move any zeros at infinity to center of stopband
    j = 2*static_cast<int> (zeros.size());
    for (int i = 0; i < static_cast<int> (poles.size() - zeros.size()); ++i)
    {
        zerosBandstop[j] = std::complex<double> (0, omega0);
        j = j + 1;
    }
    for (int i = 0; i < static_cast<int> (poles.size() - zeros.size()); ++i)
    {
        zerosBandstop[j] = std::complex<double> (0,-omega0);
        j = j + 1;
    }
    // Duplicate the poles and shift from baseband to +w0 and -w0
    std::complex<double> pProduct(1, 0);
    USignal::Vector<std::complex<double>> polesBandstop(nPolesBandstop, zero);
    j = static_cast<int> (poles.size());
    for (int i = 0; i < static_cast<int> (poles.size()); ++i)
    {
        // Invert to a highpass filter with desired bandwidth
        auto php = zbw2/poles[i];         // bw2/2/p
        auto php2 = php*php;          // php^2
        auto pdif = php2 - zw02;      // php^2 - w0^2
        auto psqrt = std::sqrt(pdif); // sqrt(php^2 - w0^2)
        polesBandstop[i] = php + psqrt; // php+sqrt(php^2-w0^2)
        polesBandstop[j] = php - psqrt; // php-sqrt(php^2-w0^2)
        // Update product
        pProduct = -poles[i]*pProduct;
        j = j + 1;
    }
    // Cancel out gain change caused by inversion
    auto gain = static_cast<double> (zpk.getGain());
    double gainBandstop = gain*std::real(zProduct/pProduct);
    USignal::FilterRepresentations::ZerosPolesGain<double>
         zpkBandstop{zerosBandstop, polesBandstop, gainBandstop};
    return zpkBandstop;
}
