#include <iostream>
#include <cmath>
#include <complex>
#include <string>
#include "uSignal/filterDesign/infiniteImpulseResponse/analogPrototype.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/vector.hpp"

USignal::FilterRepresentations::ZerosPolesGain<double> 
USignal::FilterDesign::InfiniteImpulseResponse::AnalogPrototype::butterworth(
    const int order)
{
    if (order < 0 || order > 24)
    {
        if (order < 0)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                      + " must be non-negative");
        }
        throw std::invalid_argument("order = " + std::to_string(order)
                                  + " must be less than 25");
    }
    const int nPoles{order + 1};
    const int nZeros{0};
    USignal::Vector<std::complex<double>> poles(nPoles);
    USignal::Vector<std::complex<double>> zeros(nZeros);
    const double pi_twoni = M_PI/(2.0*static_cast<double> (order + 1));
    for (int i = 0; i < nPoles; i++)
    {
        auto xm = static_cast<double> (-(order + 1) + 1 + 2*i);
        std::complex<double> arg(0.0, pi_twoni*xm);
        poles[i] =-std::exp(arg);
    }
    constexpr double gain = 1.0; //k = real(prod(-p)) which is 1
    USignal::FilterRepresentations::ZerosPolesGain<double>
        zpk{zeros, poles, gain};
    return zpk;
}

USignal::FilterRepresentations::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::
    AnalogPrototype::chebyshevTypeI(const int order, const double ripple)
{
    if (order < 0)
    {
        throw std::invalid_argument("order = " + std::to_string(order)
                                 + " must be positive");
    }
    if (ripple <= 0)
    {
        throw std::invalid_argument("ripple = " + std::to_string(ripple)
                                  + " must be positive");
    }
    // Set space
    const int nPoles{order + 1}; 
    const int nZeros{0};
    USignal::Vector<std::complex<double>> poles(nPoles);
    USignal::Vector<std::complex<double>> zeros(nZeros);
    // Ripple factor
    const double rippleTemp = std::pow(10.0, 0.1*ripple);
#ifndef NDEBUG
    assert(rippleDecibels > 1.0);
#endif
    const double eps = std::sqrt(rippleTemp - 1.0);
    const double mu = 1.0/(static_cast<double> (nPoles))*std::asinh(1.0/eps);
    // Arrange poles in an ellipse on the left half of the S-plane
    constexpr std::complex<double> one(1, 0);
    std::complex<double> product = one; // Initialize product
    double twoni = 1.0/(2.0*static_cast<double> (nPoles));
    for (int i = 0; i < nPoles; i++)
    {
        double xm = static_cast<double> (-nPoles + 1 + 2*i);
        double theta = (M_PI*xm)*twoni;
        std::complex<double> arg(mu, theta);
        poles[i] =-std::sinh(arg);
        product = -poles[i]*product;
    }
    double gain = std::real(product); //Take real
    if (nPoles%2 == 0){gain = gain/std::sqrt(1.0 + eps*eps);}
    USignal::FilterRepresentations::ZerosPolesGain<double>
        zpk{zeros, poles, gain};
    return zpk;
}

USignal::FilterRepresentations::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::
    AnalogPrototype::chebyshevTypeII(const int order, const double ripple)
{
    if (order < 0)
    {   
        throw std::invalid_argument("order = " + std::to_string(order)
                                 + " must be positive");
    }   
    if (ripple <= 0)
    {   
        throw std::invalid_argument("ripple = " + std::to_string(ripple)
                                  + " must be positive");
    }   
    // Set space
    const int nPoles{order + 1};
    const int nZeros = (order + 1)%2 == 1 ? order : order + 1;
    USignal::Vector<std::complex<double>> poles(nPoles);
    USignal::Vector<std::complex<double>> zeros(nZeros);
    // Ripple factor
    const double rippleTemp = std::pow(10.0, 0.1*ripple);
#ifndef NDEBUG
    assert(rippleDecibels > 1.0);
#endif
    const double twoni = 1.0/(2.0*static_cast<double> (order + 1));
    const double eps = 1.0/std::sqrt(rippleTemp - 1.0);
    const double mu = std::asinh(1.0/eps)/static_cast<double> (order + 1);
    // Compute zeros
    std::complex<double> denominator(1, 0);
    // Odd
    if ((order + 1)%2 == 1)
    {
        int j = 0;
        for (int i = 0; i < order/2; ++i)
        {
            auto xm = static_cast<double> (-(order + 1) + 1 + 2*i);
            zeros[j]
                = std::complex<double> (0.0, 1.0/(std::sin(xm*M_PI*twoni)));
            denominator = -zeros[j]*denominator;
            j = j + 1;
        }
        for (int i = 0; i < order/2; ++i)
        {
            auto xm = static_cast<double> (2 + 2*i);
            zeros[j]
                = std::complex<double> (0.0, 1.0/(std::sin(xm*M_PI*twoni)));
            denominator = -zeros[j]*denominator;
            j = j + 1;
        }
    }
    // Even
    else
    {
        for (int i = 0; i < nZeros; ++i)
        {
            auto xm = static_cast<double> (-(order + 1) + 1 + 2*i);
            zeros[i] = std::complex<double> (0.0, 1.0/(std::sin(xm*M_PI*twoni)));
            denominator = -zeros[i]*denominator;
        }
    }
    // Poles around unit circle like butterworth; then warp into cheby II
    const double sinhmu = 0.5*(std::exp(mu) - std::exp(-mu));
    const double coshmu = 0.5*(std::exp(mu) + std::exp(-mu));
    std::complex<double> numerator(1, 0);
    for (int i = 0; i < nPoles; i++)
    {
       //arg = 0.0 + (M_PI*((double)(-n + 1 + 2*i))*twoni)*_Complex_I;
       double temp = M_PI*static_cast<double> (-nPoles + 1 + 2*i)*twoni;
       std::complex<double> arg(0.0, temp);
       std::complex<double> polesi =-std::exp(arg);
       // Warp into chebyshev II
       poles[i] = std::complex<double> (sinhmu*std::real(polesi),
                                        coshmu*std::imag(polesi));
       poles[i] = std::complex<double> (1, 0)/poles[i]; //p = 1/p
       numerator = -poles[i]*numerator;
    }
    double gain = std::real(numerator/denominator);
    USignal::FilterRepresentations::ZerosPolesGain<double>
        zpk{zeros, poles, gain};
    return zpk;
}

USignal::FilterRepresentations::ZerosPolesGain<double>
USignal::FilterDesign::InfiniteImpulseResponse::
    AnalogPrototype::bessel(const int order)
{
    if (order < 0)
    {
        throw std::invalid_argument("order = " + std::to_string(order)
                                 + " must be positive");
    }
    const int nPoles{order + 1};
    const int nZeros{0};
    USignal::Vector<std::complex<double>> poles(nPoles);
    USignal::Vector<std::complex<double>> zeros(nZeros);
    constexpr double gain{1};
    // Load precomputed poles
    if (order == 0)
    {
        poles[0] = std::complex<double> (-1, 0);
    }
    else if (order == 1)
    {
        poles[0] = std::complex<double> (-0.866025403784438597, 0.500000000000000111);
        poles[1] = std::complex<double> (-0.866025403784438597, -0.500000000000000111);
    }
    else if (order == 2)
    {
        poles[0] = std::complex<double> (-0.745640385848076570, 0.711366624972835093);
        poles[1] = std::complex<double> (-0.941600026533206735, -0.000000000000000000);
        poles[2] = std::complex<double> (-0.745640385848076570, -0.711366624972835093);
    }
    else if (order == 3)
    {
        poles[0] = std::complex<double> (-0.657211171671882588, 0.830161435004873161);
        poles[1] = std::complex<double> (-0.904758796788244779, 0.270918733003874646);
        poles[2] = std::complex<double> (-0.904758796788244779, -0.270918733003874646);
        poles[3] = std::complex<double> (-0.657211171671882588, -0.830161435004873161);
    }
    else if (order == 4)
    {
        poles[0] = std::complex<double> (-0.590575944611919090, 0.907206756457454855);
        poles[1] = std::complex<double> (-0.851553619368839554, 0.442717463944332756);
        poles[2] = std::complex<double> (-0.926442077387759966, -0.000000000000000000);
        poles[3] = std::complex<double> (-0.851553619368839554, -0.442717463944332756);
        poles[4] = std::complex<double> (-0.590575944611919090, -0.907206756457454855);
    }
    else if (order == 5)
    {
        poles[0] = std::complex<double> (-0.538552681669310918, 0.961687688195427604);
        poles[1] = std::complex<double> (-0.799654185832828657, 0.562171734693731828);
        poles[2] = std::complex<double> (-0.909390683047227033, 0.185696439679304687);
        poles[3] = std::complex<double> (-0.909390683047227033, -0.185696439679304687);
        poles[4] = std::complex<double> (-0.799654185832828657, -0.562171734693731828);
        poles[5] = std::complex<double> (-0.538552681669310918, -0.961687688195427604);
    }
    else if (order == 6)
    {
        poles[0] = std::complex<double> (-0.496691725667232020, 1.002508508454420522);
        poles[1] = std::complex<double> (-0.752735543409321695, 0.650469630552255262);
        poles[2] = std::complex<double> (-0.880002934152337879, 0.321665276230774067);
        poles[3] = std::complex<double> (-0.919487155649029053, -0.000000000000000000);
        poles[4] = std::complex<double> (-0.880002934152337879, -0.321665276230774067);
        poles[5] = std::complex<double> (-0.752735543409321695, -0.650469630552255262);
        poles[6] = std::complex<double> (-0.496691725667232020, -1.002508508454420522);
    }
    else if (order == 7)
    {
        poles[0] = std::complex<double> (-0.462174041253212320, 1.034388681126901188);
        poles[1] = std::complex<double> (-0.711138180848539747, 0.718651731410840156);
        poles[2] = std::complex<double> (-0.847325080235933448, 0.425901753827293450);
        poles[3] = std::complex<double> (-0.909683154665291149, 0.141243797667142290);
        poles[4] = std::complex<double> (-0.909683154665291149, -0.141243797667142290);
        poles[5] = std::complex<double> (-0.847325080235933448, -0.425901753827293450);
        poles[6] = std::complex<double> (-0.711138180848539747, -0.718651731410840156);
        poles[7] = std::complex<double> (-0.462174041253212320, -1.034388681126901188);
    }
    else if (order == 8)
    {
        poles[0] = std::complex<double> (-0.433141556155362317, 1.060073670135930124);
        poles[1] = std::complex<double> (-0.674362268685476551, 0.773054621269118503);
        poles[2] = std::complex<double> (-0.814802111226901271, 0.508581568963149988);
        poles[3] = std::complex<double> (-0.891121701707975888, 0.252658093458216437);
        poles[4] = std::complex<double> (-0.915495779749903704, -0.000000000000000000);
        poles[5] = std::complex<double> (-0.891121701707975888, -0.252658093458216437);
        poles[6] = std::complex<double> (-0.814802111226901271, -0.508581568963149988);
        poles[7] = std::complex<double> (-0.674362268685476551, -0.773054621269118503);
        poles[8] = std::complex<double> (-0.433141556155362317, -1.060073670135930124);
    }
    else if (order == 9)
    {
        poles[0] = std::complex<double> (-0.408322073286886023, 1.081274842819124560);
        poles[1] = std::complex<double> (-0.641751386698832027, 0.817583616719102069);
        poles[2] = std::complex<double> (-0.783769441310144366, 0.575914753849994909);
        poles[3] = std::complex<double> (-0.868845964128476256, 0.343000823376630959);
        poles[4] = std::complex<double> (-0.909134732090050468, 0.113958313733551142);
        poles[5] = std::complex<double> (-0.909134732090050468, -0.113958313733551142);
        poles[6] = std::complex<double> (-0.868845964128476256, -0.343000823376630959);
        poles[7] = std::complex<double> (-0.783769441310144366, -0.575914753849994909);
        poles[8] = std::complex<double> (-0.641751386698832027, -0.817583616719102069);
        poles[9] = std::complex<double> (-0.408322073286886023, -1.081274842819124560);
    }
    else if (order == 10)
    {
        poles[0] = std::complex<double> (-0.386814951005509500, 1.099117466763121609);
        poles[1] = std::complex<double> (-0.612687155491519531, 0.854781389331476849);
        poles[2] = std::complex<double> (-0.754693893472230815, 0.631915005072185121);
        poles[3] = std::complex<double> (-0.845304401471296374, 0.417869691780124897);
        poles[4] = std::complex<double> (-0.896365670572116913, 0.208048037507103240);
        poles[5] = std::complex<double> (-0.912906724451898466, -0.000000000000000000);
        poles[6] = std::complex<double> (-0.896365670572116913, -0.208048037507103240);
        poles[7] = std::complex<double> (-0.845304401471296374, -0.417869691780124897);
        poles[8] = std::complex<double> (-0.754693893472230815, -0.631915005072185121);
        poles[9] = std::complex<double> (-0.612687155491519531, -0.854781389331476849);
        poles[10] = std::complex<double> (-0.386814951005509500, -1.099117466763121609);
    }
    else if (order == 11)
    {
        poles[0] = std::complex<double> (-0.367964008552631394, 1.114373575641546266);
        poles[1] = std::complex<double> (-0.586636932186147542, 0.886377275132072651);
        poles[2] = std::complex<double> (-0.727668161539516190, 0.679296117876469485);
        poles[3] = std::complex<double> (-0.821729693993907606, 0.481021211510067659);
        poles[4] = std::complex<double> (-0.880253434201683227, 0.287177950352422773);
        poles[5] = std::complex<double> (-0.908447823414068600, 0.095506365213450406);
        poles[6] = std::complex<double> (-0.908447823414068600, -0.095506365213450406);
        poles[7] = std::complex<double> (-0.880253434201683227, -0.287177950352422773);
        poles[8] = std::complex<double> (-0.821729693993907606, -0.481021211510067659);
        poles[9] = std::complex<double> (-0.727668161539516190, -0.679296117876469485);
        poles[10] = std::complex<double> (-0.586636932186147542, -0.886377275132072651);
        poles[11] = std::complex<double> (-0.367964008552631394, -1.114373575641546266);
    }
    else if (order == 12)
    {
        poles[0] = std::complex<double> (-0.351279232338981673, 1.127591548317704806);
        poles[1] = std::complex<double> (-0.563155984243019159, 0.913590033832510251);
        poles[2] = std::complex<double> (-0.702623467572127569, 0.719961189017130354);
        poles[3] = std::complex<double> (-0.798746069247097235, 0.535075212069680117);
        poles[4] = std::complex<double> (-0.862509419826055113, 0.354741373117299086);
        poles[5] = std::complex<double> (-0.899131466547519631, 0.176834295616104364);
        poles[6] = std::complex<double> (-0.911091466598418220, -0.000000000000000000);
        poles[7] = std::complex<double> (-0.899131466547519631, -0.176834295616104364);
        poles[8] = std::complex<double> (-0.862509419826055113, -0.354741373117299086);
        poles[9] = std::complex<double> (-0.798746069247097235, -0.535075212069680117);
        poles[10] = std::complex<double> (-0.702623467572127569, -0.719961189017130354);
        poles[11] = std::complex<double> (-0.563155984243019159, -0.913590033832510251);
        poles[12] = std::complex<double> (-0.351279232338981673, -1.127591548317704806);
    }
    else if (order == 13)
    {
        poles[0] = std::complex<double> (-0.336386822490203019, 1.139172297839859516);
        poles[1] = std::complex<double> (-0.541876677511230587, 0.937304368351692618);
        poles[2] = std::complex<double> (-0.679425642511922501, 0.755285730504203112);
        poles[3] = std::complex<double> (-0.776659138706362606, 0.581917067737760862);
        poles[4] = std::complex<double> (-0.844119916090985134, 0.413165382510269352);
        poles[5] = std::complex<double> (-0.886950667491644640, 0.247007917876533367);
        poles[6] = std::complex<double> (-0.907793213839649060, 0.082196399419401531);
        poles[7] = std::complex<double> (-0.907793213839649060, -0.082196399419401531);
        poles[8] = std::complex<double> (-0.886950667491644640, -0.247007917876533367);
        poles[9] = std::complex<double> (-0.844119916090985134, -0.413165382510269352);
        poles[10] = std::complex<double> (-0.776659138706362606, -0.581917067737760862);
        poles[11] = std::complex<double> (-0.679425642511922501, -0.755285730504203112);
        poles[12] = std::complex<double> (-0.541876677511230587, -0.937304368351692618);
        poles[13] = std::complex<double> (-0.336386822490203019, -1.139172297839859516);
    }
    else if (order == 14)
    {
        poles[0] = std::complex<double> (-0.322996305976644360, 1.149416154583629890);
        poles[1] = std::complex<double> (-0.522495406965833076, 0.958178726109252810);
        poles[2] = std::complex<double> (-0.657919659311100413, 0.786289550372251900);
        poles[3] = std::complex<double> (-0.755602716897072257, 0.622939635875826347);
        poles[4] = std::complex<double> (-0.825663145258714870, 0.464234875273432657);
        poles[5] = std::complex<double> (-0.873126462083498422, 0.308235247056426687);
        poles[6] = std::complex<double> (-0.900698169417697758, 0.153768119727843905);
        poles[7] = std::complex<double> (-0.909748236384906206, -0.000000000000000000);
        poles[8] = std::complex<double> (-0.900698169417697758, -0.153768119727843905);
        poles[9] = std::complex<double> (-0.873126462083498422, -0.308235247056426687);
        poles[10] = std::complex<double> (-0.825663145258714870, -0.464234875273432657);
        poles[11] = std::complex<double> (-0.755602716897072257, -0.622939635875826347);
        poles[12] = std::complex<double> (-0.657919659311100413, -0.786289550372251900);
        poles[13] = std::complex<double> (-0.522495406965833076, -0.958178726109252810);
        poles[14] = std::complex<double> (-0.322996305976644360, -1.149416154583629890);
    }
    else if (order == 15)
    {
        poles[0] = std::complex<double> (-0.310878275564538509, 1.158552841199330441);
        poles[1] = std::complex<double> (-0.504760644442476036, 0.976713747779908603);
        poles[2] = std::complex<double> (-0.637950251403906599, 0.813745353710876196);
        poles[3] = std::complex<double> (-0.735616630471311872, 0.659195087786039524);
        poles[4] = std::complex<double> (-0.807479029323600384, 0.509293375117179870);
        poles[5] = std::complex<double> (-0.858426423152132245, 0.362169727180206291);
        poles[6] = std::complex<double> (-0.891172307032364275, 0.216708965990057567);
        poles[7] = std::complex<double> (-0.907209959508700203, 0.072142113041117312);
        poles[8] = std::complex<double> (-0.907209959508700203, -0.072142113041117312);
        poles[9] = std::complex<double> (-0.891172307032364275, -0.216708965990057567);
        poles[10] = std::complex<double> (-0.858426423152132245, -0.362169727180206291);
        poles[11] = std::complex<double> (-0.807479029323600384, -0.509293375117179870);
        poles[12] = std::complex<double> (-0.735616630471311872, -0.659195087786039524);
        poles[13] = std::complex<double> (-0.637950251403906599, -0.813745353710876196);
        poles[14] = std::complex<double> (-0.504760644442476036, -0.976713747779908603);
        poles[15] = std::complex<double> (-0.310878275564538509, -1.158552841199330441);
    }
    else if (order == 16)
    {
        poles[0] = std::complex<double> (-0.299848945999007410, 1.166761272925667559);
        poles[1] = std::complex<double> (-0.488462933767270679, 0.993297195631678176);
        poles[2] = std::complex<double> (-0.619371071734213685, 0.838249725282698588);
        poles[3] = std::complex<double> (-0.716689384237234606, 0.691493628639360702);
        poles[4] = std::complex<double> (-0.789764414779970059, 0.549372440528108519);
        poles[5] = std::complex<double> (-0.843341449583612790, 0.410075928291002145);
        poles[6] = std::complex<double> (-0.880110070443862469, 0.272534715647880288);
        poles[7] = std::complex<double> (-0.901627385078727861, 0.136026799517302371);
        poles[8] = std::complex<double> (-0.908714116133638838, -0.000000000000000000);
        poles[9] = std::complex<double> (-0.901627385078727861, -0.136026799517302371);
        poles[10] = std::complex<double> (-0.880110070443862469, -0.272534715647880288);
        poles[11] = std::complex<double> (-0.843341449583612790, -0.410075928291002145);
        poles[12] = std::complex<double> (-0.789764414779970059, -0.549372440528108519);
        poles[13] = std::complex<double> (-0.716689384237234606, -0.691493628639360702);
        poles[14] = std::complex<double> (-0.619371071734213685, -0.838249725282698588);
        poles[15] = std::complex<double> (-0.488462933767270679, -0.993297195631678176);
        poles[16] = std::complex<double> (-0.299848945999007410, -1.166761272925667559);
    }
    else if (order == 17)
    {
        poles[0] = std::complex<double> (-0.289759202988048470, 1.174183010600058363);
        poles[1] = std::complex<double> (-0.473426806991615268, 1.008234300314800880);
        poles[2] = std::complex<double> (-0.602048266809064425, 0.860270896189366363);
        poles[3] = std::complex<double> (-0.698782144500527003, 0.720469650972662801);
        poles[4] = std::complex<double> (-0.772628503073955697, 0.585277816208663926);
        poles[5] = std::complex<double> (-0.828188501624283147, 0.452938569781591360);
        poles[6] = std::complex<double> (-0.868109550362883176, 0.322420492516325763);
        poles[7] = std::complex<double> (-0.893976427813245600, 0.193037464089475780);
        poles[8] = std::complex<double> (-0.906700432416277624, 0.064279241063930681);
        poles[9] = std::complex<double> (-0.906700432416277624, -0.064279241063930681);
        poles[10] = std::complex<double> (-0.893976427813245600, -0.193037464089475780);
        poles[11] = std::complex<double> (-0.868109550362883176, -0.322420492516325763);
        poles[12] = std::complex<double> (-0.828188501624283147, -0.452938569781591360);
        poles[13] = std::complex<double> (-0.772628503073955697, -0.585277816208663926);
        poles[14] = std::complex<double> (-0.698782144500527003, -0.720469650972662801);
        poles[15] = std::complex<double> (-0.602048266809064425, -0.860270896189366363);
        poles[16] = std::complex<double> (-0.473426806991615268, -1.008234300314800880);
        poles[17] = std::complex<double> (-0.289759202988048470, -1.174183010600058363);
    }
    else if (order == 18)
    {
        poles[0] = std::complex<double> (-0.280486685143936099, 1.180931628453290472);
        poles[1] = std::complex<double> (-0.459504344973098222, 1.021768776912670651);
        poles[2] = std::complex<double> (-0.585861332121782818, 0.880181713101456431);
        poles[3] = std::complex<double> (-0.681842441291243939, 0.746627235794776078);
        poles[4] = std::complex<double> (-0.756126097154162680, 0.617648391797017471);
        poles[5] = std::complex<double> (-0.813172555157820054, 0.491536503556245952);
        poles[6] = std::complex<double> (-0.855576876561842226, 0.367292589639987177);
        poles[7] = std::complex<double> (-0.884929058503438282, 0.244259075754981708);
        poles[8] = std::complex<double> (-0.902193763939065585, 0.121956838187202599);
        poles[9] = std::complex<double> (-0.907893421789939925, -0.000000000000000000);
        poles[10] = std::complex<double> (-0.902193763939065585, -0.121956838187202599);
        poles[11] = std::complex<double> (-0.884929058503438282, -0.244259075754981708);
        poles[12] = std::complex<double> (-0.855576876561842226, -0.367292589639987177);
        poles[13] = std::complex<double> (-0.813172555157820054, -0.491536503556245952);
        poles[14] = std::complex<double> (-0.756126097154162680, -0.617648391797017471);
        poles[15] = std::complex<double> (-0.681842441291243939, -0.746627235794776078);
        poles[16] = std::complex<double> (-0.585861332121782818, -0.880181713101456431);
        poles[17] = std::complex<double> (-0.459504344973098222, -1.021768776912670651);
        poles[18] = std::complex<double> (-0.280486685143936099, -1.180931628453290472);
    }
    else if (order == 19)
    {
        poles[0] = std::complex<double> (-0.271929958025165341, 1.187099379810885980);
        poles[1] = std::complex<double> (-0.446570069820514837, 1.034097702560842214);
        poles[2] = std::complex<double> (-0.570702680691571596, 0.898282906646825419);
        poles[3] = std::complex<double> (-0.665812054482992854, 0.770372170110075749);
        poles[4] = std::complex<double> (-0.740278030964676481, 0.646997523760522442);
        poles[5] = std::complex<double> (-0.798425119129060223, 0.526494238881712984);
        poles[6] = std::complex<double> (-0.842790747995666445, 0.407891732629193093);
        poles[7] = std::complex<double> (-0.874956031667333489, 0.290555929656790890);
        poles[8] = std::complex<double> (-0.895915094192576644, 0.174031717591870444);
        poles[9] = std::complex<double> (-0.906257011557676795, 0.057961780277849498);
        poles[10] = std::complex<double> (-0.906257011557676795, -0.057961780277849498);
        poles[11] = std::complex<double> (-0.895915094192576644, -0.174031717591870444);
        poles[12] = std::complex<double> (-0.874956031667333489, -0.290555929656790890);
        poles[13] = std::complex<double> (-0.842790747995666445, -0.407891732629193093);
        poles[14] = std::complex<double> (-0.798425119129060223, -0.526494238881712984);
        poles[15] = std::complex<double> (-0.740278030964676481, -0.646997523760522442);
        poles[16] = std::complex<double> (-0.665812054482992854, -0.770372170110075749);
        poles[17] = std::complex<double> (-0.570702680691571596, -0.898282906646825419);
        poles[18] = std::complex<double> (-0.446570069820514837, -1.034097702560842214);
        poles[19] = std::complex<double> (-0.271929958025165341, -1.187099379810885980);
    }
    else if (order == 20)
    {
        poles[0] = std::complex<double> (-0.264004159583402453, 1.192762031948052304);
        poles[1] = std::complex<double> (-0.434516890681526713, 1.045382255856986076);
        poles[2] = std::complex<double> (-0.556476648891856773, 0.914819840584672805);
        poles[3] = std::complex<double> (-0.650631537860946740, 0.792034934262949686);
        poles[4] = std::complex<double> (-0.725083968710661164, 0.673742606302438318);
        poles[5] = std::complex<double> (-0.784028798040834696, 0.558318634802285718);
        poles[6] = std::complex<double> (-0.829943547067444665, 0.444817773940795691);
        poles[7] = std::complex<double> (-0.864391581364320261, 0.332625851252218663);
        poles[8] = std::complex<double> (-0.888380810666444920, 0.221306921508435006);
        poles[9] = std::complex<double> (-0.902542807319269391, 0.110525257278985642);
        poles[10] = std::complex<double> (-0.907226265314296287, -0.000000000000000000);
        poles[11] = std::complex<double> (-0.902542807319269391, -0.110525257278985642);
        poles[12] = std::complex<double> (-0.888380810666444920, -0.221306921508435006);
        poles[13] = std::complex<double> (-0.864391581364320261, -0.332625851252218663);
        poles[14] = std::complex<double> (-0.829943547067444665, -0.444817773940795691);
        poles[15] = std::complex<double> (-0.784028798040834696, -0.558318634802285718);
        poles[16] = std::complex<double> (-0.725083968710661164, -0.673742606302438318);
        poles[17] = std::complex<double> (-0.650631537860946740, -0.792034934262949686);
        poles[18] = std::complex<double> (-0.556476648891856773, -0.914819840584672805);
        poles[19] = std::complex<double> (-0.434516890681526713, -1.045382255856986076);
        poles[20] = std::complex<double> (-0.264004159583402453, -1.192762031948052304);
    }
    else if (order == 21)
    {
        poles[0] = std::complex<double> (-0.256637698793931945, 1.197982433555213166);
        poles[1] = std::complex<double> (-0.423252874564262804, 1.055755605227546079);
        poles[2] = std::complex<double> (-0.543098305630630551, 0.929994782443987700);
        poles[3] = std::complex<double> (-0.636242768326783170, 0.811887504024634832);
        poles[4] = std::complex<double> (-0.710530545641879008, 0.698226626592452493);
        poles[5] = std::complex<double> (-0.770033293055681578, 0.587425542635115150);
        poles[6] = std::complex<double> (-0.817168208846272059, 0.478561949220278171);
        poles[7] = std::complex<double> (-0.853475403685168943, 0.371038931948232065);
        poles[8] = std::complex<double> (-0.879966145564017421, 0.264436303920153493);
        poles[9] = std::complex<double> (-0.897298313815353188, 0.158435191228986527);
        poles[10] = std::complex<double> (-0.905870226993087058, 0.052774908289999027);
        poles[11] = std::complex<double> (-0.905870226993087058, -0.052774908289999027);
        poles[12] = std::complex<double> (-0.897298313815353188, -0.158435191228986527);
        poles[13] = std::complex<double> (-0.879966145564017421, -0.264436303920153493);
        poles[14] = std::complex<double> (-0.853475403685168943, -0.371038931948232065);
        poles[15] = std::complex<double> (-0.817168208846272059, -0.478561949220278171);
        poles[16] = std::complex<double> (-0.770033293055681578, -0.587425542635115150);
        poles[17] = std::complex<double> (-0.710530545641879008, -0.698226626592452493);
        poles[18] = std::complex<double> (-0.636242768326783170, -0.811887504024634832);
        poles[19] = std::complex<double> (-0.543098305630630551, -0.929994782443987700);
        poles[20] = std::complex<double> (-0.423252874564262804, -1.055755605227546079);
        poles[21] = std::complex<double> (-0.256637698793931945, -1.197982433555213166);
    }
    else if (order == 22)
    {
        poles[0] = std::complex<double> (-0.249769720220895552, 1.202813187870697575);
        poles[1] = std::complex<double> (-0.412698661751014773, 1.065328794475513430);
        poles[2] = std::complex<double> (-0.530492246381019550, 0.943976036401830587);
        poles[3] = std::complex<double> (-0.622590322877134117, 0.830155830281298024);
        poles[4] = std::complex<double> (-0.696596603391270941, 0.720734137475304903);
        poles[5] = std::complex<double> (-0.756466014682988352, 0.614159485947603390);
        poles[6] = std::complex<double> (-0.804556164205317836, 0.509530591222725926);
        poles[7] = std::complex<double> (-0.842380594802112692, 0.406265794823760240);
        poles[8] = std::complex<double> (-0.870946939558741473, 0.303958199395004125);
        poles[9] = std::complex<double> (-0.890928324247125425, 0.202302469938122453);
        poles[10] = std::complex<double> (-0.902756497991250795, 0.101053433531404516);
        poles[11] = std::complex<double> (-0.906673247632499124, -0.000000000000000000);
        poles[12] = std::complex<double> (-0.902756497991250795, -0.101053433531404516);
        poles[13] = std::complex<double> (-0.890928324247125425, -0.202302469938122453);
        poles[14] = std::complex<double> (-0.870946939558741473, -0.303958199395004125);
        poles[15] = std::complex<double> (-0.842380594802112692, -0.406265794823760240);
        poles[16] = std::complex<double> (-0.804556164205317836, -0.509530591222725926);
        poles[17] = std::complex<double> (-0.756466014682988352, -0.614159485947603390);
        poles[18] = std::complex<double> (-0.696596603391270941, -0.720734137475304903);
        poles[19] = std::complex<double> (-0.622590322877134117, -0.830155830281298024);
        poles[20] = std::complex<double> (-0.530492246381019550, -0.943976036401830587);
        poles[21] = std::complex<double> (-0.412698661751014773, -1.065328794475513430);
        poles[22] = std::complex<double> (-0.249769720220895552, -1.202813187870697575);
    }
    else if (order == 23)
    {
        poles[0] = std::complex<double> (-0.243348133752487705, 1.207298683731972577);
        poles[1] = std::complex<double> (-0.402785385519751682, 1.074195196518674900);
        poles[2] = std::complex<double> (-0.518591457482031948, 0.956904838525905688);
        poles[3] = std::complex<double> (-0.609622156737833043, 0.847029243307719670);
        poles[4] = std::complex<double> (-0.683256580353652110, 0.741503269509164897);
        poles[5] = std::complex<double> (-0.743339228508853256, 0.638808421622256928);
        poles[6] = std::complex<double> (-0.792169546234348654, 0.538062849096801576);
        poles[7] = std::complex<double> (-0.831232646681324239, 0.438698593359730660);
        poles[8] = std::complex<double> (-0.861527830401635608, 0.340320211261862571);
        poles[9] = std::complex<double> (-0.883735803455570679, 0.242633523440138416);
        poles[10] = std::complex<double> (-0.898310510439786936, 0.145405613387360994);
        poles[11] = std::complex<double> (-0.905531236337277279, 0.048440066540478693);
        poles[12] = std::complex<double> (-0.905531236337277279, -0.048440066540478693);
        poles[13] = std::complex<double> (-0.898310510439786936, -0.145405613387360994);
        poles[14] = std::complex<double> (-0.883735803455570679, -0.242633523440138416);
        poles[15] = std::complex<double> (-0.861527830401635608, -0.340320211261862571);
        poles[16] = std::complex<double> (-0.831232646681324239, -0.438698593359730660);
        poles[17] = std::complex<double> (-0.792169546234348654, -0.538062849096801576);
        poles[18] = std::complex<double> (-0.743339228508853256, -0.638808421622256928);
        poles[19] = std::complex<double> (-0.683256580353652110, -0.741503269509164897);
        poles[20] = std::complex<double> (-0.609622156737833043, -0.847029243307719670);
        poles[21] = std::complex<double> (-0.518591457482031948, -0.956904838525905688);
        poles[22] = std::complex<double> (-0.402785385519751682, -1.074195196518674900);
        poles[23] = std::complex<double> (-0.243348133752487705, -1.207298683731972577);
    }
    else if (order == 24)
    {
        poles[0] = std::complex<double> (-0.237328066932203097, 1.211476658382565796);
        poles[1] = std::complex<double> (-0.393452987819108091, 1.082433927173832133);
        poles[2] = std::complex<double> (-0.507336286107846868, 0.968900630534487051);
        poles[3] = std::complex<double> (-0.597289866133556280, 0.862667633038803228);
        poles[4] = std::complex<double> (-0.670482712802955794, 0.760734885816783835);
        poles[5] = std::complex<double> (-0.730654927184996805, 0.661614964735775080);
        poles[6] = std::complex<double> (-0.780049627818650171, 0.564444121034971436);
        poles[7] = std::complex<double> (-0.820122604393687782, 0.468666857465696585);
        poles[8] = std::complex<double> (-0.851861688655402349, 0.373897787590759645);
        poles[9] = std::complex<double> (-0.875949798967786464, 0.279852132177141111);
        poles[10] = std::complex<double> (-0.892855145988355359, 0.186306896980430181);
        poles[11] = std::complex<double> (-0.902883339022802489, 0.093077131185103010);
        poles[12] = std::complex<double> (-0.906207387181171109, -0.000000000000000000);
        poles[13] = std::complex<double> (-0.902883339022802489, -0.093077131185103010);
        poles[14] = std::complex<double> (-0.892855145988355359, -0.186306896980430181);
        poles[15] = std::complex<double> (-0.875949798967786464, -0.279852132177141111);
        poles[16] = std::complex<double> (-0.851861688655402349, -0.373897787590759645);
        poles[17] = std::complex<double> (-0.820122604393687782, -0.468666857465696585);
        poles[18] = std::complex<double> (-0.780049627818650171, -0.564444121034971436);
        poles[19] = std::complex<double> (-0.730654927184996805, -0.661614964735775080);
        poles[20] = std::complex<double> (-0.670482712802955794, -0.760734885816783835);
        poles[21] = std::complex<double> (-0.597289866133556280, -0.862667633038803228);
        poles[22] = std::complex<double> (-0.507336286107846868, -0.968900630534487051);
        poles[23] = std::complex<double> (-0.393452987819108091, -1.082433927173832133);
        poles[24] = std::complex<double> (-0.237328066932203097, -1.211476658382565796);
    }
    else
    {
        throw std::invalid_argument("unsupported filter order = "
                                  + std::to_string(order));
    }
    USignal::FilterRepresentations::ZerosPolesGain<double>
        zpk{zeros, poles, gain};
    return zpk;
}
