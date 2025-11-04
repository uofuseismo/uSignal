#include <iostream>
#include <limits>
#include <complex>
#include <cmath>
#include <vector>
#include "uSignal/filterDesign/infiniteImpulseResponse/bilinearTransform.hpp"
#include "uSignal/filterDesign/infiniteImpulseResponse/analogPrototype.hpp"
#include "uSignal/filterDesign/infiniteImpulseResponse/convertBand.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace
{
template<typename T>
bool matched(const std::complex<T> value,
             const USignal::Vector<std::complex<T>> &reference,
             const T tolerance)
{
    for (const auto &referenceValue : reference)
    {
        if (std::abs(value - referenceValue) < tolerance){return true;}
    }
    return false; 
}
}

TEMPLATE_TEST_CASE("CoreTest::FilterDesign::InfiniteImpulseResponse::BilinearTransform",
                   "[TypeName][template]",
                   double, float)
{
    USignal::Vector<std::complex<TestType>> analogZeros(
        std::vector<std::complex<TestType>> {0 + 0i,
                                             0 + 0i,
                                             0 + 0i,
                                             0 + 0i} );
    USignal::Vector<std::complex<TestType>> analogPoles(
         std::vector<std::complex<TestType>> {
            -5.188311655529189 - 44.616532466296775i,
            -15.243050879954183 - 50.63131935118821i,
            -15.243050879954183 + 50.63131935118821i,
            -5.188311655529189 + 44.616532466296775i,
            -9.238513861695127 + 79.44597029196994i,
            -19.586386945718992 + 65.05814486841253i,
            -19.586386945718992 - 65.05814486841253i,
            -9.238513861695127-79.44597029196994i} );
    constexpr TestType analogGain{2019874.9116810758};

    USignal::Vector<std::complex<TestType>> digitalZerosRef(
        std::vector<std::complex<TestType>> { 1 + 0i,
                                              1 + 0i,
                                              1 + 0i,
                                              1 + 0i,
                                             -1 + 0i,
                                             -1 + 0i,
                                             -1 + 0i,
                                             -1 + 0i} );
    USignal::Vector<std::complex<TestType>> digitalPolesRef(
        std::vector<std::complex<TestType>> {
            0.861419077078615-0.40475046563708234i,
            0.7609277750082251-0.4142205574875641i,
            0.7609277750082251+0.4142205574875641i,
            0.861419077078615+0.40475046563708234i,
            0.6708198383910375+0.634395173212654i,
            0.6746102884496724+0.4961465975166494i,
            0.6746102884496724-0.4961465975166494i,
            0.6708198383910375-0.634395173212654i} );
    constexpr TestType digitalGainRef{0.0005705645409457355};

    constexpr TestType samplingFrequency{100};
    USignal::FilterRepresentations::ZerosPolesGain
        analogZPK{analogZeros, analogPoles, analogGain};

    auto digitalZPK 
        = USignal::FilterDesign::InfiniteImpulseResponse::bilinearTransform(
             analogZPK, samplingFrequency);
    auto digitalZeros = digitalZPK.getZeros();
    auto digitalPoles = digitalZPK.getPoles();
    auto digitalGain = digitalZPK.getGain();
    REQUIRE(digitalZerosRef.size() == digitalZeros.size());
    for (int i = 0; i < static_cast<int> (digitalZeros.size()); ++i)
    {
        auto residual = std::abs(digitalZerosRef.at(i) - digitalZeros.at(i));
        CHECK(Catch::Approx(0).margin(
              std::numeric_limits<TestType>::epsilon()*100) == residual);
    }
    REQUIRE(digitalPolesRef.size() == digitalPoles.size());
    for (int i = 0; i < static_cast<int> (digitalPoles.size()); ++i)
    {
        auto residual = std::abs(digitalPolesRef.at(i) - digitalPoles.at(i));
        CHECK(Catch::Approx(0).margin(
              std::numeric_limits<TestType>::epsilon()*100) == residual);
    }
    CHECK(Catch::Approx(digitalGainRef).margin(
             std::numeric_limits<TestType>::epsilon()*100) == digitalGain);
}

TEST_CASE("CoreTest::FilterDesign::InfiniteImpulseResponse::AnalogPrototype::butteworth")
{
    namespace UAnalogPrototype 
        = USignal::FilterDesign::InfiniteImpulseResponse::AnalogPrototype;
    SECTION("Order 1")
    {
        constexpr int order{1};
        auto zpk = UAnalogPrototype::butterworth(order); 
        constexpr double gain{1};
        //REQUIRE(zpk.getOrder() == order);
        REQUIRE(Catch::Approx(gain).margin(
                  std::numeric_limits<double>::epsilon()) == gain);
        REQUIRE(zpk.getZeros().size() == 0);
        REQUIRE(zpk.getPoles().size() == 1);
        constexpr std::complex<double> pole{-1 + 0i};
        REQUIRE(std::abs(pole - zpk.getPoles().at(0)) <
                std::numeric_limits<double>::epsilon());
    }
    SECTION("Order 5")
    {
        constexpr int order{5};
        constexpr double gain{1};
        USignal::Vector<std::complex<double>> polesRef
        {
            std::vector<std::complex<double>>
            {
               -0.3090169943749474  + 0.95105651629515353i,
               -0.80901699437494745 + 0.58778525229247303i,
               -1 + 0i,
               -0.80901699437494745 + -0.58778525229247303i,
               -0.30901699437494751 + -0.95105651629515353i
            }
        };
        auto zpk = UAnalogPrototype::butterworth(order);
        REQUIRE(Catch::Approx(gain).margin(
                  std::numeric_limits<double>::epsilon()) == gain);
        REQUIRE(zpk.getZeros().size() == 0);
        auto poles = zpk.getPoles();
        REQUIRE(poles.size() == polesRef.size());
        for (const auto &p : poles)
        {
            bool matched{false};
            for (const auto &pr : polesRef)
            {
                if (std::abs(p - pr) < std::numeric_limits<double>::epsilon()*10)
                {
                    matched = true;
                    break;
                }
            }
            REQUIRE(matched);
        }
    }
}

TEST_CASE("CoreTest::FilterDesign::InfiniteImpulseResponse::AnalogPrototype::chebyshevI")
{
    namespace UAnalogPrototype 
        = USignal::FilterDesign::InfiniteImpulseResponse::AnalogPrototype;
    SECTION("Order 1")
    {
        const int order{1};
        const double ripple{2.2};
        auto zpk = UAnalogPrototype::chebyshevTypeI(order, ripple);
        REQUIRE(zpk.getZeros().size() == 0);
        REQUIRE(zpk.getPoles().size() == 1);
        constexpr std::complex<double> pole{-1.2313003041963828 + 0i};
        constexpr double gainRef{1.2313003041963828};
        REQUIRE(std::abs(zpk.getGain() - gainRef) <
                std::numeric_limits<double>::epsilon()*100);
        REQUIRE(std::abs(pole - zpk.getPoles().at(0)) <
                std::numeric_limits<double>::epsilon()*100);
    }
    SECTION("Order 6")
    {
        const int order{6};
        const double ripple{0.994};
        auto zpk = UAnalogPrototype::chebyshevTypeI(order, ripple);
        USignal::Vector<std::complex<double>> polesRef
        {
            std::vector<std::complex<double>>
            {
                -0.062314231644038744 + 0.9935274525619241i,
                -0.17024564688613014  + 0.72731257398980598i,
                -0.23255987853016891  + 0.26621487857211812i,
                -0.23255987853016891  +-0.26621487857211795i,
                -0.17024564688613017  +-0.72731257398980587i,
                -0.062314231644038813 +-0.99352745256192398i
            }
        };
        constexpr double gainRef{0.061620501119488615};
        REQUIRE(zpk.getZeros().size() == 0); 
        auto poles = zpk.getPoles();
        REQUIRE(poles.size() == polesRef.size()); 
        REQUIRE(std::abs(zpk.getGain() - gainRef) <
                std::numeric_limits<double>::epsilon()*100);
        REQUIRE(poles.size() == polesRef.size());
        for (const auto &p : poles)
        {
            bool matched{false};
            for (const auto &pr : polesRef)
            {
                if (std::abs(p - pr) < std::numeric_limits<double>::epsilon()*10)
                {
                    matched = true;
                    break;
                }
            }
            REQUIRE(matched);
        }
    }
}

TEST_CASE("CoreTest::FilterDesign::InfiniteImpulseResponse::AnalogPrototype::chebyshevII")
{
    namespace UAnalogPrototype 
        = USignal::FilterDesign::InfiniteImpulseResponse::AnalogPrototype;
    SECTION("Order 1")
    {   
        const int order{1};
        constexpr double ripple{1.1};
        auto zpk = UAnalogPrototype::chebyshevTypeII(order, ripple);
        REQUIRE(zpk.getZeros().size() == 0);
        REQUIRE(zpk.getPoles().size() == 1);
        constexpr std::complex<double> pole{-1.862583192806328 + 0i};
        constexpr double gainRef{1.862583192806328};
        REQUIRE(std::abs(zpk.getGain() - gainRef) <
                std::numeric_limits<double>::epsilon()*100);
        REQUIRE(std::abs(pole - zpk.getPoles().at(0)) <
                std::numeric_limits<double>::epsilon()*100);
    }
    SECTION("Order 3")
    {
        const int order{3};
        constexpr double ripple{1.2};
        auto zpk = UAnalogPrototype::chebyshevTypeII(order, ripple);
        REQUIRE(zpk.getZeros().size() == 2);
        REQUIRE(zpk.getPoles().size() == 3); 
        constexpr double gainRef{5.317805523446352};
        REQUIRE(std::abs(zpk.getGain() - gainRef) <
                std::numeric_limits<double>::epsilon()*100);
        USignal::Vector<std::complex<double>> zerosRef
        {
            std::vector<std::complex<double>>
            {
                -1.1547005383792517i,
                 1.1547005383792517i,
            }
        };
        USignal::Vector<std::complex<double>> polesRef
        {
            std::vector<std::complex<double>>
            {
                -0.11517150279344834-1.1245944747171477i,
                -5.548148529033251+0i,
                -0.11517150279344834+1.1245944747171477i
            }
        };
        auto zeros = zpk.getZeros();
        REQUIRE(zeros.size() == zerosRef.size());
        for (const auto &z : zeros)
        {
            bool matched{false};
            for (const auto &zr : zerosRef)
            {
                if (std::abs(z - zr) < std::numeric_limits<double>::epsilon()*10)
                {
                    matched = true;
                    break;
                }
            }
            REQUIRE(matched);
        }
        auto poles = zpk.getPoles();
        REQUIRE(poles.size() == polesRef.size());
        for (const auto &p : poles)
        {
            bool matched{false};
            for (const auto &pr : polesRef)
            {
                if (std::abs(p - pr) < std::numeric_limits<double>::epsilon()*10)
                {
                    matched = true;
                    break;
                }
            }
            REQUIRE(matched);
        }
    }
    SECTION("Order 6")
    {
        const int order{6};
        constexpr double ripple{1.2};
        auto zpk = UAnalogPrototype::chebyshevTypeII(order, ripple);
        constexpr double gainRef{0.8709635899560811};
        REQUIRE(std::abs(zpk.getGain() - gainRef) <
                std::numeric_limits<double>::epsilon()*100);
        USignal::Vector<std::complex<double>> zerosRef
        {
            std::vector<std::complex<double>>
            { 
                0 - 1.035276180410083i,
                0 - 1.4142135623730951i,
                0 - 3.8637033051562737i,
                0 + 3.8637033051562737i,
                0 + 1.4142135623730951i,
                0 + 1.035276180410083i
            }
        };
        USignal::Vector<std::complex<double>> polesRef
        {
            std::vector<std::complex<double>>
            { 
                -0.024686186266327684 -1.0305393933278832i,
                -0.12492582633346083 -1.3973824335027194i,
                -1.1553327165440157 -3.462761441343769i,
                -1.1553327165440157 +3.462761441343769i,
                -0.12492582633346083 +1.3973824335027194i,
                -0.024686186266327684 +1.0305393933278832i
            }
        };
        auto zeros = zpk.getZeros();
        REQUIRE(zeros.size() == zerosRef.size());
        for (const auto &z : zeros)
        {
            REQUIRE(::matched(z, zerosRef,
                              std::numeric_limits<double>::epsilon()*10));
        }

        auto poles = zpk.getPoles();
        REQUIRE(poles.size() == polesRef.size());
        for (const auto &p : poles)
        {
            REQUIRE(::matched(p, polesRef,
                              std::numeric_limits<double>::epsilon()*10));
        }
    }
}

TEST_CASE("Convert lowpass to lowpass")
{
    namespace UIIRDesign = USignal::FilterDesign::InfiniteImpulseResponse;
    const int order{3};
    const double ripple{3};
    auto zpkLowpassAnalogPrototype
        = UIIRDesign::AnalogPrototype::chebyshevTypeI(order, ripple);
    constexpr double cutoffFrequency{2};
    constexpr double cutoffFrequencyRadSeconds{cutoffFrequency*(2*M_PI)};
    auto zpkLowpass
        = UIIRDesign::convertLowpassAnalogPrototypeToLowpass(
            zpkLowpassAnalogPrototype, cutoffFrequencyRadSeconds);
    // [z, p, k] = signal.cheb1ap(3, 3.)
    // [z1, p1, k1] = signal.lp2lp_zpk(z, p, k, 2*2*pi)
    USignal::Vector<std::complex<double>> polesRef
    {   
        std::vector<std::complex<double>>
        {
            -1.876286105033209+11.35766707809726i,
            -3.752572210066417-0i,
            -1.876286105033209-11.35766707809726i
        }
    };
    constexpr double gainReference{497.279802960698};
    auto poles = zpkLowpass.getPoles();
    REQUIRE(std::abs(zpkLowpass.getGain() - gainReference) < 1.e-10);
    REQUIRE(zpkLowpass.getZeros().size() == 0);
    REQUIRE(poles.size() == polesRef.size());  
    for (const auto &p : poles)
    {
        REQUIRE(::matched(p, polesRef,
                          std::numeric_limits<double>::epsilon()*100));
    }
}

TEST_CASE("Convert lowpass to highpass")
{
    namespace UIIRDesign = USignal::FilterDesign::InfiniteImpulseResponse;
    const int order{3};
    const double ripple{3};
    auto zpkLowpassAnalogPrototype
        = UIIRDesign::AnalogPrototype::chebyshevTypeI(order, ripple);
    constexpr double cutoffFrequency{2};
    constexpr double cutoffFrequencyRadSeconds{cutoffFrequency*(2*M_PI)};
    auto zpkHighpass
        = UIIRDesign::convertLowpassAnalogPrototypeToHighpass(
            zpkLowpassAnalogPrototype, cutoffFrequencyRadSeconds);
    // [z, p, k] = signal.cheb1ap(3, 3.)
    // [z1, p1, k1] = signal.lp2hp_zpk(z, p, k, 2*2*pi)
    USignal::Vector<std::complex<double>> zerosRef
    {
        std::vector<std::complex<double>> {0 + 0i, 0 + 0i, 0 + 0i}
    };
    USignal::Vector<std::complex<double>> polesRef
    {   
        std::vector<std::complex<double>>
        {
            -2.235872465862575-13.534340540192852i,
            -42.081447491888454-0i,
            -2.235872465862575+13.534340540192852i
        }
    };  
    constexpr double gainReference{1};
    auto zeros = zpkHighpass.getZeros();
    auto poles = zpkHighpass.getPoles();
    REQUIRE(std::abs(zpkHighpass.getGain() - gainReference) < 1.e-10);
    REQUIRE(zeros.size() == zerosRef.size()); 
    for (const auto &z : zeros)
    {
        REQUIRE(::matched(z, zerosRef,
                          std::numeric_limits<double>::epsilon()*100));
    }
    REQUIRE(poles.size() == polesRef.size());  
    for (const auto &p : poles)
    {   
        REQUIRE(::matched(p, polesRef,
                          std::numeric_limits<double>::epsilon()*100));
    }   
}

TEST_CASE("Convert lowpass to bandpass")
{
    namespace UIIRDesign = USignal::FilterDesign::InfiniteImpulseResponse;
    const int order{4};
    auto zpkLowpassAnalogPrototype
        = UIIRDesign::AnalogPrototype::butterworth(order);
    constexpr double lowCornerRadSeconds{4*(2*M_PI)};
    constexpr double bandwidthRadSeconds{10*(2*M_PI)};
    constexpr double highCornerRadSeconds{lowCornerRadSeconds + bandwidthRadSeconds};
    constexpr std::pair<double, double>
        cutoffFrequenciesRadSeconds{lowCornerRadSeconds, highCornerRadSeconds};
    auto zpkBandpass
        = UIIRDesign::convertLowpassAnalogPrototypeToBandpass(
            zpkLowpassAnalogPrototype, cutoffFrequenciesRadSeconds);
    USignal::Vector<std::complex<double>> zerosRef
    {   
        std::vector<std::complex<double>> {0 + 0i, 0 + 0i, 0 + 0i, 0 + 0i} 
    };
    USignal::Vector<std::complex<double>> polesRef
    {
        std::vector<std::complex<double>>
        {
            -2.747844608990089-8.59936563405789i,
            -9.43813420914492-5.793234681905291i,
            -9.43813420914492+5.793234681905291i,
            -2.747844608990089+8.59936563405789i,
           -21.296864586383766+66.6484286768465i,
           -48.6109288336437+29.837943877279145i,
           -48.6109288336437-29.837943877279145i,
           -21.296864586383766-66.6484286768465i,
        }
    };
    constexpr double gainReference{15585454.565440388}; 
    REQUIRE(std::abs(zpkBandpass.getGain() - gainReference) < 1.e-5);
    auto zeros = zpkBandpass.getZeros();
    auto poles = zpkBandpass.getPoles();
    REQUIRE(zerosRef.size() == zeros.size());
    for (const auto &z : zeros)
    {
        REQUIRE(::matched(z, zerosRef,
                          std::numeric_limits<double>::epsilon()*100));
    }
    REQUIRE(polesRef.size() == poles.size());
    for (const auto &p : poles)
    {   
        REQUIRE(::matched(p, polesRef,
                          std::numeric_limits<double>::epsilon()*200));
    }    
}

TEST_CASE("Convert lowpass to bandstop")
{
    namespace UIIRDesign = USignal::FilterDesign::InfiniteImpulseResponse;
    const int order{4};
    auto zpkLowpassAnalogPrototype
        = UIIRDesign::AnalogPrototype::butterworth(order);
    constexpr double lowCornerRadSeconds{4*(2*M_PI)};
    constexpr double bandwidthRadSeconds{10*(2*M_PI)};
    constexpr double highCornerRadSeconds{lowCornerRadSeconds + bandwidthRadSeconds};
    constexpr std::pair<double, double>
        cutoffFrequenciesRadSeconds{lowCornerRadSeconds, highCornerRadSeconds};
    auto zpkBandstop
        = UIIRDesign::convertLowpassAnalogPrototypeToBandstop(
            zpkLowpassAnalogPrototype, cutoffFrequenciesRadSeconds);
    USignal::Vector<std::complex<double>> zerosRef
    {
        std::vector<std::complex<double>> 
        {
            0 + 25.132741228718345i,
            0 + 25.132741228718345i,
            0 + 25.132741228718345i,
            0 + 25.132741228718345i,
            0 - 25.132741228718345i,     
            0 - 25.132741228718345i,     
            0 - 25.132741228718345i,     
            0 - 25.132741228718345i
        }
    };
    USignal::Vector<std::complex<double>> polesRef
    {   
        std::vector<std::complex<double>>
        {   
            -2.747844608990089+8.59936563405789i,
            -9.43813420914492+5.79323468190529i,
            -9.43813420914492-5.79323468190529i,
            -2.747844608990089-8.59936563405789i,
            -21.296864586383766-66.6484286768465i,
            -48.6109288336437-29.837943877279145i,
            -48.6109288336437+29.837943877279145i,
            -21.296864586383766+66.6484286768465i
        }
    };
    constexpr double gainReference{1};
    REQUIRE(std::abs(zpkBandstop.getGain() - gainReference) < 1.e-14);
    auto zeros = zpkBandstop.getZeros();
    auto poles = zpkBandstop.getPoles();
    REQUIRE(zerosRef.size() == zeros.size());
    for (const auto &z : zeros)
    {
        REQUIRE(::matched(z, zerosRef,
                          std::numeric_limits<double>::epsilon()*200));
    }   
    REQUIRE(polesRef.size() == poles.size());
    for (const auto &p : poles)
    {
        REQUIRE(::matched(p, polesRef,
                          std::numeric_limits<double>::epsilon()*200));
    }    
}
