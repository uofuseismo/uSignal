#include <iostream>
#include <limits>
#include <complex>
#include <cmath>
#include <vector>
#include "uSignal/filterDesign/finiteImpulseResponse/windowBased.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace UFIRDesign = USignal::FilterDesign::FiniteImpulseResponse;

TEMPLATE_TEST_CASE(
    "CoreTest::FilterDesign::FiniteImpulseResponse::WindowBase::lowpass",
    "[TypeName][template]",
    double, float)
{
    // b = fir1(13, 0.6); %default hamming
    USignal::Vector<TestType> ref( {-0.001213632932920, -0.006228195759795,
                                     0.015988039487026,  0.013651616475176,
                                    -0.089746334361735,  0.058133494212324,
                                     0.509415012879924,  0.509415012879924,
                                     0.058133494212324, -0.089746334361735,
                                     0.013651616475176,  0.015988039487026,
                                    -0.006228195759795, -0.001213632932920} );
    constexpr TestType lowPass{0.6};
    constexpr int order{13};
    constexpr auto window{UFIRDesign::WindowBased::Window::Hamming};
    auto firFilter = UFIRDesign::WindowBased::lowPass(order, lowPass, window);
    REQUIRE(firFilter.getOrder() == order);
    auto bs = firFilter.getFilterCoefficients();
    REQUIRE(bs.size() == ref.size());
    for (int i = 0; i < static_cast<int> (bs.size()); ++i)
    {
        if constexpr (std::is_same<TestType, double>::value)
        {
             REQUIRE(std::abs(bs.at(i) - ref.at(i)) < 1.e-12);
        }
        else if constexpr (std::is_same<TestType, float>::value)
        {
             REQUIRE(std::abs(bs.at(i) - ref.at(i)) < 1.e-7);
        }
    }
}

TEMPLATE_TEST_CASE(
    "CoreTest::FilterDesign::FiniteImpulseResponse::WindowBase::highpass",
    "[TypeName][template]",
    double, float)
{
    // fir1(16, 0.45, 'high', hann(16+1));
    USignal::Vector<TestType> ref( {0,                  0.000785256397236,
                                   -0.006281697554181, -0.013886226413618,
                                    0.023373298499376,  0.065319627573863,
                                   -0.041954095979927, -0.302244991506668,
                                    0.549672322171091, -0.302244991506668,
                                   -0.041954095979927,  0.065319627573863,
                                    0.023373298499376, -0.013886226413618,
                                   -0.006281697554181,  0.000785256397236,
                                    0
                                    } );
    constexpr TestType highPass{0.45};
    constexpr int order{16};
    constexpr auto window{UFIRDesign::WindowBased::Window::Hanning};
    auto firFilter = UFIRDesign::WindowBased::highPass(order, highPass, window);
    REQUIRE(firFilter.getOrder() == order);
    auto bs = firFilter.getFilterCoefficients();
    REQUIRE(bs.size() == ref.size());
    for (int i = 0; i < static_cast<int> (bs.size()); ++i)
    {   
        if constexpr (std::is_same<TestType, double>::value)
        {
             REQUIRE(std::abs(bs.at(i) - ref.at(i)) < 1.e-12);
        }
        else if constexpr (std::is_same<TestType, float>::value)
        {
             REQUIRE(std::abs(bs.at(i) - ref.at(i)) < 1.e-7);
        }
    }   
}

TEMPLATE_TEST_CASE(
    "CoreTest::FilterDesign::FiniteImpulseResponse::WindowBase::bandpass",
    "[TypeName][template]",
    double, float)
{
    // fir1(11, [0.2, 0.8], 'bandpass', bartlett(11+1));
    USignal::Vector<TestType> ref( { 0,                 -0.018100893326756,
                                    -0.008171956220958, -0.077570731498438,
                                    -0.240779937850341,  0.415028961385482,
                                     0.415028961385482, -0.240779937850341,
                                    -0.077570731498438, -0.008171956220958,
                                    -0.018100893326756,  0} );
    constexpr std::pair<TestType, TestType> passband{0.2, 0.8};
    constexpr int order{11};
    constexpr auto window{UFIRDesign::WindowBased::Window::Bartlett};
    auto firFilter = UFIRDesign::WindowBased::bandPass(order, passband, window);
    REQUIRE(firFilter.getOrder() == order);
    auto bs = firFilter.getFilterCoefficients();
    REQUIRE(bs.size() == ref.size());
    for (int i = 0; i < static_cast<int> (bs.size()); ++i)
    {   
        if constexpr (std::is_same<TestType, double>::value)
        {
             REQUIRE(std::abs(bs.at(i) - ref.at(i)) < 1.e-12);
        }
        else if constexpr (std::is_same<TestType, float>::value)
        {
             REQUIRE(std::abs(bs.at(i) - ref.at(i)) < 1.e-7);
        }
    }   
}

TEMPLATE_TEST_CASE(
    "CoreTest::FilterDesign::FiniteImpulseResponse::WindowBase::bandstop",
    "[TypeName][template]",
    double, float)
{
    // IPP uses optimal blackman
    // alpha = -0.5/(1 + cos(2*pi/20));
    // fir1(20, [0.15,0.85], 'stop', (alpha+1)/2 - 0.5*cos(2*pi*linspace(0,20,21)/20) - alpha/2*cos(4*pi*linspace(0,20,21)/20))
    USignal::Vector<TestType> ref( { 0.000000000000000,  0,
                                    -0.000380306566042,  0.000000000000000,
                                     0.004359748553493,  0.000000000000000,
                                     0.074832232270797,  0.000000000000000,
                                     0.245754965804324,                  0,
                                     0.350866719874856,                  0,
                                     0.245754965804324,  0.000000000000000,
                                     0.074832232270797,  0.000000000000000,
                                     0.004359748553493,  0.000000000000000,
                                    -0.000380306566042, -0.000000000000000,
                                     0} );
    const std::pair<TestType, TestType> stopband{0.15, 0.85};
    constexpr int order{20};
    constexpr auto window{UFIRDesign::WindowBased::Window::OptimalBlackman};
    auto firFilter = UFIRDesign::WindowBased::bandStop(order, stopband, window);
    REQUIRE(firFilter.getOrder() == order);
    auto bs = firFilter.getFilterCoefficients();
    REQUIRE(bs.size() == ref.size());
    for (int i = 0; i < static_cast<int> (bs.size()); ++i)
    {   
        if constexpr (std::is_same<TestType, double>::value)
        {
             REQUIRE(std::abs(bs.at(i) - ref.at(i)) < 1.e-12);
        }
        else if constexpr (std::is_same<TestType, float>::value)
        {
             REQUIRE(std::abs(bs.at(i) - ref.at(i)) < 1.e-7);
        }
    }   
}

