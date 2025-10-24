#include <iostream>
#include <limits>
#include <complex>
#include <cmath>
#include <vector>
#include "uSignal/filterDesign/finiteImpulseResponse/windowBased.hpp"
#include "uSignal/filterDesign/finiteImpulseResponse/hilbertTransformer.hpp"
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

TEMPLATE_TEST_CASE(
    "CoreTest::FilterDesign::FiniteImpulseResponse::hilbertTransformer",
    "[TypeName][template]",
    double)
{
    TestType tol{std::numeric_limits<float>::epsilon()*10};
    if constexpr (std::is_same<TestType, double>::value)
    {
        tol = std::numeric_limits<double>::epsilon()*100;
    } 

    namespace UFIRFilterDesign = USignal::FilterDesign::FiniteImpulseResponse;
    SECTION("Order 0")
    {
        constexpr int order{0};
        constexpr double beta{8};
        auto firFilter
            = UFIRFilterDesign::hilbertTransformer<TestType> (order, beta);
        auto taps = firFilter.getFilterCoefficients();
        REQUIRE(taps.size() == 1);
        REQUIRE(std::abs(taps.at(0) - std::complex<TestType> (1, 0)) < tol);
    }
    SECTION("Order 16")
    {
        USignal::Vector<std::complex<TestType>> ref
        ( 
            std::vector<std::complex<TestType>>
                                   {0 + -0.000000000000000i,
                                    0 + -0.002154949073886i,
                                    0 + -0.000000000000000i,
                                    0 + -0.025048021740503i,
                                    0 + -0.000000000000000i,
                                    0 + -0.123110554263836i,
                                    0 + -0.000000000000000i,
                                    0 + -0.600346453738133i,
                                    1 +  0i,
                                    0 +  0.600346453738133i,
                                    0 +  0.000000000000000i,
                                    0 +  0.123110554263836i,
                                    0 +  0.000000000000000i,
                                    0 +  0.025048021740503i,
                                    0 +  0.000000000000000i,
                                    0 +  0.002154949073886i,
                                    0 +  0.000000000000000}
        );
        constexpr int order{16};
        constexpr double beta{8};
        auto firFilter
            = UFIRFilterDesign::hilbertTransformer<TestType> (order, beta);
        auto bs = firFilter.getFilterCoefficients();
        REQUIRE(ref.size() == bs.size());
        for (int i = 0; i < static_cast<int> (ref.size()); ++i)
        {
            auto residual = std::abs(ref.at(i) - bs.at(i));
            CHECK(residual < tol);
        }
    }
    SECTION("Order 19")
    {
        USignal::Vector<std::complex<TestType>> ref 
        (
            std::vector<std::complex<TestType>>
                                   {
                                      -0.000078362668030 + -0.000078362668030i,
                                       0.000687065652501 + -0.000687065652501i,
                                      -0.002495774178985 + -0.002495774178985i,
                                       0.006621258126556 + -0.006621258126556i,
                                      -0.014693594258722 + -0.014693594258722i,
                                       0.029092392815675 + -0.029092392815675i,
                                      -0.053803945082584 + -0.053803945082584i,
                                       0.097839237445465 + -0.097839237445464i,
                                      -0.193197503788858 + -0.193197503788858i,
                                       0.630029225936982 + -0.630029225936981i,
                                       0.630029225936982 +  0.630029225936981i,
                                      -0.193197503788858 +  0.193197503788858i,
                                       0.097839237445465 +  0.097839237445464i,
                                      -0.053803945082584 +  0.053803945082584i,
                                       0.029092392815675 +  0.029092392815675i,
                                      -0.014693594258722 +  0.014693594258722i,
                                       0.006621258126556 +  0.006621258126556i,
                                      -0.002495774178985 +  0.002495774178985i,
                                       0.000687065652501 +  0.000687065652501i,
                                      -0.000078362668030 +  0.000078362668030i
                                   }
        );
        constexpr int order{19};
        constexpr double beta{8};
        auto firFilter
            = UFIRFilterDesign::hilbertTransformer<TestType> (order, beta);
        auto bs = firFilter.getFilterCoefficients();
        REQUIRE(ref.size() == bs.size());
        for (int i = 0; i < static_cast<int> (ref.size()); ++i)
        {
            auto residual = std::abs(ref.at(i) - bs.at(i));
            CHECK(residual < tol);
        }

    }
    
}

