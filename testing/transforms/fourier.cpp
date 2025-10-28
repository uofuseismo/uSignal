#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include "uSignal/vector.hpp"
#include "uSignal/transforms/fourier/forward.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;
namespace UFT = USignal::Transforms::Fourier;

TEST_CASE("CoreTest::Transforms::Fourier::ForwardOptions")
{
    SECTION("DFT")
    {
        UFT::ForwardOptions options{
           UFT::ForwardOptions::Implementation::Discrete
        };
        REQUIRE(options.getImplementation() ==
                UFT::ForwardOptions::Implementation::Discrete);
    }
    SECTION("FFT")
    {
        UFT::ForwardOptions options;
        options.setImplementation(
            UFT::ForwardOptions::Implementation::Fast);
        REQUIRE(options.getImplementation() ==
                UFT::ForwardOptions::Implementation::Fast);
        UFT::ForwardOptions copy{options};
        REQUIRE(copy.getImplementation() ==
                UFT::ForwardOptions::Implementation::Fast);
    }
}

TEMPLATE_TEST_CASE("CoreTest::Transforms::Fourier::Forward::DFT",
                  "[TypeName][template]", float, double)
{
    UFT::ForwardOptions options{UFT::ForwardOptions::Implementation::Discrete};
    UFT::Forward<TestType> dft{options};
    SECTION("Real Length 5")
    {
        // ft = fft.rfft([1, 2, 3, 4, 5])
        USignal::Vector<TestType> x(std::vector<TestType> {1, 2, 3, 4, 5});
        USignal::Vector<std::complex<TestType>> yRef{
            std::vector<std::complex<TestType>> { 15  + 0i, 
                                                 -2.5 + 3.4409548011779334i,
                                                 -2.5 + 0.8122992405822659i} };
        REQUIRE_NOTHROW(dft.setInput(x));
        REQUIRE_NOTHROW(dft.apply());
        auto result = dft.getOutput();
        REQUIRE(result.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            auto residual = std::abs(yRef.at(i) - result[i]);
            REQUIRE(residual < std::numeric_limits<TestType>::epsilon()*10);
        } 
    }
    SECTION("Real Length 6")
    {
        // ft = fft.rfft([1, 2, 3, -1, -2, -3])
        USignal::Vector<TestType> x(std::vector<TestType> {1, 2, 3, -1, -2, -3.1});
        USignal::Vector<std::complex<TestType>> yRef{
           std::vector<std::complex<TestType>> {-0.09999999999999964 + 0i,
                                                 0.9500000000000002 - 8.746856578222829i,
                                                 0.04999999999999982 - 0.0866025403784434i,
                                                 4.1 + 0i}};
        REQUIRE_NOTHROW(dft.setInput(x));
        REQUIRE_NOTHROW(dft.apply());
        auto result = dft.getOutput();
        REQUIRE(result.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            auto residual = std::abs(yRef.at(i) - result[i]);
            REQUIRE(residual < std::numeric_limits<TestType>::epsilon()*10);
        }
    }

    SECTION("Real Length 8")
    {
        USignal::Vector<TestType> x(std::vector<TestType> {1, 2, 3, -1, -2, -3.1, 0, 0});
        USignal::Vector<std::complex<TestType>> yRef{
           std::vector<std::complex<TestType>> {-0.10000000000000009 + 0i, 
                                                 7.313351365237939 - 5.8991378028648445i,
                                                -4 + 0.10000000000000009i,
                                                 -1.3133513652379394 + 0.1008621971351551i,
                                                 4.1 + 0i}};
        REQUIRE_NOTHROW(dft.setInput(x));
        REQUIRE_NOTHROW(dft.apply());
        auto result = dft.getOutput();
        REQUIRE(result.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {   
            auto residual = std::abs(yRef.at(i) - result.at(i));
            REQUIRE(residual < std::numeric_limits<TestType>::epsilon()*10);
        }
    }
}

TEMPLATE_TEST_CASE("CoreTest::Transforms::Fourier::Forward::FFT",
                  "[TypeName][template]", double)
{
    UFT::ForwardOptions options{UFT::ForwardOptions::Implementation::Fast};
    UFT::Forward<TestType> dft{options};
    SECTION("Real Length 5")
    {   
        // ft = fft.rfft([1, 2, 3, -1, -2, -3], 8)
        USignal::Vector<TestType> x(std::vector<TestType> {1, 2, 3, -1, -2, -3.1});
        USignal::Vector<std::complex<TestType>> yRef{
           std::vector<std::complex<TestType>> {-0.10000000000000009 + 0i,
                                                 7.313351365237939 - 5.8991378028648445i,
                                                -4 + 0.10000000000000009i,
                                                 -1.3133513652379394 + 0.1008621971351551i,
                                                 4.1 + 0i}};
        REQUIRE_NOTHROW(dft.setInput(x));
        REQUIRE_NOTHROW(dft.apply());
        auto result = dft.getOutput();
        REQUIRE(result.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            auto residual = std::abs(yRef.at(i) - result.at(i));
            REQUIRE(residual < std::numeric_limits<TestType>::epsilon()*10);
        }
    }
}
