#include <iostream>
#include <complex>
#include <cmath>
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;

TEMPLATE_TEST_CASE("CoreTest::FilterRepresentations::FiniteImpulseResponse",
                   "[Real][TypeName][template]",
                   int, float, double)
{
    int nCoefficients{5}; 
    USignal::Vector<TestType> filterCoefficients;
    for (int ic = 0; ic < nCoefficients; ++ic)
    {
        filterCoefficients.push_back(ic);
    }

    USignal::FilterRepresentations::FiniteImpulseResponse<TestType>
         fir{filterCoefficients};
    REQUIRE(fir.getOrder() == nCoefficients - 1);
 
    auto copy = fir.getFilterCoefficientsReference();
    REQUIRE(copy.size() == filterCoefficients.size());
    for (int i = 0; i < static_cast<int> (copy.size()); ++i)
    {   
        CHECK(copy[i] == Catch::Approx(filterCoefficients[i]));
    }   

/*
    // Want to make sure these extra coefficients get removed
    SECTION("Zero Pad")
    {
        USignal::Vector<TestType> filterCoefficientsPadded;
        filterCoefficientsPadded.push_back(0);
        filterCoefficientsPadded.push_back(0);
        for (int ic = 0; ic < nCoefficients; ++ic)
        {
            filterCoefficients.push_back(ic + 1);
        }
        filterCoefficientsPadded.push_back(0);
        filterCoefficientsPadded.push_back(0);
        filterCoefficientsPadded.push_back(0);
        USignal::FilterRepresentations::FiniteImpulseResponse<TestType>
            firPadded{filterCoefficientsPadded};
        REQUIRE(firPadded.getOrder() == nCoefficients - 1); 
        auto copy2 = firPadded.getFilterCoefficientsReference();
        REQUIRE(copy2.size() == filterCoefficients.size());
        for (int i = 0; i < static_cast<int> (copy2.size()); ++i)
        {
            CHECK(copy2[i] == Catch::Approx(filterCoefficients[i]));
        }   
    }
*/
}

TEMPLATE_TEST_CASE("CoreTest::FilterRepresentations::FiniteImpulseResponse",
                   "[Complex][TypeName][template]",
                   std::complex<float>,
                   std::complex<double>)
{
    int nCoefficients{5}; 
    USignal::Vector<TestType> filterCoefficients;
    for (int ic = 0; ic < nCoefficients - 1; ++ic)
    {
        TestType b = ic + 1i;
        filterCoefficients.push_back( std::exp(b) );
    }
    filterCoefficients.push_back(0 + 1i); // Make sure 0-equality check works

/*
    USignal::FilterRepresentations::FiniteImpulseResponse<TestType>
         fir{filterCoefficients};
    REQUIRE(fir.getOrder() == nCoefficients - 1); 
    auto copy = fir.getFilterCoefficients();
    REQUIRE(copy.size() == filterCoefficients.size());
    for (int i = 0; i < static_cast<int> (copy.size()); ++i)
    {
        auto residual 
            = static_cast<double> (std::abs(copy[i] - filterCoefficients[i]));
        CHECK(residual < std::numeric_limits<double>::epsilon());
    }


    SECTION("Zero Pad")
    {
        auto filterCoefficientsPadded = filterCoefficients;
        filterCoefficientsPadded.push_back(0 + 0i);
        filterCoefficientsPadded.push_back(0 + 0i);
        USignal::FilterRepresentations::FiniteImpulseResponse<TestType>
            firPadded{filterCoefficientsPadded};
        REQUIRE(firPadded.getOrder() == nCoefficients - 1); 
        auto copy2 = firPadded.getFilterCoefficientsReference();
        REQUIRE(copy2.size() == filterCoefficients.size());
        for (int i = 0; i < static_cast<int> (copy.size()); ++i)
        {
            auto residual 
                = static_cast<double>
                  (std::abs(copy2[i] - filterCoefficients[i]));
            CHECK(residual < std::numeric_limits<double>::epsilon());
        }
    }   
*/
}

