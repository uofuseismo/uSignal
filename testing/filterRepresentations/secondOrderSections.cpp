#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <array>
#include "uSignal/filterRepresentations/secondOrderSections.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;

TEMPLATE_TEST_CASE("CoreTest::FilterRepresentations::SecondOrderSections",
                   "[TypeName][template]",
                   int, float, double)
{
    USignal::Vector<TestType> numeratorCoefficients({3, 6, 9, 2, 5, 8, 1, 4, 7});
    USignal::Vector<TestType> denominatorCoefficients({-3, -6, -9, -2, -5, -8, -1, -4, -7});

    SECTION("Real Vector Interface")
    {
        USignal::FilterRepresentations::SecondOrderSections<TestType> 
           sos{numeratorCoefficients, denominatorCoefficients};
        REQUIRE(sos.getNumberOfSections() == 3);
        auto bs = sos.getNumeratorFilterCoefficientsReference();
        auto as = sos.getDenominatorFilterCoefficientsReference();
        REQUIRE(bs.size() == numeratorCoefficients.size());
        REQUIRE(as.size() == denominatorCoefficients.size());
        for (int i = 0; i < static_cast<int> (bs.size()); ++i)
        {
            CHECK(bs[i] == Catch::Approx(numeratorCoefficients[i]));
        }
        for (int i = 0; i < static_cast<int> (as.size()); ++i)
        {
            CHECK(as[i] == Catch::Approx(denominatorCoefficients[i]));
        }
    }
    SECTION("Array Interface")
    {
        std::vector<std::array<TestType, 3>> 
           numeratorCoefficientsArray( {{3, 6, 9},
                                        {2, 5, 8},
                                        {1, 4, 7}} );
        std::vector<std::array<TestType, 3>>
           denominatorCoefficientsArray( {{-3, -6, -9},
                                          {-2, -5, -8},
                                          {-1, -4, -7}} );
        USignal::FilterRepresentations::SecondOrderSections<TestType> 
           sos{numeratorCoefficientsArray, denominatorCoefficientsArray};
 
        auto bs = sos.getNumeratorFilterCoefficients();
        auto as = sos.getDenominatorFilterCoefficients();
        REQUIRE(bs.size() == numeratorCoefficients.size());
        REQUIRE(as.size() == denominatorCoefficients.size());
        for (int i = 0; i < static_cast<int> (bs.size()); ++i)
        {
            CHECK(bs[i] == Catch::Approx(numeratorCoefficients[i]));
        }
        for (int i = 0; i < static_cast<int> (as.size()); ++i)
        {
            CHECK(as[i] == Catch::Approx(denominatorCoefficients[i]));
        }
    }
}
