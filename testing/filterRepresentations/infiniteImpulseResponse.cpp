#include <iostream>
#include <complex>
#include <cmath>
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;

TEMPLATE_TEST_CASE("CoreTest::FilterRepresentations::InfiniteImpulseResponse",
                   "[TypeName][template]",
                   int, float, double)
{
    int nNumeratorCoefficients{5}; 
    USignal::Vector<TestType> numeratorCoefficients;
    for (int ic = 0; ic < nNumeratorCoefficients; ++ic)
    {
        numeratorCoefficients.push_back(ic);
    }
    int nDenominatorCoefficients{6};
    USignal::Vector<TestType> denominatorCoefficients;
    for (int ic = 0; ic < nDenominatorCoefficients; ++ic)
    {
        denominatorCoefficients.push_back(-ic - 1);
    }
    int order = std::max(nNumeratorCoefficients, nDenominatorCoefficients) - 1;


    USignal::FilterRepresentations::InfiniteImpulseResponse<TestType>
         iir{numeratorCoefficients, denominatorCoefficients};
    auto b = iir.getNumeratorFilterCoefficients(); 
    auto a = iir.getDenominatorFilterCoefficients();
    REQUIRE(order == iir.getOrder());
    REQUIRE(b.size() == numeratorCoefficients.size());
    REQUIRE(a.size() == denominatorCoefficients.size());
    for (int i = 0; i < static_cast<int> (b.size()); ++i)
    {
        CHECK(b[i] == Catch::Approx(numeratorCoefficients[i]));
    }
    for (int i = 0; i < static_cast<int> (a.size()); ++i)
    {
        CHECK(a[i] == Catch::Approx(denominatorCoefficients[i]));
    }

    SECTION("Copy")
    {
        auto copy = iir;
        auto bCopy = copy.getNumeratorFilterCoefficientsReference();
        auto aCopy = copy.getDenominatorFilterCoefficientsReference();
        REQUIRE(order == copy.getOrder());
        REQUIRE(bCopy.size() == numeratorCoefficients.size());
        REQUIRE(aCopy.size() == denominatorCoefficients.size());
        for (int i = 0; i < static_cast<int> (b.size()); ++i)
        {
            CHECK(b[i] == Catch::Approx(numeratorCoefficients[i]));
        }
        for (int i = 0; i < static_cast<int> (a.size()); ++i)
        {
            CHECK(a[i] == Catch::Approx(denominatorCoefficients[i]));
        }
    }
}

