#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <array>
#include "uSignal/filterRepresentations/secondOrderSections.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;

TEMPLATE_TEST_CASE("CoreTest::FilterRepresentations::SecondOrderSections",
                   "[TypeName][template]",
                   double, float)
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
    SECTION("From Transfer Function")
    {
        USignal::Vector<TestType> bs
        {
            std::vector<TestType> {0.9877613892768257,
                                  -3.9510455571073027,
                                   5.926568335660954,
                                  -3.9510455571073027,
                                   0.9877613892768257}
        };
        USignal::Vector<TestType> as
        {
            std::vector<TestType> {1.0,
                                  -3.97537191256092,
                                   5.926418555965423,
                                  -3.926719197756784,
                                   0.975672562146085}
        };
        USignal::Vector<TestType> bsRef
        {
            std::vector<TestType> {
                                   0.9877613892768257, 
                                  -1.9755228135715934,
                                   0.987761379675909,
                                   1.,
                                  -1.9999999645481794,
                                   1.000000009719877
                                  }
        };
        USignal::Vector<TestType> asRef
        {
            std::vector<TestType> {1.0, -1.982647799569223, 0.9827358586551028,
                                   1.0, -1.9927241129916946,0.9928126195387987}
        };
        USignal::FilterRepresentations::InfiniteImpulseResponse<TestType>
            ba{bs, as};
        const auto pairingStrategy
        {
            USignal::FilterRepresentations::
              SecondOrderSections<TestType>::PairingStrategy::None //earest
        };
        USignal::FilterRepresentations::SecondOrderSections<TestType>
            sos(ba, pairingStrategy); 
        auto bsos = sos.getNumeratorFilterCoefficients();
        auto asos = sos.getDenominatorFilterCoefficients();
        REQUIRE(bsos.size() == bsRef.size());
        for (int i = 0; i < static_cast<int> (bsRef.size()); ++i)
        {
            REQUIRE(std::abs(bsos[i] - bsRef[i]) < 1.e-7);
        }
        REQUIRE(asos.size() == asRef.size());
        const auto tol = std::same_as<TestType, double> ? 1.e-8 : 1.e-1;
        for (int i = 0; i < static_cast<int> (asRef.size()); ++i)
        {
            REQUIRE(std::abs(asos[i] - asRef[i]) < tol);
            //std::cout << "compare " << std::setprecision(12) << std::abs(asos[i] - asRef[i]) << std::endl; //" " << asRef[i] << std::endl; //std::abs(asos[i] - asRef[i]) << std::endl;
        }
    }
}
