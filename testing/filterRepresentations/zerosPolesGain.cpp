#include <vector>
#include <complex>
#include <cmath>
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEMPLATE_TEST_CASE("CoreTest::FilterRepresentations::ZerosPolesGain",
                   "[TypeName][template]",
                   double, float)
{
    USignal::Vector<std::complex<TestType>> zerosRef(
        std::vector<std::complex<TestType>> { 1 + 0i, 
                                              1 + 0i, 
                                              1 + 0i, 
                                              1 + 0i, 
                                             -1 + 0i, 
                                             -1 + 0i, 
                                             -1 + 0i, 
                                             -1 + 0i} );
    USignal::Vector<std::complex<TestType>> polesRef(
        std::vector<std::complex<TestType>> {
            0.861419077078615-0.40475046563708234i,
            0.7609277750082251-0.4142205574875641i,
            0.7609277750082251+0.4142205574875641i,
            0.861419077078615+0.40475046563708234i,
            0.6708198383910375+0.634395173212654i,
            0.6746102884496724+0.4961465975166494i,
            0.6746102884496724-0.4961465975166494i,
            0.6708198383910375-0.634395173212654i} );
    constexpr TestType gainRef{0.0005705645409457355};
    USignal::FilterRepresentations::ZerosPolesGain
        zpk{zerosRef, polesRef, gainRef};

    SECTION("Preliminary")
    {
        auto zeros = zpk.getZeros();
        auto poles = zpk.getPoles();
        auto gain = zpk.getGain();
        REQUIRE(zerosRef.size() == zeros.size());
        for (int i = 0; i < static_cast<int> (zeros.size()); ++i)
        {
            auto residual = std::abs(zerosRef.at(i) - zeros.at(i));
            CHECK(Catch::Approx(0).margin(
                  std::numeric_limits<TestType>::epsilon()*100) == residual);
        }
        REQUIRE(polesRef.size() == poles.size());
        for (int i = 0; i < static_cast<int> (poles.size()); ++i)
        {
            auto residual = std::abs(polesRef.at(i) - poles.at(i));
            CHECK(Catch::Approx(0).margin(
                  std::numeric_limits<TestType>::epsilon()*100) == residual);
        }   
        CHECK(Catch::Approx(gainRef).margin(
             std::numeric_limits<TestType>::epsilon()*100) == gain);
     }

    SECTION("Copy")
    {   
        USignal::FilterRepresentations::ZerosPolesGain copy{zpk};
        const auto zeros = copy.getZerosReference();
        const auto poles = copy.getPolesReference();
        auto gain = zpk.getGain();
        REQUIRE(zerosRef.size() == zeros.size());
        for (int i = 0; i < static_cast<int> (zeros.size()); ++i)
        {
            auto residual = std::abs(zerosRef.at(i) - zeros.at(i));
            CHECK(Catch::Approx(0).margin(
                  std::numeric_limits<TestType>::epsilon()*100) == residual);
        }
        REQUIRE(polesRef.size() == poles.size());
        for (int i = 0; i < static_cast<int> (poles.size()); ++i)
        {
            auto residual = std::abs(polesRef.at(i) - poles.at(i));
            CHECK(Catch::Approx(0).margin(
                  std::numeric_limits<TestType>::epsilon()*100) == residual);
        }
        CHECK(Catch::Approx(gainRef).margin(
             std::numeric_limits<TestType>::epsilon()*100) == gain);
     }   
}
