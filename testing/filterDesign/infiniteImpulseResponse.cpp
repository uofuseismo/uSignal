#include <iostream>
#include <limits>
#include <complex>
#include <cmath>
#include <vector>
#include "uSignal/filterDesign/infiniteImpulseResponse/bilinearTransform.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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
