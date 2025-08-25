#include <iostream>
#include <filesystem>
#include <fstream>
#include <complex>
#include <vector>
#include <chrono>
#include <cmath>
#include "uSignal/filterImplementations/finiteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/finiteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include "loadSignal.hpp"
#include "infinityNorm.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;

TEMPLATE_TEST_CASE("CoreTest::FilterImplementations::FiniteImpulseResponse",
                   "[TypeName][template]",
                   double)
{
    std::filesystem::path dataDirectory{"data"};
    std::filesystem::path inputSignalFileName{dataDirectory/"gse2.txt"};
    std::filesystem::path firFileName{dataDirectory/"firReference.gse2.txt"};

    auto inputSignal = ::loadSignal<TestType> (inputSignalFileName);
    SECTION("Post Processing")
    {
        Vector<TestType> bs( 
            std::vector<TestType> {-0.000000000000000, -0.001056235801065,
                                   -0.000769341020163,  0.000956323223723,
                                    0.001976082742122, -0.000000000000000,
                                   -0.003265384800345, -0.002568519852901,
                                    0.003234633130890,  0.006519908075213,
                                   -0.000000000000000, -0.009825739114867,
                                   -0.007365685405410,  0.008881348924986,
                                    0.017256056989442, -0.000000000000000,
                                   -0.024784271698734, -0.018417666768131,
                                    0.022299534288278,  0.044222443880910,
                                   -0.000000000000000, -0.071469809226860,
                                   -0.060430328816090,  0.092317626953209,
                                    0.302027315266443,  0.400523418058701,
                                    0.302027315266443,  0.092317626953209,
                                   -0.060430328816090, -0.071469809226860,
                                   -0.000000000000000,  0.044222443880910,
                                    0.022299534288278, -0.018417666768131,
                                   -0.024784271698734, -0.000000000000000,
                                    0.017256056989442,  0.008881348924986,
                                   -0.007365685405410, -0.009825739114867,
                                   -0.000000000000000,  0.006519908075213,
                                    0.003234633130890, -0.002568519852901,
                                   -0.003265384800345, -0.000000000000000,
                                    0.001976082742122,  0.000956323223723,
                                   -0.000769341020163, -0.001056235801065,
                                   -0.000000000000000});
        FilterRepresentations::FiniteImpulseResponse<TestType>
            filterCoefficients{bs};

        auto yRef = ::loadSignal<TestType> (firFileName);

        auto startTime = std::chrono::high_resolution_clock::now();
        FilterImplementations::FiniteImpulseResponse<TestType>
            fir{filterCoefficients,
                FilterImplementations::FiniteImpulseResponse<TestType>::Implementation::Direct};
        REQUIRE_NOTHROW(fir.setInput(inputSignal));
        REQUIRE_NOTHROW(fir.apply());
        auto outputSignal = fir.getOutputReference();

        auto endTime = std::chrono::high_resolution_clock::now();
        auto elapsedTime
            = std::chrono::duration_cast<std::chrono::microseconds>
              (endTime - startTime).count()*1.e-6;
        std::cout << "FIR processing time: " << elapsedTime << std::endl;

        REQUIRE(outputSignal.size() == yRef.size());
        auto l8Norm = ::computeInfinityNorm(yRef, outputSignal);
        REQUIRE(l8Norm < 1.e-12);
    }
}
