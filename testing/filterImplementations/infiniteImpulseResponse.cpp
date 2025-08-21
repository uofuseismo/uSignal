#include <iostream>
#include <filesystem>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include "uSignal/filterImplementations/infiniteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;

namespace
{

std::filesystem::path dataDirectory{"data"};
std::filesystem::path inputSignalFileName{dataDirectory/"gse2.txt"};
std::filesystem::path highOrderIIRFileName{dataDirectory/"iirReference1.txt"};
std::filesystem::path lowOrderIIRFileName{dataDirectory/"iirReference2.txt"};

template<typename T>
Vector<T> loadSignal(const std::filesystem::path &fileName)
{
    if (!std::filesystem::exists(fileName))
    {
        throw std::invalid_argument(std::string {fileName} + " does not exist");
    }
    std::ifstream signalFile;
    signalFile.open(fileName);
    Vector<T> signal;
    signal.reserve(12000);
    std::string line;
    while (std::getline(signalFile, line))
    {
        T yi;
        std::sscanf(line.c_str(), "%lf\n", &yi);
        signal.push_back(yi);
    }
    signalFile.close();
    return signal;
}
}

TEMPLATE_TEST_CASE("CoreTest::FilterImplementations::InfiniteImpulseResponse",
                   "[TypeName][template]",
                   double)
{
    Vector<TestType> b({0.000401587491686, 
                        0.0,
                       -0.001606349966746,
                        0.0,
                        0.002409524950119,
                        0.0,
                       -0.001606349966746,
                        0.0,
                        0.000401587491686});
    Vector<TestType> a({ 1.000000000000000,
                        -7.185226122700763,
                        22.615376628798678,
                       -40.733465892344896,
                        45.926605646620146,
                       -33.196326377161412,
                        15.023103545324197,
                        -3.891997997268024,
                         0.441930568732716});
    USignal::FilterRepresentations::InfiniteImpulseResponse<TestType>
         iirFilterCoefficients{b, a};
    USignal::FilterImplementations::InfiniteImpulseResponse<TestType>
         iir{iirFilterCoefficients};
    auto inputSignal = ::loadSignal<TestType> (inputSignalFileName);
    auto yRef = ::loadSignal<TestType> (highOrderIIRFileName);
    REQUIRE_NOTHROW(iir.setInput(inputSignal));
    REQUIRE_NOTHROW(iir.apply());
    auto outputSignal = iir.getOutput();
/*
    REQUIRE(outputSignal.size() == yRef.size());
    double l2Error{0};
    for (int i = 0; i < yRef.size(); ++i)
    {
        auto residual = std::abs(yRef[i] - outputSignal[i]);
        l2Error = l2Error + residual*residual;
    }
    l2Error = std::sqrt(l2Error);
    auto rmse = l2Error/yRef.size();
    REQUIRE(rmse < 0.05);
//std::cout << "rmse: " << std::setprecision(12) << l2Error/yRef.size() << std::endl;
*/
 for (int i = yRef.size() - 15; i < yRef.size(); ++i)
// for (int i = 0; i < 25; ++i)
 {
 //    std::cout << outputSignal[i] << " " << yRef[i] << std::endl;
 }
}
