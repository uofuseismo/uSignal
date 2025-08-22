#include <iostream>
#include <filesystem>
#include <fstream>
#include <complex>
#include <vector>
#include <chrono>
#include <cmath>
#include "uSignal/filterImplementations/transposeDirectForm2.hpp"
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
std::filesystem::path highOrderIIRFileName{dataDirectory/"transposeDF2.gse2.txt"};
//std::filesystem::path lowOrderIIRFileName{dataDirectory/"iirReference2.txt"};

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
        double yi;
        std::sscanf(line.c_str(), "%lf\n", &yi);
        signal.push_back(static_cast<T> (yi));
    }
    signalFile.close();
    return signal;
}
}

TEMPLATE_TEST_CASE("CoreTest::FilterImplementations::TransposeDirectForm2",
                   "[TypeName][template]",
                   double)
{
    auto inputSignal = ::loadSignal<TestType> (inputSignalFileName);

    SECTION("High Order")
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

        auto yRef = ::loadSignal<TestType> (highOrderIIRFileName);

        auto startTime = std::chrono::high_resolution_clock::now();
        USignal::FilterImplementations::TransposeDirectForm2<TestType>
            iir{iirFilterCoefficients};
        REQUIRE_NOTHROW(iir.setInput(inputSignal));
        //auto startTime = std::chrono::high_resolution_clock::now();
        REQUIRE_NOTHROW(iir.apply());
        //auto endTime = std::chrono::high_resolution_clock::now();
        auto outputSignal = iir.getOutput();

        auto endTime = std::chrono::high_resolution_clock::now(); 
        auto elapsedTime
            = std::chrono::duration_cast<std::chrono::microseconds>
              (endTime - startTime).count()*1.e-6;
        std::cout << "High order processing time: " << elapsedTime << std::endl;

        REQUIRE(outputSignal.size() == yRef.size());
        TestType l1Error{0};
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            auto residual = std::abs(yRef[i] - outputSignal[i]);
            l1Error = std::max(l1Error, residual);
        }
        REQUIRE(l1Error < std::numeric_limits<TestType>::epsilon()*100);
    }

    SECTION("Order 0 Test")
    {
        Vector<TestType> b; b.push_back(10);
        Vector<TestType> a; a.push_back(5);
        constexpr TestType ratio{2}; // 10/5
        //a[0] y[i] = b[0] x[i] -> y[i] = b[0]/a[0] x[i]  
        USignal::FilterRepresentations::InfiniteImpulseResponse<TestType>
            iirFilterCoefficients{b, a};
        auto startTime = std::chrono::high_resolution_clock::now();
        USignal::FilterImplementations::TransposeDirectForm2<TestType>
            iir{iirFilterCoefficients};
        REQUIRE_NOTHROW(iir.setInput(inputSignal));
        //auto startTime = std::chrono::high_resolution_clock::now();
        REQUIRE_NOTHROW(iir.apply());
        //auto endTime = std::chrono::high_resolution_clock::now();
        auto outputSignal = iir.getOutput();
        auto endTime = std::chrono::high_resolution_clock::now();
        auto elapsedTime
            = std::chrono::duration_cast<std::chrono::microseconds>
              (endTime - startTime).count()*1.e-6;
        std::cout << "High order processing time: " << elapsedTime << std::endl;

        REQUIRE(outputSignal.size() == inputSignal.size());

        TestType l1Error{0};
        for (int i = 0; i < static_cast<int> (inputSignal.size()); ++i)
        {
            auto residual = std::abs(outputSignal[i] - ratio*inputSignal[i]);
            l1Error = std::max(l1Error, residual);
        }
        REQUIRE(l1Error < std::numeric_limits<TestType>::epsilon()*100);
    }
}
