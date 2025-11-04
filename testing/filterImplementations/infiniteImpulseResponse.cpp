#include <iostream>
#include <filesystem>
#include <fstream>
#include <complex>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include "uSignal/filterImplementations/transposeDirectForm2.hpp"
#include "uSignal/filterImplementations/secondOrderSections.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/filterRepresentations/secondOrderSections.hpp"
#include "uSignal/vector.hpp"
#include "loadSignal.hpp"
#include "infinityNorm.hpp"
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
std::filesystem::path sosFileName{dataDirectory/"sosReference.gse2.txt"};
//std::filesystem::path lowOrderIIRFileName{dataDirectory/"iirReference2.txt"};

/*
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

template<typename T>
T computeInfinityNorm(const USignal::Vector<T> &yTrue,
                      const USignal::Vector<T> &yEstimate)
{
    if (yTrue.size() != yEstimate.size())
    {
        throw std::invalid_argument("Inconsistent input sizes");
    }
    T error{0};
    for (int i = 0; i < static_cast<int> (yTrue.size()); ++i)
    {
        error = std::max(error, std::abs(yTrue[i] - yEstimate[i]));
    }
    return error;
}
*/

}

TEMPLATE_TEST_CASE("CoreTest::FilterImplementations::SecondOrderSections",
                   "[TypeName][template]",
                   double)
{
    SECTION("Impulse Response PostProcessing")
    {
        Vector<TestType> bs(
            std::vector<TestType> {6.37835424e-05,  6.37835424e-05,  0.00000000e+00,
                                   1.00000000e+00, -1.78848938e+00,  1.00000000e+00,
                                   1.00000000e+00, -1.93118487e+00,  1.00000000e+00,
                                   1.00000000e+00, -1.95799864e+00,  1.00000000e+00,
                                   1.00000000e+00, -1.96671846e+00,  1.00000000e+00,
                                   1.00000000e+00, -1.97011885e+00,  1.00000000e+00,
                                   1.00000000e+00, -1.97135784e+00,  1.00000000e+00} );
        Vector<TestType> as(
             std::vector<TestType> {1.00000000e+00, -9.27054679e-01,  0.00000000e+00,
                                    1.00000000e+00, -1.87008942e+00,  8.78235919e-01,
                                    1.00000000e+00, -1.90342568e+00,  9.17455718e-01,
                                    1.00000000e+00, -1.93318668e+00,  9.52433552e-01,
                                    1.00000000e+00, -1.95271141e+00,  9.75295685e-01,
                                    1.00000000e+00, -1.96423610e+00,  9.88608056e-01,
                                    1.00000000e+00, -1.97157693e+00,  9.96727086e-01} );
        Vector<TestType> yRef(
             std::vector<TestType> {6.37835424e-05,  1.23511272e-04,  1.34263690e-04,
                                    1.78634911e-04,  2.50312740e-04,  3.46332848e-04,
                                    4.66239952e-04,  6.11416691e-04,  7.84553129e-04,
                                    9.89232232e-04,  1.22960924e-03,  1.51016546e-03,
                                    1.83551947e-03,  2.21028135e-03,  2.63893773e-03,
                                    3.12575784e-03,  3.67471270e-03,  4.28940130e-03,
                                    4.97297977e-03,  5.72809028e-03,  6.55678845e-03,
                                    7.46046851e-03,  8.43978671e-03,  9.49458408e-03,
                                    1.06238101e-02,  1.18254496e-02,  1.30964547e-02,
                                    1.44326848e-02,  1.58288573e-02,  1.72785101e-02,
                                    1.87739799e-02,  2.03063976e-02,  2.18657022e-02,
                                    2.34406756e-02,  2.50189979e-02,  2.65873261e-02,
                                    2.81313940e-02,  2.96361349e-02,  3.10858256e-02,
                                    3.24642512e-02} );
        USignal::FilterRepresentations::SecondOrderSections
            filterCoefficients{ bs, as };
        USignal::FilterImplementations::SecondOrderSections
            sos{filterCoefficients};
        // Create an impulse
        Vector<TestType> impulse(40, 0); 
        impulse.at(0) = 1;  
        REQUIRE(sos.isInitialized());
        REQUIRE_NOTHROW(sos.setInput(impulse));
        REQUIRE_NOTHROW(sos.apply());
        auto y = sos.getOutput();
        REQUIRE(y.size() == yRef.size());
        auto tolerance = std::max(std::numeric_limits<TestType>::epsilon()*10,
                                  1.e-8);
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        { 
            CHECK(std::abs(y[i] - yRef[i]) < tolerance);
        }
        REQUIRE_NOTHROW(sos.apply()); 
        y = sos.getOutput();
        REQUIRE(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(std::abs(y[i] - yRef[i]) < tolerance);
        }
    }

    auto inputSignal = ::loadSignal<TestType> (inputSignalFileName);
    SECTION("Post Processing")
    {
        Vector<TestType> bs(
            std::vector<TestType> {0.000401587491686,  0.000803175141692,  0.000401587491549,
                                   1.000000000000000, -2.000000394412897,  0.999999999730209,
                                   1.000000000000000,  1.999999605765104,  1.000000000341065,
                                   1.000000000000000, -1.999999605588274,  1.000000000269794} );
        Vector<TestType> as(
            std::vector<TestType> {1.000000000000000, -1.488513049541281,  0.562472929601870,
                                   1.000000000000000, -1.704970593447777,  0.792206889942566,
                                   1.000000000000000, -1.994269533089365,  0.994278822534674,
                                   1.000000000000000, -1.997472946622339,  0.997483252685326} );
        USignal::FilterRepresentations::SecondOrderSections
            filterCoefficients{ bs, as };

        auto yRef = ::loadSignal<TestType> (sosFileName);

        auto startTime = std::chrono::high_resolution_clock::now();
        USignal::FilterImplementations::SecondOrderSections
            sos{filterCoefficients};
        REQUIRE_NOTHROW(sos.setInput(inputSignal));
        REQUIRE_NOTHROW(sos.apply());
        const auto outputSignal = sos.getOutputReference();

        auto endTime = std::chrono::high_resolution_clock::now();
        auto elapsedTime
            = std::chrono::duration_cast<std::chrono::microseconds>
              (endTime - startTime).count()*1.e-6;
        std::cout << "SOS processing time: " << elapsedTime << std::endl;
        auto l8Norm = ::computeInfinityNorm<TestType> (yRef, outputSignal);
        REQUIRE(l8Norm < 1.e-9);
    }

    SECTION("Real Time")
    {
        Vector<TestType> bs( 
            std::vector<TestType> {0.000401587491686,  0.000803175141692,  0.000401587491549,
                                   1.000000000000000, -2.000000394412897,  0.999999999730209,
                                   1.000000000000000,  1.999999605765104,  1.000000000341065,
                                   1.000000000000000, -1.999999605588274,  1.000000000269794} );
        Vector<TestType> as( 
            std::vector<TestType> {1.000000000000000, -1.488513049541281,  0.562472929601870,
                                   1.000000000000000, -1.704970593447777,  0.792206889942566,
                                   1.000000000000000, -1.994269533089365,  0.994278822534674,
                                   1.000000000000000, -1.997472946622339,  0.997483252685326} );
        USignal::FilterRepresentations::SecondOrderSections
            filterCoefficients{ bs, as };

        auto yRef = ::loadSignal<TestType> (sosFileName);

        USignal::FilterImplementations::SecondOrderSections
            sos{filterCoefficients, true};

        std::mt19937 generator(83823);
        std::uniform_int_distribution<int> distribution{1, 100};
        USignal::Vector<TestType> outputSignal(yRef.size());
        int kStart{0};
        auto nSamples = static_cast<int> (inputSignal.size());
        for (int k = 0; k < nSamples; ++k)
        {   
            int nSamplesInPacket = distribution(generator);
            auto kEnd = std::min(kStart + nSamplesInPacket, nSamples);
            nSamplesInPacket = kEnd - kStart; // Fix packet length 
            if (nSamplesInPacket <= 0){break;}
            USignal::Vector<TestType> packet(nSamplesInPacket);
            std::copy(inputSignal.data() + kStart, inputSignal.data() + kEnd,
                      packet.data());
            REQUIRE_NOTHROW(sos.setInput(packet));
            REQUIRE_NOTHROW(sos.apply());
            auto outputPacket = sos.getOutput();
            std::copy(outputPacket.data(),
                      outputPacket.data() + nSamplesInPacket,
                      outputSignal.data() + kStart);

            kStart = kStart + nSamplesInPacket;
            if (kStart >= nSamples){break;}
        }   
        auto l8Norm = ::computeInfinityNorm<TestType> (yRef, outputSignal);
        REQUIRE(l8Norm < 1.e-9);


    }
}

///--------------------------------------------------------------------------///

TEMPLATE_TEST_CASE("CoreTest::FilterImplementations::TransposeDirectForm2",
                   "[TypeName][template]",
                   double)
{
    auto inputSignal = ::loadSignal<TestType> (inputSignalFileName);

    SECTION("High Order Post Processing")
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
        auto outputSignal = iir.getOutputReference();

        auto endTime = std::chrono::high_resolution_clock::now(); 
        auto elapsedTime
            = std::chrono::duration_cast<std::chrono::microseconds>
              (endTime - startTime).count()*1.e-6;
        std::cout << "High order processing time: " << elapsedTime << std::endl;

        REQUIRE(outputSignal.size() == yRef.size());
        auto l8Error = ::computeInfinityNorm<TestType> (yRef, outputSignal);
        REQUIRE(l8Error < std::numeric_limits<TestType>::epsilon()*100);
    }

    SECTION("Order 0 Post Processing")
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
        const auto outputSignal = iir.getOutput();
        auto endTime = std::chrono::high_resolution_clock::now();
        auto elapsedTime
            = std::chrono::duration_cast<std::chrono::microseconds>
              (endTime - startTime).count()*1.e-6;
        std::cout << "Order 0 processing time: " << elapsedTime << std::endl;

        REQUIRE(outputSignal.size() == inputSignal.size());

        TestType l8Error{0};
        for (int i = 0; i < static_cast<int> (inputSignal.size()); ++i)
        {
            auto residual = std::abs(outputSignal[i] - ratio*inputSignal[i]);
            l8Error = std::max(l8Error, residual);
        }
        REQUIRE(l8Error < std::numeric_limits<TestType>::epsilon()*100);
    }
}
