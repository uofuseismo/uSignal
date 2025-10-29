#include <cmath>
#include <complex>
#include <random>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include "uSignal/vector.hpp"
#include "uSignal/transforms/hilbert/finiteImpulseResponse.hpp"
#include "uSignal/filterDesign/finiteImpulseResponse/hilbertTransformer.hpp"
#include "uSignal/filterImplementations/finiteImpulseResponse.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;
namespace UHT = USignal::Transforms::Hilbert;

namespace
{

template<typename T>
USignal::Vector<std::complex<T>>
    createReference(const USignal::Vector<T> &x,
                    const int order, const double beta,
                    const bool isRealTime = false)
{
    
    // Now build the filter the slow way
    namespace UFD = USignal::FilterDesign::FiniteImpulseResponse;
    auto firFilter
        = UFD::hilbertTransformer<T> (order, beta).getFilterCoefficients();
    USignal::Vector<T> realTaps(firFilter.size());
    USignal::Vector<T> imagTaps(firFilter.size());
    for (int i = 0; i < static_cast<int> (firFilter.size()); ++i)
    {
        realTaps[i] = std::real(firFilter[i]);
        imagTaps[i] = std::imag(firFilter[i]);
    }
    namespace UFR = USignal::FilterRepresentations;
    UFR::FiniteImpulseResponse realFIR{realTaps};
    UFR::FiniteImpulseResponse imagFIR{imagTaps};
    namespace UFI = USignal::FilterImplementations;
    constexpr auto implementation
    {
        UFI::FiniteImpulseResponse<T>::Implementation::Direct
    };
    UFI::FiniteImpulseResponse<T>
        realFIRFilter(realFIR, implementation, isRealTime);
    realFIRFilter.setInput(x);
    realFIRFilter.apply();
    const auto yReal = realFIRFilter.getOutput();

    UFI::FiniteImpulseResponse<T>
        imagFIRFilter(imagFIR, implementation, isRealTime);
    imagFIRFilter.setInput(x);
    imagFIRFilter.apply();
    const auto yImag = imagFIRFilter.getOutput();
 
    USignal::Vector<std::complex<T>> y(yReal.size());
    for (int i = 0; i < static_cast<int> (yReal.size()); ++i)
    {
        y[i] = std::complex<T> ( yReal[i], yImag[i] );
    }
    return y;
}

}

TEST_CASE("CoreTest::Transforms::Hilbert::FiniteImpulseResponse",
          "[options]")
{
    SECTION("Defaults")
    {
        UHT::FiniteImpulseResponseOptions options;
        REQUIRE(options.getOrder() == 300);
        REQUIRE(std::abs(options.getBeta() - 8) < 1.e-14);
    }

    SECTION("Options")
    {
        const int order{59};
        constexpr double beta{6};
        UHT::FiniteImpulseResponseOptions options;
        options.setOrder(order);
        options.setBeta(beta);
        REQUIRE(options.getOrder() == order);
        REQUIRE(std::abs(options.getBeta() - beta) < 1.e-14);
    }
}

TEMPLATE_TEST_CASE("CoreTest::Transforms::Hilbert::FiniteImpulseResponse",
                  "[TypeName][template]", double)
{
    // Create a random signal
    const int nSamples{15000};
    // Apply the signal over the filter length
    std::mt19937 generator(26342);
    std::uniform_real_distribution<TestType> realDistribution{-10, 10};
    USignal::Vector<TestType> x(nSamples);
    for (int i = 0; i < nSamples; ++i)
    {
        x[i] = realDistribution(generator);
    }

    UHT::FiniteImpulseResponseOptions options;
    SECTION("Type III")
    {
        const bool isRealTime{false};
        const int order{100};
        const double beta{8};

        // Create reference
        auto yRef = createReference(x, order, beta, false);

        options.setOrder(order);
        options.setBeta(beta);

        UHT::FiniteImpulseResponse<double> hilbert{options, isRealTime};
        REQUIRE(hilbert.isInitialized());
        hilbert.setInput(x);
        hilbert.apply();
        auto yOutput = hilbert.getOutput();
        REQUIRE(yRef.size() == yOutput.size()); 
        for (int i = 0; i < yRef.size(); ++i)
        {
            auto residual = std::abs(yRef[i] - yOutput[i]);
            CHECK(residual < std::numeric_limits<TestType>::epsilon()*100);
        }
    }
 
    SECTION("Type IV")
    {
        const bool isRealTime{false};
        const int order{101};
        const double beta{6};

        // Create reference
        auto yRef = createReference(x, order, beta, false);

        options.setOrder(order);
        options.setBeta(beta);

        UHT::FiniteImpulseResponse<double> hilbert{options, isRealTime};
        REQUIRE(hilbert.isInitialized());
        hilbert.setInput(x);
        hilbert.apply();
        auto yOutput = hilbert.getOutput();
        REQUIRE(yRef.size() == yOutput.size()); 
        for (int i = 0; i < yRef.size(); ++i)
        {
            auto residual = std::abs(yRef[i] - yOutput[i]);
            CHECK(residual < std::numeric_limits<TestType>::epsilon()*100);
        }
    }

    SECTION("Type III Real Time")
    {   
        const bool isRealTime{true};
        const int order{250};
        const double beta{8};

        // Create reference
        auto yRef = createReference(x, order, beta, false);

        options.setOrder(order);
        options.setBeta(beta);

        UHT::FiniteImpulseResponse<double> hilbert{options, isRealTime};
        REQUIRE(hilbert.isInitialized());

        USignal::Vector<std::complex<TestType>> outputSignal(yRef.size());
        int kStart{0};
        auto nSamples = static_cast<int> (x.size());
        std::uniform_int_distribution<int> distribution{1, 500};
        for (int k = 0; k < nSamples; ++k)
        {
            int nSamplesInPacket = distribution(generator);
            auto kEnd = std::min(kStart + nSamplesInPacket, nSamples);
            nSamplesInPacket = kEnd - kStart; // Fix packet length 
            if (nSamplesInPacket <= 0){break;}
            USignal::Vector<TestType> packet(nSamplesInPacket);
            std::copy(x.data() + kStart, x.data() + kEnd,
                      packet.data());
            REQUIRE_NOTHROW(hilbert.setInput(packet));
            REQUIRE_NOTHROW(hilbert.apply());
            auto outputPacket = hilbert.getOutput();
            std::copy(outputPacket.data(),
                      outputPacket.data() + nSamplesInPacket,
                      outputSignal.data() + kStart);

            kStart = kStart + nSamplesInPacket;
            if (kStart >= nSamples){break;}
        }

        REQUIRE(yRef.size() == outputSignal.size()); 
        for (int i = 0; i < yRef.size(); ++i)
        {
            auto residual = std::abs(yRef[i] - outputSignal[i]);
            CHECK(residual < std::numeric_limits<TestType>::epsilon()*100);
        }
    }

    SECTION("Type IV Real Time")
    {   
        const bool isRealTime{true};
        const int order{301};
        const double beta{8};

        // Create reference
        auto yRef = createReference(x, order, beta, false);

        options.setOrder(order);
        options.setBeta(beta);

        UHT::FiniteImpulseResponse<double> hilbert{options, isRealTime};
        REQUIRE(hilbert.isInitialized());

        USignal::Vector<std::complex<TestType>> outputSignal(yRef.size());
        int kStart{0};
        auto nSamples = static_cast<int> (x.size());
        std::uniform_int_distribution<int> distribution{1, 500};
        for (int k = 0; k < nSamples; ++k)
        {
            int nSamplesInPacket = distribution(generator);
            auto kEnd = std::min(kStart + nSamplesInPacket, nSamples);
            nSamplesInPacket = kEnd - kStart; // Fix packet length 
            if (nSamplesInPacket <= 0){break;}
            USignal::Vector<TestType> packet(nSamplesInPacket);
            std::copy(x.data() + kStart, x.data() + kEnd,
                      packet.data());
            REQUIRE_NOTHROW(hilbert.setInput(packet));
            REQUIRE_NOTHROW(hilbert.apply());
            auto outputPacket = hilbert.getOutput();
            std::copy(outputPacket.data(),
                      outputPacket.data() + nSamplesInPacket,
                      outputSignal.data() + kStart);

            kStart = kStart + nSamplesInPacket;
            if (kStart >= nSamples){break;}
        }

        REQUIRE(yRef.size() == outputSignal.size()); 
        for (int i = 0; i < yRef.size(); ++i)
        {
            auto residual = std::abs(yRef[i] - outputSignal[i]);
            CHECK(residual < std::numeric_limits<TestType>::epsilon()*100);
        }
    }
}

