#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include "uSignal/vector.hpp"
#include "uSignal/waterLevelTrigger.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;

TEMPLATE_TEST_CASE("CoreTest::WaterLevelTrigger",
                   "[TypeName][template]", float, double)
{
    SECTION("Options")
    {
        constexpr double onThreshold{0.5};
        constexpr double offThreshold{0.4};
        USignal::WaterLevelTriggerOptions
            options{std::pair {onThreshold, offThreshold}};
        auto onAndOffThreshold = options.getOnAndOffThreshold();
        REQUIRE(std::abs(onAndOffThreshold.first  - onThreshold) < 1.e-14);
        REQUIRE(std::abs(onAndOffThreshold.second - offThreshold) < 1.e-14);
    }
    constexpr int nSamples{1000};
    USignal::Vector<TestType> offOnSignal{nSamples, 0};
    USignal::Vector<TestType> onOffSignal{nSamples, 0};
    USignal::Vector<int> offOnReferenceSignal{nSamples, 0};
    USignal::Vector<int> onOffReferenceSignal{nSamples, 0};
    bool isOff{false};
    for (int i = 0; i < nSamples; ++i)
    {
        if (isOff)
        {
            offOnSignal[i] = 1; offOnReferenceSignal[i] = 1;
            onOffSignal[i] =-1; onOffReferenceSignal[i] = 0;
        }
        else
        {
            offOnSignal[i] =-1; offOnReferenceSignal[i] = 0;
            onOffSignal[i] = 1; onOffReferenceSignal[i] = 1;
        }
        if (i%10 == 0){isOff = !isOff;}
    }
    SECTION("Post processing start off")
    {
        constexpr bool isRealTime{false};
        constexpr double onThreshold{0.5};
        constexpr double offThreshold{0.4};
        USignal::WaterLevelTriggerOptions
            options{std::pair {onThreshold, offThreshold}};
        USignal::WaterLevelTrigger<TestType> trigger{options, isRealTime};
        REQUIRE(trigger.isInitialized());
        trigger.setInput(offOnSignal); 
        trigger.apply();
        auto y = trigger.getOutput();
        REQUIRE(y.size() == offOnReferenceSignal.size()); 
        for (int i = 0; i < offOnReferenceSignal.size(); ++i)
        {
            CHECK(y.at(i) == offOnReferenceSignal.at(i));
        }
    }
    SECTION("Post processing start on")
    {
        constexpr bool isRealTime{false};
        constexpr double onThreshold{0.5};
        constexpr double offThreshold{0.4};
        USignal::WaterLevelTriggerOptions
            options{std::pair {onThreshold, offThreshold}};
        USignal::WaterLevelTrigger<TestType> trigger{options, isRealTime};
        REQUIRE(trigger.isInitialized());
        trigger.setInput(onOffSignal); 
        trigger.apply();
        auto y = trigger.getOutput();
        REQUIRE(y.size() == onOffReferenceSignal.size()); 
        for (int i = 0; i < onOffReferenceSignal.size(); ++i)
        {
            CHECK(y.at(i) == onOffReferenceSignal.at(i));
        }
    } 
    SECTION("Real time start off")
    {
        constexpr bool isRealTime{false};
        constexpr double onThreshold{0.5};
        constexpr double offThreshold{0.4};
        USignal::WaterLevelTriggerOptions
            options{std::pair {onThreshold, offThreshold}};
        USignal::WaterLevelTrigger<TestType> trigger{options, isRealTime};
        REQUIRE(trigger.isInitialized());

        // Apply the signal over the filter length
        std::mt19937 generator(83823);
        std::uniform_int_distribution<int> distribution{1, 10};

        USignal::Vector<int> outputSignal(offOnReferenceSignal.size());
        int kStart{0};
        auto nSamples = static_cast<int> (offOnSignal.size());
        for (int k = 0; k < nSamples; ++k)
        {
            int nSamplesInPacket = distribution(generator);
            auto kEnd = std::min(kStart + nSamplesInPacket, nSamples);
            nSamplesInPacket = kEnd - kStart; // Fix packet length 
            if (nSamplesInPacket <= 0){break;}
            USignal::Vector<TestType> packet(nSamplesInPacket);
            std::copy(offOnSignal.data() + kStart, offOnSignal.data() + kEnd,
                      packet.data());
            REQUIRE_NOTHROW(trigger.setInput(packet));
            REQUIRE_NOTHROW(trigger.apply());
            auto outputPacket = trigger.getOutput();
            std::copy(outputPacket.data(),
                      outputPacket.data() + nSamplesInPacket,
                      outputSignal.data() + kStart);

            kStart = kStart + nSamplesInPacket;
            if (kStart >= nSamples){break;}
        }
        for (int i = 0; i < offOnReferenceSignal.size(); ++i)
        {
            CHECK(outputSignal.at(i) == offOnReferenceSignal.at(i));
        }
    }

    SECTION("Real time start on")
    {
        constexpr bool isRealTime{false};
        constexpr double onThreshold{0.5};
        constexpr double offThreshold{0.4};
        USignal::WaterLevelTriggerOptions
            options{std::pair {onThreshold, offThreshold}};
        USignal::WaterLevelTrigger<TestType> trigger{options, isRealTime};
        REQUIRE(trigger.isInitialized());

        // Apply the signal over the filter length
        std::mt19937 generator(83823);
        std::uniform_int_distribution<int> distribution{1, 10};

        USignal::Vector<int> outputSignal(onOffReferenceSignal.size());
        int kStart{0};
        auto nSamples = static_cast<int> (onOffSignal.size());
        for (int k = 0; k < nSamples; ++k)
        {
            int nSamplesInPacket = distribution(generator);
            auto kEnd = std::min(kStart + nSamplesInPacket, nSamples);
            nSamplesInPacket = kEnd - kStart; // Fix packet length 
            if (nSamplesInPacket <= 0){break;}
            USignal::Vector<TestType> packet(nSamplesInPacket);
            std::copy(onOffSignal.data() + kStart, onOffSignal.data() + kEnd,
                      packet.data());
            REQUIRE_NOTHROW(trigger.setInput(packet));
            REQUIRE_NOTHROW(trigger.apply());
            auto outputPacket = trigger.getOutput();
            std::copy(outputPacket.data(),
                      outputPacket.data() + nSamplesInPacket,
                      outputSignal.data() + kStart);

            kStart = kStart + nSamplesInPacket;
            if (kStart >= nSamples){break;}
        }
        for (int i = 0; i < onOffReferenceSignal.size(); ++i)
        {
            CHECK(outputSignal.at(i) == onOffReferenceSignal.at(i));
        }
    }


}

