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
                   float, double)
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

    SECTION("From ZPK")
    {
        /*
        fs = 100  # sampling frequency
        om_c = 2 * np.pi * np.array([7, 13])  # corner frequencies
        zz_s, pp_s, k_s = signal.butter(4, om_c, btype='bandpass', analog=True, output='zpk')
        print("in z")
        for v in zz_s:
            print(v)
        print("in p")
        for v in pp_s:
            print(v)
        print(k_s)
        print("analog transfer function")
        b_s, a_s = signal.zpk2tf(zz_s, pp_s, k_s)
        print("b_s")
        for v in b_s:
            print(v)
        print("a_s")
        for v in a_s:
            print(v)
        """
        */
        USignal::Vector<std::complex<TestType>> zeros(
           std::vector<std::complex<TestType>> { 0 + 0i, 
                                                 0 + 0i, 
                                                 0 + 0i, 
                                                 0 + 0i });
        USignal::Vector<std::complex<TestType>> poles(
            std::vector<std::complex<TestType>> {
                -5.188311655529189 - 44.616532466296775i,
                -15.243050879954183 - 50.63131935118821i,
                -15.243050879954183 + 50.63131935118821i,
                -5.188311655529189 + 44.616532466296775i,
                -9.238513861695127 + 79.44597029196994i,
                -19.586386945718992 + 65.05814486841253i,
                -19.586386945718992 - 65.05814486841253i,
                -9.238513861695127-79.44597029196994i});
        constexpr TestType gain{2019874.9116810758};
        USignal::FilterRepresentations::ZerosPolesGain
            zpk{zeros, poles, gain};

        USignal::FilterRepresentations::InfiniteImpulseResponse<TestType>
            ba{zpk};
       
        USignal::Vector<TestType> bRef(
            std::vector<TestType> {2019874.9116810758,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0} );
        USignal::Vector<TestType> aRef(
            std::vector<TestType> {1.0,
                                   98.51252668579498,
                                   19222.50296499584,
                                   1201737.668338526,
                                   114322312.96086375,
                                   4317285838.461515,
                                   248091676925.3078,
                                   4567671318148.967,
                                   166572964959828.47});
        auto b = ba.getNumeratorFilterCoefficients();
        auto a = ba.getDenominatorFilterCoefficients();
        REQUIRE(b.size() == bRef.size());
        REQUIRE(a.size() == aRef.size());
        REQUIRE(static_cast<int> (std::max(b.size(), a.size())) - 1
                == ba.getOrder());
        for (int i = 0; i < static_cast<int> (b.size()); ++i)
        {
            CHECK(std::abs(bRef[i] - b[i]) < 1.e-1); 
        }
        CHECK(std::abs(aRef[0] - a[0]) < 1.e-5); 
        CHECK(std::abs(aRef[1] - a[1]) < 1.e-5);
        CHECK(std::abs(aRef[2] - a[2]) < 1.e-1);
        //CHECK(std::abs(aRef[8] - a[8]) < 100);
    }
}

