#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include "uSignal/vector.hpp"
#include "utilities/math/polynomial.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEMPLATE_TEST_CASE("CoreTest::Utilities::Math::Polynomial",
                   "[TypeName][template]",
                   float, double)//, std::complex<std::double>)
{

    USignal::Vector<TestType> x({5, 7, 9, -2});
    SECTION("constant")
    {
        USignal::Vector<TestType> p0(1, 3); // 1 element filled with 3
        USignal::Vector<TestType> yRef({3, 3, 3, 3});
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p0, x);
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y[i] == Catch::Approx(yRef[i]));
        }
    }
    SECTION("linear") 
    {
        USignal::Vector<TestType> p1(std::vector<TestType> {3, 2});
        USignal::Vector<TestType> yRef({17, 23,  29, -4});
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p1, x);
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y[i] == Catch::Approx(yRef[i]));
        }
    }
    SECTION("quadratic")
    {
        USignal::Vector<TestType> p2({3, 2, 1});
        USignal::Vector<TestType> yRef({86, 162, 262, 9});
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p2, x);
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y[i] == Catch::Approx(yRef[i]));
        }
    }
    SECTION("cubic")
    {
        USignal::Vector<TestType> p3({3, 2, 1, -2});
        USignal::Vector<TestType> yRef({428, 1132, 2356, -20});
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p3, x);
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y[i] == Catch::Approx(yRef[i]));
        }
    } 
    SECTION("quartic")
    {
        USignal::Vector<TestType> p4({3, 2, 1, -2, -3});
        USignal::Vector<TestType> yRef({2137, 7921, 21201,  37});
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p4, x);
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y[i] == Catch::Approx(yRef[i]));
        }
    }

    USignal::Vector<std::complex<TestType>>  zx;
    zx.push_back(std::complex<TestType> (5, 0) );
    zx.push_back(std::complex<TestType> (0, 7) );
    zx.push_back(std::complex<TestType> (2, 9) );
    zx.push_back(std::complex<TestType> (4,-2) );
    SECTION("complex constant")
    {
        USignal::Vector<std::complex<TestType>> p0;
        p0.push_back( std::complex<TestType> (3,0));

        USignal::Vector<std::complex<TestType>>
            yRef( { 3 + 0i, 3 + 0i, 3 + 0i, 3 + 0i } );
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p0, zx);
        constexpr TestType zero{0};
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(std::abs(y[i] - yRef[i]) == Catch::Approx(zero));
        }
    }
    SECTION("complex linear")
    {
        USignal::Vector<std::complex<TestType>>
            p1( std::vector<std::complex<TestType>> {3 + 0i, 0 + 2i} );
        USignal::Vector<std::complex<TestType>>
            yRef( {15 + 2i, 0 + 23i, 6 + 29i, 12 - 4i} );
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p1, zx);
        constexpr TestType zero{0};
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(std::abs(y[i] - yRef[i]) == Catch::Approx(zero));
        }
    }
    SECTION("complex quadratic")
    {
        USignal::Vector<std::complex<TestType>>
            p2( {3 + 0i, 0 + 2i, 1 + 0i} );
        USignal::Vector<std::complex<TestType>>
            yRef( {76 + 10i, -160 + 0i, -248 + 112i, 41 - 40i} );
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p2, zx);
        constexpr TestType zero{0};
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(std::abs(y[i] - yRef[i]) == Catch::Approx(zero));
        }
    }
    SECTION("complex cubic")
    {
        USignal::Vector<std::complex<TestType>>
            p3( {3 + 0i, 0 + 2i, 1 + 0i, 1 - 2i} );
        USignal::Vector<std::complex<TestType>> 
            yRef( {381 + 48i, 1 - 1122i, -1503 - 2010i, 85 - 244i} );
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p3, zx);
        constexpr TestType zero{0};
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(std::abs(y[i] - yRef[i]) == Catch::Approx(zero));
        }
    }
    SECTION("complex quartic")
    {
        USignal::Vector<std::complex<TestType>>
             p4( {3 + 0i, 0 + 2i, 1 + 0i, 1 - 2i, -3 + 0i } );
        USignal::Vector<std::complex<TestType>> 
            yRef( {1902 + 240i, 7851 + 7i, 15081 - 17547i, -151 - 1146i} );
        auto y = USignal::Utilities::Math::Polynomial::evaluate(p4, zx);
        constexpr TestType zero{0};
        CHECK(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(std::abs(y[i] - yRef[i]) == Catch::Approx(zero));
        }
    }
}

TEMPLATE_TEST_CASE("CoreTest::Utilities::Math::Polynomial::Roots",
                   "[TypeName][template]",
                   double, float)
{
    
    constexpr TestType zero{0};
    const TestType tolerance{50*std::numeric_limits<TestType>::epsilon()};

    SECTION("Real")
    {
        USignal::Vector<TestType> coefficients( std::vector<TestType> {1, -6, -72, -27} );
        auto roots = USignal::Utilities::Math::Polynomial::computeRoots(coefficients);
  
        USignal::Vector<std::complex<TestType>> rootsRef(
               {12.122893784632392  + 0i,
                -5.7345099422250705 + 0i,
                -0.3883838424073199 + 0i});
        CHECK(roots.size() == rootsRef.size());
        for (int i = 0; i < static_cast<int> (rootsRef.size()); ++i)
        {
            CHECK(std::abs(roots[i] - rootsRef[i]) ==
                  Catch::Approx(zero).margin(tolerance));
        }
    }

    SECTION("Complex")
    {
        USignal::Vector<TestType> coefficients( std::vector<TestType> {1, -6, 72, 27} );
        auto roots = USignal::Utilities::Math::Polynomial::computeRoots(coefficients);

        USignal::Vector<std::complex<TestType>> rootsRef(
           {3.1816664666582546 + 8.011804223473874i,
            3.1816664666582546 - 8.011804223473874i,
           -0.3633329333165073 + 0i});
        CHECK(roots.size() == rootsRef.size());
        for (int i = 0; i < static_cast<int> (rootsRef.size()); ++i)
        {
            bool found{false};
            for (int j = 0; j < static_cast<int> (rootsRef.size()); ++j)
            {
                if (std::abs(roots[i] - rootsRef[j]) < tolerance)
                {
                    found = true;
                    break;
                }
            }
            CHECK(found); 
        }
    }
}

/*
TEMPLATE_TEST_CASE("CoreTest::Utilities::Math::ComplexPolynomial",
                   "[TypeName][template]",
                   std::complex<double>, std::complex<float> )
{

    std::vector<std::complex<double>> zx; 
    zx.push_back(std::complex<double> (5,0) );
    zx.push_back(std::complex<double> (0,7) );
    zx.push_back(std::complex<double> (2,9) );
    zx.push_back(std::complex<double> (4,-2) );

}
*/

TEMPLATE_TEST_CASE("CoreTest::Utilities::Math::ExpandPolynomial",
                   "[TypeName][template]",
                   double, float)
{
    const TestType tolerance{50*std::numeric_limits<TestType>::epsilon()};
    constexpr TestType zero{0};
    constexpr TestType one{1};
    SECTION("Edge Case Order 0")
    {
        USignal::Vector<std::complex<TestType>> roots;
        auto coefficients
            = USignal::Utilities::Math::Polynomial::expand(roots);
        REQUIRE(coefficients.size() == 1);
        CHECK(std::abs(coefficients.at(0)) ==
              Catch::Approx(one).margin(tolerance));
    }
    SECTION("Edge Case Order 1")
    {
        USignal::Vector<std::complex<TestType>> roots;
        roots.push_back( -1 + 2i );
        auto coefficients
            = USignal::Utilities::Math::Polynomial::expand(roots);
        REQUIRE(coefficients.size() == 2); 
        CHECK(std::abs(coefficients.at(0)) ==
              Catch::Approx(one).margin(tolerance));
        // Coefficient 1 =-root -> 1 + root = 0
        CHECK(std::abs(coefficients.at(1) + roots.at(0)) ==
              Catch::Approx(zero).margin(tolerance));
    }
    SECTION("Complex Test")
    {
        USignal::Vector<std::complex<TestType>> reference;
        reference.push_back( 1    + 0i);
        reference.push_back(-5.3  + 1i);
        reference.push_back( 6.3  - 5.3i);
        reference.push_back(-10.6 + 4.3i);
        reference.push_back( 8.6  + 0i);
        USignal::Vector<std::complex<TestType>> roots;
        roots.push_back(1 + 0i);
        roots.push_back(0 + 1i);
        roots.push_back(4.3 + 0i);
        roots.push_back(0 - 2i);
        auto coefficients
            = USignal::Utilities::Math::Polynomial::expand(roots);
        REQUIRE(coefficients.size() == reference.size());
        for (int i = 0; i < static_cast<int> (reference.size()); ++i)
        {
            CHECK(std::abs(reference.at(i) - coefficients.at(i)) ==
                  Catch::Approx(zero).margin(tolerance));
        }
    }
    SECTION("Real Test")
    {
        USignal::Vector<TestType> reference( { 1., -4.3, -3., 14.9, -8.6} );
        USignal::Vector<std::complex<TestType>> roots;
        roots.push_back(   1 + 0i);
        roots.push_back(   1 + 0i);
        roots.push_back( 4.3 + 0i);
        roots.push_back(-2.0 + 0i);
        auto coefficients
            = USignal::Utilities::Math::Polynomial::expandToRealCoefficients(roots);
        REQUIRE(coefficients.size() == reference.size());
        for (int i = 0; i < static_cast<int> (reference.size()); ++i)
        {
            CHECK(coefficients.at(i) ==
                  Catch::Approx(reference.at(i)).margin(tolerance));
        }
    }
}
