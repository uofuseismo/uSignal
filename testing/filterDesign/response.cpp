#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include "uSignal/filterDesign/response.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/vector.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace USignal;

TEMPLATE_TEST_CASE("CoreTest::FilterDesign::Response::AnalogResponse",
                   "[TypeName][template]",
                   double, float)
{
    Vector<TestType> bs(
        std::vector<TestType> {1611.7315706,  0.,  0.,  0.,  0.} );
    Vector<TestType> as( 
        std::vector<TestType> {1.00000000e+00,   8.57530241e+00,   1.57767906e+02,
                               7.98628595e+02,   4.76375068e+03,   7.98628595e+03,
                               1.57767906e+04,   8.57530241e+03,   1.00000000e+04} );
    Vector<std::complex<TestType>> hRef (
        std::vector<std::complex<TestType>> {
             +0.03488577413550 - 0.02699025100928i,
             +0.05348973603688 - 0.05152025023515i,
             +0.07995180799247 - 0.10404220047713i,
             +0.10278736966951 - 0.22594753568789i,
             +0.00848918600791 - 0.50531182099015i,
             -0.62270477129993 - 0.66247486038261i,
             -0.98269958366215 - 0.00383789980467i,
             -0.80405788502503 + 0.42332819774255i,
             -0.60615492623016 + 0.65515371685383i,
             -0.39658757122898 + 0.83087903105862i,
             -0.13970324532235 + 0.95304016024704i,
             +0.15497598941797 + 0.98179877285411i,
             +0.43290796318079 + 0.89967494452506i,
             +0.64293416921331 + 0.73874321858729i,
             +0.77441733458820 + 0.54927348734769i,
             +0.84574252312641 + 0.36381879901862i,
             +0.87952203818244 + 0.19215869533035i,
             +0.89095834680252 + 0.03087283491702i,
             +0.88609894228890 - 0.12857036725930i,
             +0.86230825135038 - 0.29537391207665i,
             +0.80798359282187 - 0.47566481744457i,
             +0.70276882681403 - 0.66602904861868i,
             +0.52426157690959 - 0.84353059262555i,
             +0.26865301465227 - 0.96221985462147i,
             -0.02714612573254 - 0.97737881098598i,
             -0.30297302526253 - 0.88677677485093i,
             -0.52857685878989 - 0.72830931064414i,
             -0.72551867160719 - 0.52450849935703i,
             -0.92781800750244 - 0.20388405324537i,
             -0.88469336203972 + 0.44215940372882i,
             -0.15481996440243 + 0.64414889509246i,
             +0.09562472100872 + 0.31084903769236i,
             +0.09120258029299 + 0.13942195037108i,
             +0.06291517417591 + 0.06720800519134i,
             +0.04122820755716 + 0.03451868927268i,
             +0.02685204970343 + 0.01856698156418i,
             +0.01756607670780 + 0.01032696165320i,
             +0.01156809207517 + 0.00588806178611i,
             +0.00766837172005 + 0.00342072639359i,
             +0.00511306514726 + 0.00201622740725i,
             +0.00342630435821 + 0.00120189184245i,
             +0.00230562492063 + 0.00072288480789i,
             +0.00155691392926 + 0.00043789173559i,
             +0.00105437478338 + 0.00026677845884i,
             +0.00071575337674 + 0.00016328428015i,
             +0.00048684418034 + 0.00010031589633i,
             +0.00033168506020 + 0.00006181943458i,
             +0.00022628105324 + 0.00003819156150i,
             +0.00015454510105 + 0.00002364289385i,
             +0.00010564846092 + 0.00001466105711i
        } );
    REQUIRE(hRef.size() == 50);
    constexpr double x1{-1};
    constexpr double x2{+1};
    const double dx = (x2 - x1)/static_cast<double> (hRef.size() - 1);
    USignal::Vector<TestType> frequencies(hRef.size());
    for (int i = 0; i < static_cast<int> (hRef.size()); ++i)
    {
        frequencies[i]
            = static_cast<TestType> (2*M_PI*std::pow(10.0, x1 + i*dx));
    }

    USignal::FilterRepresentations::InfiniteImpulseResponse iir{bs, as};
    auto h
         = USignal::FilterDesign::Response::computeAnalog(iir, frequencies);
    REQUIRE(h.size() == hRef.size());
    for (size_t i = 0; i < hRef.size(); ++i)
    {
        TestType residual = std::abs(hRef.at(i) - h.at(i));
        CHECK(Catch::Approx(0).margin(
                 std::numeric_limits<TestType>::epsilon()*100) == residual);
    }
}

TEMPLATE_TEST_CASE("CoreTest::FilterDesign::Response::DigitalResponse",
                   "[TypeName][template]",
                   double, float)
{
    Vector<TestType> bz( 
        std::vector<TestType> {0.056340000000000, -0.000935244000000,
                              -0.000935244000000,  0.056340000000000} );
    Vector<TestType> az( 
        std::vector<TestType> {1.000000000000000, -2.129100000000000,
                               1.783386300000000, -0.543463100000000} );
    Vector<std::complex<TestType>> hRef (
        std::vector<std::complex<TestType>> {
            +0.99987648795559 + 0.00000000000000i,
            +0.95601820836871 - 0.24635538126662i,
            +0.84182473416863 - 0.45234494448287i,
            +0.69218147536286 - 0.60548822657119i,
            +0.53073903024459 - 0.72046222089002i,
            +0.35389724345662 - 0.82273951176198i,
            +0.11450731879359 - 0.92765406032472i,
            -0.30356986746249 - 0.94904121363739i,
            -0.74535194596939 - 0.48906787904061i,
            -0.52345405359234 + 0.03883239252739i,
            -0.23022242683644 + 0.11340130414619i,
            -0.09247116985374 + 0.07419998418252i,
            -0.03163529655289 + 0.03366165116973i,
            -0.00365030620855 + 0.00475347063792i,
            +0.00945480661491 - 0.01445905891768i,
            +0.01535492883430 - 0.02690992479787i,
            +0.01759436581746 - 0.03478100394682i,
            +0.01792468519841 - 0.03954002587212i,
            +0.01725815161918 - 0.04216037245807i,
            +0.01608849455770 - 0.04329285168849i,
            +0.01468989855489 - 0.04337885254581i,
            +0.01321670410931 - 0.04272230968619i,
            +0.01175579491817 - 0.04153550307532i,
            +0.01035521965306 - 0.03996857117507i,
            +0.00904033577320 - 0.03812884046514i,
            +0.00782315309451 - 0.03609370155880i,
            +0.00670785329979 - 0.03391932843019i,
            +0.00569410163719 - 0.03164667442229i,
            +0.00477905564053 - 0.02930565474428i,
            +0.00395859088629 - 0.02691810244844i,
            +0.00322804947486 - 0.02449988289502i,
            +0.00258269445972 - 0.02206242321480i,
            +0.00201798188745 - 0.01961383022046i,
            +0.00152971946493 - 0.01715971571506i,
            +0.00111415501448 - 0.01470381187907i,
            +0.00076802196457 - 0.01224843497540i,
            +0.00048855920133 - 0.00979483895752i,
            +0.00027351634834 - 0.00734348911731i,
            +0.00012115155313 - 0.00489427799722i,
            +0.00003022628902 - 0.00244670032556i,
            -0.00000000000000 - 0.00000000000000i
        } );
    REQUIRE(hRef.size() == 41);
    const double df = (M_PI - 0)/static_cast<double> (hRef.size() - 1);
    USignal::Vector<TestType> frequencies(hRef.size());
    for (int i = 0; i < static_cast<int> (hRef.size()); ++i)
    {
        frequencies[i] = static_cast<TestType> (0 + static_cast<double> (i)*df);
    }
    USignal::FilterRepresentations::InfiniteImpulseResponse iir{bz, az};
    auto h
         = USignal::FilterDesign::Response::computeDigital(iir, frequencies);
    REQUIRE(h.size() == hRef.size());
    for (size_t i = 0; i < hRef.size(); ++i)
    {   
        TestType residual = std::abs(hRef.at(i) - h.at(i));
        CHECK(Catch::Approx(0).margin(
                 std::numeric_limits<TestType>::epsilon()*100) == residual);
    }   
}

