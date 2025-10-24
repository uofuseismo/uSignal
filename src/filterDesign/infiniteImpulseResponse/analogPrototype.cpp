#include <cmath>
#include <complex>
#include <string>
#include "uSignal/filterDesign/infiniteImpulseResponse/analogPrototype.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/vector.hpp"

USignal::FilterRepresentations::ZerosPolesGain<double> 
USignal::FilterDesign::InfiniteImpulseResponse::AnalogPrototype::butterworth(
    const int order)
{
    if (order < 0 || order > 24)
    {
        if (order < 0)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                      + " must be non-negative");
        }
        throw std::invalid_argument("order = " + std::to_string(order)
                                  + " must be less than 25");
    }
    const int nPoles{order + 1};
    const int nZeros{0};
    USignal::Vector<std::complex<double>> poles(nPoles);
    USignal::Vector<std::complex<double>> zeros(nZeros);
    const double pi_twoni = M_PI/(2.0*static_cast<double> (order + 1));
    for (int i = 0; i < nPoles; i++)
    {
        auto xm = static_cast<double> (-(order + 1) + 1 + 2*i);
        std::complex<double> arg(0.0, pi_twoni*xm);
        poles[i] =-std::exp(arg);
    }
    constexpr double gain = 1.0; //k = real(prod(-p)) which is 1
    USignal::FilterRepresentations::ZerosPolesGain<double>
        zpk{zeros, poles, gain};
    return zpk;
}

/*
ZPK AnalogPrototype::cheb1ap(const int n, const double rp)
{
    if (n < 1)
    {
        throw std::invalid_argument("order = " + std::to_string(n)
                                 + " must be positive");
    }
    if (rp <= 0)
    {
        throw std::invalid_argument("rp = " + std::to_string(rp)
                                  + " must be positive");
    }
    double rpdb = std::pow(10.0, 0.1*rp);
#ifdef DEBUG
    assert(rpdb > 1.0);
#endif
    // Set space
    size_t npoles = static_cast<size_t> (n);
    size_t nzeros = 0;
    std::vector<std::complex<double>> poles(npoles, 0); 
    std::vector<std::complex<double>> zeros(nzeros, 0);
    // Ripple factor
    double eps = std::sqrt(rpdb - 1.0);
    double xmu = 1.0/(static_cast<double> (n))*std::asinh(1.0/eps);
    // Arrange poles in an ellipse on the left half of the S-plane
    std::complex<double> zone(1, 0);
    std::complex<double> zprod = zone; // Initialize product
    double twoni = 1.0/(2.0*static_cast<double> (n));
    for (int i = 0; i < n; i++)
    {
        double xm = static_cast<double> (-n + 1 + 2*i);
        double theta = (M_PI*xm)*twoni;
        std::complex<double> arg(xmu, theta);
        poles[i] =-std::sinh(arg);
        zprod = (-poles[i])*zprod;
    }
    double k = std::real(zprod); //Take real
    if (n%2 == 0){k = k/std::sqrt(1.0 + eps*eps);}
    ZPK zpk(zeros, poles, k);
    return zpk;
}

ZPK AnalogPrototype::cheb2ap(const int n, const double rs)
{
    if (n < 1)
    {
        throw std::invalid_argument("order = " + std::to_string(n)
                                 + " must be positive");
    }
    if (rs <= 0)
    {
        throw std::invalid_argument("rs = " + std::to_string(rs)
                                  + " must be positive");
    }
    // Figure out size
    int ntarg = n;
    if (n%2 == 1){ntarg = n - 1;}
    size_t npoles = static_cast<size_t> (n);
    size_t nzeros = static_cast<size_t> (ntarg);
    std::vector<std::complex<double>> poles(npoles, 0); 
    std::vector<std::complex<double>> zeros(nzeros, 0);
    // Ripple factor check
    double rdb = pow(10.0, 0.1*rs);
#ifdef DEBUG
    assert(rdb > 1);
#endif
    // Ripple factor
    int j = 0;
    double twoni = 1.0/(2.0*static_cast<double> (n)); 
    double eps = 1.0/std::sqrt(rdb - 1.0);
    double xmu = std::asinh(1.0/eps)/static_cast<double> (n);
    // Compute zeros
    std::complex<double> zone(1, 0);
    std::complex<double> zden = zone;
    j = 0;
    // Odd
    if (n%2 == 1)
    {
        for (int i = 0; i <n/2; i++)
        {
            double xm = static_cast<double> (-n + 1 + 2*i);
            zeros[j] = std::complex<double> (0.0, 1.0/(std::sin(xm*M_PI*twoni)));
            zden = (-zeros[j])*zden;
            j = j + 1;
        }
        for (int i = 0; i <n/2; i++)
        {
            double xm = static_cast<double> (2 + 2*i);
            zeros[j] = std::complex<double> (0.0, 1.0/(sin(xm*M_PI*twoni)));
            zden = (-zeros[j])*zden;
            j = j + 1;
        }
    }
    // Even
    else
    {
        for (int i = 0; i < n; i++)
        {
            double xm = static_cast<double> (-n + 1 + 2*i);
            zeros[j] = std::complex<double> (0.0, 1.0/(std::sin(xm*M_PI*twoni)));
            zden = (-zeros[j])*zden;
            j = j + 1;
        }
    }
    // Poles around unit circle like butterworth; then warp into cheby II
    double sinhmu = 0.5*(std::exp(xmu) - std::exp(-xmu));
    double coshmu = 0.5*(std::exp(xmu) + std::exp(-xmu));
    std::complex<double> znum = zone;
    for (int i = 0; i < n; i++)
    {
       //arg = 0.0 + (M_PI*((double)(-n + 1 + 2*i))*twoni)*_Complex_I;
       double temp = M_PI*static_cast<double> (-n + 1 + 2*i)*twoni;
       std::complex<double> arg = std::complex<double> (0.0, temp);
       std::complex<double> polesi =-std::exp(arg);
       // Warp into cheby II
       poles[i] = std::complex<double> (sinhmu*std::real(polesi),
                                        coshmu*std::imag(polesi));
       poles[i] = zone/poles[i]; //p = 1/p
       znum = (-poles[i])*znum;
    }
    double k = std::real(znum/zden);
    ZPK zpk(zeros, poles, k);
    return zpk;
}
*/
