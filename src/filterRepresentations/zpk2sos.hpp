#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <cmath>
#include <complex>
namespace
{

int getNextWorstRealPoleIndex(const std::vector<bool> &pmask,
                                     const std::vector<bool> &isreal_pole,
                                     const std::vector<std::complex<double>> &p);
int getWorstPoleIndex(const std::vector<bool> &pmask,
                      const std::vector<std::complex<double>> &p);
int isreal_p_sum(const std::vector<bool> &pmask,
                        const std::vector<bool> &isreal_pole);
int nearest_real_complex_idx(const std::vector<bool> &zmask,
                             const std::vector<bool> &isreal_zero,
                             const std::vector<std::complex<double>> &z,
                             const std::complex<double> p1,
                             const bool lreal);
int nearest_realOrComplex_idx(const std::vector<bool> &zmask,
                              const std::vector<std::complex<double>> &z, 
                              const std::complex<double> p1);
int cmplxreal(const std::vector<std::complex<double>> &z, 
              std::vector<std::complex<double>> &p,
              std::vector<bool> &isreal,
              const double *tolIn);

std::pair
<
  std::vector<std::array<double, 3>>,
  std::vector<std::array<double, 3>>
>
zpk2sos(const std::vector<std::complex<double>> &zin,
        const std::vector<std::complex<double>> &pin,
        const double gainIn,
        const bool pairingStrategyNearest)
{    
    const double k{gainIn}; //zpkWork.getGain();
    int nzeros = static_cast<int> (zin.size());
    int npoles = static_cast<int> (pin.size());
    if (npoles <= 0 || nzeros <= 0)
    {
        std::cerr << "sos will simply scale data by " << k << std::endl;
        std::vector<std::array<double, 3>> bs( {{k, 0, 0}} );
        std::vector<std::array<double, 3>> as( {{1, 0, 0}} );
        //sos = RTSeis::FilterRepresentations::SOS(1, bs, as);
        return std::pair {bs, as}; //sos;
    }
    // Get sizes and handles to zeros and poles
    int np = npoles + std::max(nzeros - npoles, 0);
    int nz = nzeros + std::max(npoles - nzeros, 0);
    int nSections = std::max(np, nz + 1)/2;
#ifndef NDEBUG
    assert(nSections > 0);
#endif
    //std::vector<std::complex<double>> zin = zpkWork.getZeros();
    //std::vector<std::complex<double>> pin = zpkWork.getPoles();
    // Base case
    if (nSections == 1)
    {
        // (z - z_1)*(z - 0) = z^2 - z*z_1 + 0
        std::array<double, 3> bs({1, 0, 0});
        if (nzeros == 1)
        {
            std::cerr << "nzeros = 1 not yet tested" << std::endl;
            if (std::abs(std::imag(zin[0])) > 1.e-15)
            {
                std::cerr << "Taking real part of zin" << std::endl;
            }
            bs[1] =-std::real(zin[0]);
        }
        else if (nzeros == 2)
        {
            bs[1] =-2.0*std::real(zin[0]);
            bs[2] = std::pow(std::abs(zin[0]), 2);
        }
        // Introduce gain
        bs[0] = bs[0]*k;
        bs[1] = bs[1]*k;
        bs[2] = bs[2]*k;
        // (p - p_1)*(p - 0) = p^2 - p*p_1 + 0
        std::array<double, 3> as({1, 0, 0});
        if (npoles == 1)
        {
            std::cerr << "npoles = 1 not yet tested" << std::endl;
            if (std::abs(std::imag(pin[0])) > 1.e-15)
            {
                std::cerr << "Taking real part of pin" << std::endl;
            }
            as[1] =-std::real(pin[0]);
        }
        // Conjugate pairs: (p - p_1)*(p - conj(p_1)) 
        else if (npoles == 2)
        {
            as[1] =-2.0*std::real(pin[0]);
            as[2] = std::pow(std::abs(pin[0]), 2);
        }
        std::vector<std::array<double, 3>> bsOut;
        bsOut.push_back(std::move(bs));
        std::vector<std::array<double, 3>> asOut;
        asOut.push_back(std::move(as));
        return std::pair {bsOut, asOut};
    }
    if (np%2 == 1 && pairingStrategyNearest) //pairing == SOSPairing::NEAREST)
    {
        np = np + 1;
        nz = nz + 1;
    }
    if (np != nz || np < 1 || nz < 1)
    {
        if (np != nz)
        {
            throw std::invalid_argument("Inconsistent size: nz = "
                                      + std::to_string(nz) + " != np = "
                                      + std::to_string(np));
        }
        if (nz < 1){throw std::invalid_argument("No zeros in use");}
        throw std::invalid_argument("No poles in use");
    }
    std::vector<std::complex<double>> z;
    std::vector<bool> isreal_zero;
#ifndef NDEBUG
    assert(::cmplxreal(zin, z, isreal_zero, nullptr) == 0);
#else
    ::cmplxreal(zin, z, isreal_zero, nullptr);
#endif
    std::vector<std::complex<double>> p;
    std::vector<bool> isreal_pole;
#ifndef NDEBUG
    assert(::cmplxreal(pin, p, isreal_pole, nullptr) == 0);
#else
    ::cmplxreal(pin, p, isreal_pole, nullptr); 
#endif
    std::vector<bool> zmask(z.size(), false);
    std::vector<bool> pmask(p.size(), false);
    //std::vector<BA> basAll;
    std::vector<std::array<double, 3>> bsAll;
    std::vector<std::array<double, 3>> asAll;
    for (int is = 0; is < nSections; ++is)
    {
        std::complex<double> z1 = std::complex<double> (0, 0);
        // Select the `worst' pole
        int p1_idx = getWorstPoleIndex(pmask, p);
#ifndef NDEBUG
        assert(p1_idx !=-1);
#endif
        std::complex<double> p1 = p[p1_idx];
        pmask[p1_idx] = true;
        int psum = isreal_p_sum(pmask, isreal_pole);
        // Initialize complex conjugates for real case
        std::complex<double> p2 = std::complex<double> (0, 0);
        std::complex<double> z2 = std::complex<double> (0, 0);
        // Pair that pole with a zero
        if (isreal_pole[p1_idx] && psum == 0)
        {
            // Special case to set a first order section
            int z1_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                  z, p1, true);
#ifndef NDEBUG
            assert(z1_idx !=-1);
#endif
            z1 = z[z1_idx];
            zmask[z1_idx] = true;
        }
        else
        {
            int zsum = isreal_p_sum(zmask, isreal_zero);
            int z1_idx =-1;
            if (!isreal_pole[p1_idx] && zsum == 1)
            {
                // Special case to ensure choose a complex zero to pair
                // with so later (setting up a first-order section)
                z1_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                  z, p1, false);
#ifndef NDEBUG
                assert(z1_idx !=-1);
                assert(isreal_zero[z1_idx]);
#endif
            }
            else
            {
                // Pair that pole with the closest zero (real or complex)
                z1_idx = nearest_realOrComplex_idx(zmask, z, p1);
#ifndef NDEBUG
                assert(z1_idx !=-1);
#endif
            }
            z1 = z[z1_idx];
            zmask[z1_idx] = true;
            // Now that we have p1 and z1, figure out what p2 and z2 need to be
            if (!isreal_pole[p1_idx])
            {
                // Complex pole, complex zero -> complex conjugates
                if (!isreal_zero[z1_idx])
                {
                    p2 = std::conj(p1);
                    z2 = std::conj(z1);
                }
                // Complex pole and real zero
                else
                {
                    p2 = std::conj(p1);
                    int z2_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                          z, p1, true);
#ifndef NDEBUG
                    assert(z2_idx !=-1);
                    assert(isreal_zero[z2_idx]);
#endif
                    z2 = z[z2_idx];
                    zmask[z2_idx] = true;
                }
            }
            // Real pole
            else
            {
                // Real pole, complex zero
                if (!isreal_zero[z1_idx])
                {
                    z2 = conj(z1);
                    int p2_idx = nearest_real_complex_idx(pmask, isreal_pole,
                                                          p, z1, true);
#ifndef NDEBUG
                    assert(p2_idx !=-1);
                    assert(isreal_pole[p2_idx]);
#endif
                    p2 = p[p2_idx];
                    pmask[p2_idx] = true;
                }
                // Real pole, real zero
                else
                {
                    // Pick the next `worst' pole to use
                    int p2_idx = getNextWorstRealPoleIndex(pmask,
                                                           isreal_pole, p);
#ifndef NDEBUG
                    assert(p2_idx !=-1);
                    assert(isreal_pole[p2_idx]);
#endif
                    p2 = p[p2_idx];
                    // Find a real zero to match the added pole
                    int z2_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                          z, p2, true);
#ifndef NDEBUG
                    assert(z2_idx !=-1);
                    assert(isreal_zero[z2_idx]);
#endif
                    z2 = z[z2_idx];
                    zmask[z2_idx] = true; 
                    pmask[p2_idx] = true;
                }
            }
        }
        USignal::Vector<std::complex<double>>
            ztemp( std::vector<std::complex<double>> {z1, z2} );
        USignal::Vector<std::complex<double>>
            ptemp( std::vector<std::complex<double>> {p1, p2} );
        double ktemp{1};
        if (is == nSections - 1){ktemp = k;}
        USignal::FilterRepresentations::ZerosPolesGain<double>
            zpk{ztemp, ptemp, ktemp};
        USignal::FilterRepresentations::InfiniteImpulseResponse<double>
            ba{zpk}; 
        auto bWork = ba.getNumeratorFilterCoefficientsReference();
        auto aWork = ba.getDenominatorFilterCoefficientsReference();
        std::array<double, 3> bArray{ bWork[0], bWork[1], bWork.at(2) };
        std::array<double, 3> aArray{ aWork[0], aWork[1], aWork.at(2) };
        bsAll.push_back(std::move(bArray));
        asAll.push_back(std::move(aArray));
/*
        RTSeis::FilterRepresentations::ZPK zpkTemp
             = RTSeis::FilterRepresentations::ZPK(ztemp, ptemp, ktemp);
        BA baTemp;
        baTemp = zpk2tf(zpkTemp);
        // Save the transfer function of the secon dorder section
        basAll.push_back(baTemp);
*/
    } // Loop on sections
    // Reality check
    for (size_t ip=0; ip<pmask.size(); ip++)
    {
        if (!pmask[ip])
        {
            std::cerr << "Failed to find pole " << ip << " nSections = "
                      << nSections << std::endl;
#ifndef NDEBUG
            assert(false);
#else
            throw std::runtime_error("Failed to find pole " 
                                   + std::to_string(ip) + " nSections = "
                                   + std::to_string(nSections)); 
#endif
            //return sos;
        }
    }
    for (size_t iz=0; iz<zmask.size(); iz++)
    {
        if (!zmask[iz])
        {
            std::cerr << "Failed to find zero " << iz << std::endl;
#ifndef NDEBUG
            assert(false);
#else
            throw std::runtime_error("Failed to find zero " 
                                   + std::to_string(iz) + " nSections = "
                                   + std::to_string(nSections));
#endif
            //return sos;
        }
    }
    // Construct the cascading transfer function with the `worst' pole last
    std::reverse(bsAll.begin(), bsAll.end());
    std::reverse(asAll.begin(), asAll.end());
/*
    int js = nSections;
    std::vector<double> bs;
    std::vector<double> as;
    for (int is=0; is<nSections; is++)
    {
        js = js - 1;
        std::vector<double> bsTemp = basAll[js].getNumeratorCoefficients(); 
        std::vector<double> asTemp = basAll[js].getDenominatorCoefficients();
        for (size_t k=0; k<3; k++)
        {
            bs.push_back(bsTemp[k]);
            as.push_back(asTemp[k]);
        }
    }
*/
    // Finally pack the second order sections
    return std::pair {bsAll, asAll};
}

int nearest_real_complex_idx(const std::vector<bool> &zmask,
                             const std::vector<bool> &isreal_zero,
                             const std::vector<std::complex<double>> &z,
                             const std::complex<double> p1,
                             const bool lreal)
{
    double difMin = std::numeric_limits<double>::max(); //DBL_MAX;
    int indx =-1;
    // Look for the closest real zero to the pole p1
    if (lreal)
    {
        for (size_t i=0; i<z.size(); i++)
        {
            if (isreal_zero[i] && !zmask[i])
            {
                double dif = std::abs(z[i] - p1);
                if (dif < difMin)
                {
                    difMin = dif;
                    indx = i;
                }
            }
        }
    }
    // Look for the closest complex zero to the pole p1
    else
    {
        for (size_t i=0; i<z.size(); i++)
        {
            if (!isreal_zero[i] && !zmask[i])
            {
                double dif = std::abs(z[i] - p1);
                if (dif < difMin)
                {
                    difMin = dif;
                    indx = i;
                }
            }
        }
    }
    return indx;
}

int getNextWorstRealPoleIndex(const std::vector<bool> &pmask,
                              const std::vector<bool> &isreal_pole,
                              const std::vector<std::complex<double>> &p)
{
    const std::complex<double> zone(1, 0);
    double difMin = std::numeric_limits<double>::max(); //DBL_MAX;
    int indx =-1;
    for (size_t i=0; i<p.size(); i++)
    {
        if (!pmask[i] && isreal_pole[i])
        {
            double dif = std::abs(p[i] - zone);
            if (dif < difMin)
            {
                dif = difMin;
                indx = i;
            }
        }
    }
    return indx;
}

int getWorstPoleIndex(const std::vector<bool> &pmask,
                      const std::vector<std::complex<double>> &p)
{
    double difMin = std::numeric_limits<double>::max(); //DBL_MAX;
    int indx =-1;
    for (size_t i=0; i<pmask.size(); i++)
    {
        if (!pmask[i])
        {
            double dif = std::abs(1.0 - std::abs(p[i]));
            if (dif < difMin)
            {
                difMin = dif;
                indx = i;
            }
        }
    }
    return indx;
}

int isreal_p_sum(const std::vector<bool> &pmask,
                 const std::vector<bool> &isreal_pole)
{
    int xsum = 0;
    for (size_t i=0; i<isreal_pole.size(); i++)
    {
        if (!pmask[i] && isreal_pole[i]){xsum = xsum + 1;}
    }
    return xsum;
}


int nearest_realOrComplex_idx(const std::vector<bool> &zmask,
                              const std::vector<std::complex<double>> &z,
                              const std::complex<double> p1)
{
    double difMin = std::numeric_limits<double>::max(); //DBL_MAX;
    int indx =-1;
    for (size_t i=0; i<z.size(); i++)
    {
        if (!zmask[i])
        {
            double dif = std::abs(p1 - z[i]);
            if (dif < difMin)
            {
                difMin = dif;
                indx = i;
            }
        }
    }
    return indx;
}

int cmplxreal(const std::vector<std::complex<double>> &z,
              std::vector<std::complex<double>> &p,
              std::vector<bool> &isreal,
              const double *tolIn)
{
    p.resize(0);
    isreal.resize(0);
    if (z.size() < 1)
    {
        std::cerr << "No elements in z" << std::endl;
        return -1;
    }
    // Set the tolerance
    double tol = std::numeric_limits<double>::epsilon()*100; //100.0*DBL_EPSILON;
    if (tolIn != nullptr){tol = *tolIn;}
    // Get the real poles and complex poles 
    std::vector<double> reals;
    std::vector<std::complex<double>> cmplxs;
    for (size_t i=0; i<z.size(); i++)
    {
        if (std::imag(z[i]) == 0)
        {
            reals.push_back(std::real(z[i]));
        }
        else
        {
            cmplxs.push_back(z[i]);
        }
    }
    size_t nreal  = reals.size();
    size_t ncmplx = cmplxs.size();
    if (nreal + ncmplx != z.size())
    {
        std::cerr << "Algorithmic error" << std::endl;
        return -1;
    }
    // Sort the reals into ascending order and save them 
    nreal = 0;
    if (reals.size() > 0)
    {
        std::sort(reals.begin(), reals.end());
        for (size_t i=0; i<reals.size(); i++)
        {
            p.push_back(std::complex<double> (reals[i], 0));
            isreal.push_back(true);
            nreal = nreal + 1;
        }
    }
    // Pair the complex conjugate poles up and save them
    size_t nconj = 0;
    if (cmplxs.size() > 0)
    {
        // Sort complex numbers based on magnitudes
        std::sort(cmplxs.begin(), cmplxs.end(),
                  [](std::complex<double> a,  
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
        // Verify everyone has a buddy (conjugate pair)
        std::vector<int> lskip(cmplxs.size(), 0);
        for (size_t i=0; i<cmplxs.size(); i++)
        {
            if (lskip[i] == 1){continue;}
            // Find it's complex conjugate
            for (size_t j=i+1; j<cmplxs.size(); j++)
            {
                if (lskip[j] == 1){continue;}
                if (std::abs(cmplxs[i] - std::conj(cmplxs[j])) < tol)
                {
                    // Balance poles for numerical accuracy
                    std::complex<double> zhalf
                         = (cmplxs[i] + std::conj(cmplxs[j]))/2.0;
                    // Retain the `positive' pole for consistency 
                    if (std::imag(zhalf) >= 0)
                    {
                        p.push_back(zhalf); 
                        //p.push_back(std::conj(zhalf));
                    }
                    else
                    {
                        p.push_back(std::conj(zhalf)); 
                        //p.push_back(zhalf);
                    }
                    nconj = nconj + 1;
                    isreal.push_back(false);
                    //isreal.push_back(false);
                    lskip[i] = 1;
                    lskip[j] = 1;
                    break;
                }
            }
        }
        // I'm not dealing with unpaired
        int npaired = std::accumulate(lskip.begin(), lskip.end(), 0);
        if (npaired != static_cast<int> (cmplxs.size()))
        {
            std::cerr << "All poles need complex conjugate pairs: "
                      << npaired << " " << cmplxs.size() << std::endl;
            std::cerr << "Input" << std::endl;
            for (size_t m=0; m<z.size(); m++)
            {
                std::cerr << z[m] << std::endl;//fprintf(stderr, "%e + %+ei\n", z[m].real(), z[m].imag());
            }
            std::cerr << "Found complexes" << std::endl;
            for (size_t m=0; m<cmplxs.size(); m++) 
            {
                std::cerr << cmplxs[m] << std::endl;// fprintf(stderr, "%e + %+ei\n", cmplxs[m].real(), cmplxs[m].imag());
            }
            p.resize(0);
            std::runtime_error("Failed to match conjugate pairs"); //getchar();
            return -1;
        }
    }
    if (2*nconj + nreal != z.size())
    {
        std::cerr << "Algorithmic failure" << std::endl;
        return -1;
    }
/*
printf("Input\n");
for (size_t i=0; i<z.size(); i++)
{
 printf("%+lf%+lfi\n", std::real(z[i]), std::imag(z[i]));
}
printf("Output\n");
for (size_t i=0; i<p.size(); i++)
{
 printf("%+lf%+lfi\n", std::real(p[i]), std::imag(p[i]));
}
*/
    return 0;
}
}
