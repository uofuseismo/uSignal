#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include "uSignal/filterRepresentations/secondOrderSections.hpp"
#include "uSignal/filterRepresentations/zerosPolesGain.hpp"
#include "uSignal/filterRepresentations/infiniteImpulseResponse.hpp"
#include "uSignal/vector.hpp"

using namespace USignal::FilterRepresentations;

namespace
{
template<typename T>
USignal::Vector<T> reshape(const std::vector<std::array<T, 3>> &x,
                           const bool isNumerator)
{
    if (x.empty())
    {
        if (isNumerator)
        {
            throw std::invalid_argument("Numerator coefficient array empty");
        }
        else
        {
            throw std::invalid_argument("Denominator coefficient array empty");
        }
    }
    USignal::Vector<T> result;
    result.reserve(3*x.size());
    for (int section = 0; section < static_cast<int> (x.size()); ++section)
    {
        auto coefficient0 = x[section][0];
        if (coefficient0 == 0)
        {
            if (isNumerator)
            {
                throw std::invalid_argument(
                   "Numerator Leading coefficient zero for section "
                  + std::to_string(section + 1));
            }
            else
            {
                throw std::invalid_argument(
                   "Denominator leading coefficient zero for section "
                  + std::to_string(section + 1));
            }
        }
        result.push_back(coefficient0);
        result.push_back(x[section][1]);
        result.push_back(x[section][2]);
    }
    return result;
}

int cmplxreal(const USignal::Vector<std::complex<double>> &z,
              USignal::Vector<std::complex<double>> &p,
              std::vector<bool> &isreal,
              const double tol)
{
    p.resize(0);
    isreal.resize(0);
    if (z.size() < 1)
    {
        std::cerr << "No elements in z" << std::endl;
        return -1;
    }
    // Set the tolerance
    //double tol = 100.0*DBL_EPSILON;
    //if (tolIn != nullptr){tol = *tolIn;}
    // Get the real poles and complex poles 
    USignal::Vector<double> reals;
    USignal::Vector<std::complex<double>> cmplxs;
    for (int i = 0; i < static_cast<int> (z.size()); ++i)
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
    auto nreal  = reals.size();
    auto ncmplx = cmplxs.size();
    if (nreal + ncmplx != z.size())
    {
        throw std::runtime_error("Algorithmic error - missed");
    }
    // Sort the reals into ascending order and save them 
    nreal = 0;
    if (reals.size() > 0)
    {
        std::sort(reals.begin(), reals.end());
        for (int i = 0; i < static_cast<int> (reals.size()); ++i)
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
                     return std::abs(a) < std::abs(b);
                  });
        // Verify everyone has a buddy (conjugate pair)
        std::vector<bool> lskip(cmplxs.size(), false);
        for (size_t i=0; i<cmplxs.size(); i++)
        {
            if (lskip[i]){continue;}
            // Find it's complex conjugate
            for (size_t j=i+1; j<cmplxs.size(); j++)
            {
                if (lskip[j]){continue;}
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
                    lskip[i] = true;
                    lskip[j] = true;
                    break;
                }
            }
        }
        // I'm not dealing with unpaired
        int npaired = 0;
        for (int i = 0; i < static_cast<int> (lskip.size()); ++i)
        {
            if (lskip[i]){npaired = npaired + 1;}
        }
        //int npaired = std::accumulate(lskip.begin(), lskip.end(), 0);
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
            return -1;
        }
    }
    if (2*nconj + nreal != z.size())
    {
        std::cerr << "Algorithmic failure" << std::endl;
        return -1;
    }
    return 0;
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
                              const USignal::Vector<std::complex<double>> &z,
                              const std::complex<double> p1)
{
    double difMin = std::numeric_limits<double>::max();
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

int getNextWorstRealPoleIndex(const std::vector<bool> &pmask,
                              const std::vector<bool> &isreal_pole,
                              const USignal::Vector<std::complex<double>> &p)
{
    const  std::complex<double> zone(1, 0);
    double difMin = std::numeric_limits<double>::max();
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
                      const USignal::Vector<std::complex<double>> &p)
{
    double difMin = std::numeric_limits<double>::max();
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

int nearest_real_complex_idx(const std::vector<bool> &zmask,
                             const std::vector<bool> &isRealZero,
                             const USignal::Vector<std::complex<double>> &z,
                             const std::complex<double> p1,
                             const bool lreal)
{
    double difMin = std::numeric_limits<double>::max();
    int indx =-1;
    // Look for the closest real zero to the pole p1
    if (lreal)
    {
        for (size_t i=0; i<z.size(); i++)
        {
            if (isRealZero[i] && !zmask[i])
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
            if (!isRealZero[i] && !zmask[i])
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

SecondOrderSections<double> 
zpk2sos(const ZerosPolesGain<double> &zpk,
        const SecondOrderSections<double>::PairingStrategy pairing)
{
    double gain = static_cast<double> (zpk.getGain());
    const auto zeros = zpk.getZerosReference();
    const auto poles = zpk.getPolesReference();
    int nZeros = static_cast<int> (zeros.size());
    int nPoles = static_cast<int> (poles.size());
    if (nPoles <= 0 || nZeros <= 0)
    {
        //std::cerr << "sos will simply scale data by " << gain << std::endl;
        USignal::Vector<double> bs( std::vector<double> {gain, 0, 0} );
        USignal::Vector<double> as( std::vector<double> {1, 0, 0} );
        SecondOrderSections<double> sos(bs, as);
        return sos;
    }
    // Get sizes and handles to zeros and poles
    int np = nPoles + std::max(nZeros - nPoles, 0);
    int nz = nZeros + std::max(nPoles - nZeros, 0);
    int nSections = std::max(np, nz + 1)/2;
#ifndef NDEBUG
    assert(nSections > 0);
#endif
    // Base case
    if (nSections == 1)
    {
        // (z - z_1)*(z - 0) = z^2 - z*z_1 + 0
        USignal::Vector<double> bs({1, 0, 0});
        if (nZeros == 1)
        {
            std::cerr << "nzeros = 1 not yet tested" << std::endl;
            if (std::abs(std::imag(zeros[0])) > 1.e-15)
            {
                std::cerr << "Taking real part of zeros" << std::endl;
            }
            bs[1] =-std::real(zeros[0]);
        }
        else if (nZeros == 2)
        {
            bs[1] =-2.0*std::real(zeros[0]);
            bs[2] = std::pow(std::abs(zeros[0]), 2);
        }
        // Introduce gain
        bs = bs*gain;
        //bs[0] = bs[0]*gain;
        //bs[1] = bs[1]*gain;
        //bs[2] = bs[2]*gain;
        // (p - p_1)*(p - 0) = p^2 - p*p_1 + 0
        USignal::Vector<double> as({1, 0, 0});
        if (nPoles == 1)
        {
            std::cerr << "npoles = 1 not yet tested" << std::endl;
            if (std::abs(std::imag(poles[0])) > 1.e-15)
            {
                std::cerr << "Taking real part of poles" << std::endl;
            }
            as[1] =-std::real(poles[0]);
        }
        // Conjugate pairs: (p - p_1)*(p - conj(p_1)) 
        else if (nPoles == 2)
        {
            as[1] =-2.0*std::real(poles[0]);
            as[2] = std::pow(std::abs(poles[0]), 2);
        }
        SecondOrderSections<double> sos(bs, as); 
        //std::cout <<  "Summary of design" << std::endl;
        //sos.print(stdout);
        return sos;
    }
    if (np%2 == 1 &&
        pairing == SecondOrderSections<double>::PairingStrategy::Nearest)
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
    USignal::Vector<std::complex<double>> z;
    std::vector<bool> isRealZero;
    constexpr double tol{std::numeric_limits<double>::epsilon()*100};
#ifndef NDEBUG
    assert(::cmplxreal(zeros, z, isRealZero, tol) == 0);
#else
    ::cmplxreal(zeros, z, isRealZero, tol);
#endif
    USignal::Vector<std::complex<double>> p;
    std::vector<bool> isRealPole;
#ifdef NDEBUG
    assert(::cmplxreal(poles, p, isRealPole, tol) == 0);
#else
    ::cmplxreal(poles, p, isRealPole, tol); 
#endif
    std::vector<std::array<double, 3>> bsSections;
    std::vector<std::array<double, 3>> asSections;
    std::vector<bool> zmask(z.size(), false);
    std::vector<bool> pmask(p.size(), false);
    for (int is=0; is<nSections; is++)
    {
        std::complex<double> z1(0, 0);
        // Select the `worst' pole
        int p1Index = getWorstPoleIndex(pmask, p);
#ifndef NDEBUG
        assert(p1Index !=-1);
#endif
        std::complex<double> p1 = p[p1Index];
        pmask[p1Index] = true;
        int psum = ::isreal_p_sum(pmask, isRealPole);
        // Initialize complex conjugates for real case
        std::complex<double> p2(0, 0);
        std::complex<double> z2(0, 0);
        // Pair that pole with a zero
        if (isRealPole[p1Index] && psum == 0)
        {
            // Special case to set a first order section
            int z1Index = nearest_real_complex_idx(zmask, isRealZero,
                                                  z, p1, true);
#ifndef NDEBUG
            assert(z1Index !=-1);
#endif
            z1 = z[z1Index];
            zmask[z1Index] = true;
        }
        else
        {
            int zsum = isreal_p_sum(zmask, isRealZero);
            int z1Index =-1;
            if (!isRealPole[p1Index] && zsum == 1)
            {
                // Special case to ensure choose a complex zero to pair
                // with so later (setting up a first-order section)
                z1Index = nearest_real_complex_idx(zmask, isRealZero,
                                                  z, p1, false);
#ifndef NDEBUG
                assert(z1Index !=-1);
                assert(isRealZero[z1Index]);
#endif
            }
            else
            {
                // Pair that pole with the closest zero (real or complex)
                z1Index = nearest_realOrComplex_idx(zmask, z, p1);
#ifndef NDEBUG
                assert(z1Index !=-1);
#endif
            }
            z1 = z[z1Index];
            zmask[z1Index] = true;
            // Now that we have p1 and z1, figure out what p2 and z2 need to be
            if (!isRealPole[p1Index])
            {
                // Complex pole, complex zero -> complex conjugates
                if (!isRealZero[z1Index])
                {
                    p2 = std::conj(p1);
                    z2 = std::conj(z1);
                }
                // Complex pole and real zero
                else
                {
                    p2 = std::conj(p1);
                    int z2Index = nearest_real_complex_idx(zmask, isRealZero,
                                                          z, p1, true);
#ifndef NDEBUG
                    assert(z2Index !=-1);
                    assert(isRealZero[z2Index]);
#endif
                    z2 = z[z2Index];
                    zmask[z2Index] = true;
                }
            }
            // Real pole
            else
            {
                // Real pole, complex zero
                if (!isRealZero[z1Index])
                {
                    z2 = conj(z1);
                    int p2Index = nearest_real_complex_idx(pmask, isRealPole,
                                                          p, z1, true);
#ifndef NDEBUG
                    assert(p2Index !=-1);
                    assert(isRealPole[p2Index]);
#endif
                    p2 = p[p2Index];
                    pmask[p2Index] = true;
                }
                // Real pole, real zero
                else
                {
                    // Pick the next `worst' pole to use
                    int p2Index = getNextWorstRealPoleIndex(pmask,
                                                            isRealPole, p);
#ifndef NDEBUG
                    assert(p2Index !=-1);
                    assert(isRealPole[p2Index]);
#endif
                    p2 = p[p2Index];
                    // Find a real zero to match the added pole
                    int z2Index = nearest_real_complex_idx(zmask, isRealZero,
                                                          z, p2, true);
#ifndef NDEBUG
                    assert(z2Index !=-1);
                    assert(isRealZero[z2Index]);
#endif
                    z2 = z[z2Index];
                    zmask[z2Index] = true; 
                    pmask[p2Index] = true;
                }
            }
        }
        USignal::Vector<std::complex<double>>
            ptemp( std::vector<std::complex<double>> {p1, p2} );
        USignal::Vector<std::complex<double>>
            ztemp( std::vector<std::complex<double>> {z1, z2} );
        double ktemp = 1;
        if (is == nSections - 1){ktemp = gain;}
        USignal::FilterRepresentations::ZerosPolesGain<double>
            zpkTemp(ztemp, ptemp, ktemp);
        USignal::FilterRepresentations::InfiniteImpulseResponse<double>
            baTemp(zpkTemp);
        // Push it back
        auto bsWork = baTemp.getNumeratorFilterCoefficients();
        auto asWork = baTemp.getDenominatorFilterCoefficients();
        std::array<double, 3>
            bSection{bsWork.at(0), bsWork.at(1), bsWork.at(2)};
        std::array<double, 3>
            aSection{asWork.at(0), asWork.at(1), asWork.at(2)};
        bsSections.push_back(std::move(bSection));
        asSections.push_back(std::move(aSection));
/*
        baTemp = zpk2tf(zpkTemp);
        // Save the transfer function of the secon dorder section
        basAll.push_back(baTemp);
*/
    } // Loop on sections 
    // Reality check
    for (int iz = 0; iz < zmask.size(); ++iz)
    {   
        if (!zmask[iz])
        {
#ifndef NDEBUG
            assert(false);
#else
            throw std::runtime_error(
               "Failed to find zero " + std::to_string(iz)
             + " nSections = " + std::to_string(nSections));
#endif
        }
    }   
    for (int ip = 0; ip < pmask.size(); ++ip)
    {
        if (!pmask[ip])
        {
#ifndef NDEBUG
            assert(false);
#else
            throw std::runtime_error(
               "Failed to find pole " + std::to_string(ip)
             + " nSections = " + std::to_string(nSections));
#endif
        }
    }
    // We want the `worst' poles last so reverse
    std::reverse(bsSections.begin(), bsSections.end());
    std::reverse(asSections.begin(), asSections.end());
    // Pack and return the result
    SecondOrderSections<double> result{bsSections, asSections};
    return result;
}

}

template<typename T>
class SecondOrderSections<T>::SecondOrderSectionsImpl
{
public:
    USignal::Vector<T> mNumeratorCoefficients;
    USignal::Vector<T> mDenominatorCoefficients;
    int mSections{0};
};

/// Destructor
template<typename T>
SecondOrderSections<T>::~SecondOrderSections() = default;

/// Constructor
template<typename T>
SecondOrderSections<T>::SecondOrderSections(
    const std::vector<std::array<T, 3>> &numeratorCoefficients,
    const std::vector<std::array<T, 3>> &denominatorCoefficients)
{
    auto b = ::reshape(numeratorCoefficients, true);
    auto a = ::reshape(denominatorCoefficients, false);
    SecondOrderSections<T> sos{b, a};
    *this = std::move(sos);
}

template<>
SecondOrderSections<double>::SecondOrderSections(
    const ZerosPolesGain<double> &zpk,
    const SecondOrderSections::PairingStrategy strategy)
{
    *this = ::zpk2sos(zpk, strategy); 
}

template<>
SecondOrderSections<float>::SecondOrderSections(
    const ZerosPolesGain<float> &zpk,
    const SecondOrderSections::PairingStrategy strategy)
{
    double gain = static_cast<double> (zpk.getGain());
    auto zerosFloat = zpk.getZeros();
    auto polesFloat = zpk.getPoles();
    USignal::Vector<std::complex<double>> zeros(zerosFloat.size());
    for (int i = 0; i < static_cast<int> (zerosFloat.size()); ++i)
    {
        zeros[i] = zerosFloat[i];
    }
    USignal::Vector<std::complex<double>> poles(polesFloat.size());
    for (int i = 0; i < static_cast<int> (polesFloat.size()); ++i)
    {
        poles[i] = polesFloat[i];
    }
    ZerosPolesGain<double> zpkDouble(zeros, poles, gain);
    auto strategyDouble
        = static_cast<SecondOrderSections<double>::PairingStrategy> (strategy);
    auto sos = ::zpk2sos(zpkDouble, strategyDouble); 
    auto bs = sos.getNumeratorFilterCoefficients();
    auto as = sos.getDenominatorFilterCoefficients();
    USignal::Vector<float> bsFloat(bs.size());
    for (int i = 0; i < static_cast<int> (bs.size()); ++i)
    {
        bsFloat[i] = static_cast<float> (bs[i]);
    }
    USignal::Vector<float> asFloat(as.size());
    for (int i = 0; i < static_cast<int> (as.size()); ++i)
    {   
        asFloat[i] = static_cast<float> (as[i]);
    }   
    *this = SecondOrderSections<float> (bsFloat, asFloat);
}


/// Constructor
template<typename T>
SecondOrderSections<T>::SecondOrderSections(
    const USignal::Vector<T> &numeratorCoefficients,
    const USignal::Vector<T> &denominatorCoefficients) :
    pImpl(std::make_unique<SecondOrderSectionsImpl> ())
{
    if (numeratorCoefficients.size() !=
        denominatorCoefficients.size())
    {
        throw std::invalid_argument(
           "numeratorCoefficients.size() != denominatorCoefficients.size()");
    }
    if (numeratorCoefficients.size() < 3)
    {
        throw std::invalid_argument(
           "At least 3 numerator and denominator coefficients required");
    }
    if (numeratorCoefficients.size()%3 != 0)
    {
        throw std::invalid_argument( 
           "Number of coefficients must be multiple of 3");
    }
    int nSections = static_cast<int> (numeratorCoefficients.size())/3;
    for (int section = 0; section < nSections; ++section)
    {
        auto numeratorCoefficient = numeratorCoefficients[3*section];
        if (numeratorCoefficient == 0)
        {
            throw std::invalid_argument(
                "Leading numerator coefficient for "
                + std::to_string(section + 1)
                + " section is zero");
        }
        auto denominatorCoefficient = denominatorCoefficients[3*section];
        if (denominatorCoefficient == 0)
        {
            throw std::invalid_argument(
                "Leading denominator coefficient for "
                + std::to_string(section + 1)
                + " section is zero");
        }
    }
    pImpl->mSections = nSections;
    pImpl->mNumeratorCoefficients = numeratorCoefficients;
    pImpl->mDenominatorCoefficients = denominatorCoefficients; 
}

/// Copy constructor
template<typename T>
SecondOrderSections<T>::SecondOrderSections(const SecondOrderSections<T> &sos)
{
    *this = sos;
}

/// Move constructor
template<typename T>
SecondOrderSections<T>::SecondOrderSections(
    SecondOrderSections<T> &&sos) noexcept
{
    *this = std::move(sos);
}

/// Copy assignment
template<typename T>
SecondOrderSections<T>& SecondOrderSections<T>::operator=(
    const SecondOrderSections &sections)
{
    if (&sections == this){return *this;}
    pImpl = std::make_unique<SecondOrderSectionsImpl> (*sections.pImpl);
    return *this;
}

/// Move assignment
template<typename T>
SecondOrderSections<T>& SecondOrderSections<T>::operator=(
    SecondOrderSections &&sections) noexcept
{
    if (&sections == this){return *this;}
    pImpl = std::move(sections.pImpl);
    return *this;
}

/// Numerator coefficients
template<class T>
const USignal::Vector<T> &
SecondOrderSections<T>::getNumeratorFilterCoefficientsReference() const noexcept
{
    return *&pImpl->mNumeratorCoefficients;
}

template<class T>
USignal::Vector<T> 
SecondOrderSections<T>::getNumeratorFilterCoefficients() const noexcept
{
    return pImpl->mNumeratorCoefficients;
}

/// Numerator coefficients
template<class T>
const USignal::Vector<T> &
SecondOrderSections<T>::getDenominatorFilterCoefficientsReference() const noexcept
{
    return *&pImpl->mDenominatorCoefficients;
}

template<class T>
USignal::Vector<T> 
SecondOrderSections<T>::getDenominatorFilterCoefficients() const noexcept
{
    return pImpl->mDenominatorCoefficients;
}

/// Number of sections
template<class T>
int SecondOrderSections<T>::getNumberOfSections() const noexcept
{
    return pImpl->mSections;
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class USignal::FilterRepresentations::SecondOrderSections<double>;
template class USignal::FilterRepresentations::SecondOrderSections<float>;
template class USignal::FilterRepresentations::SecondOrderSections<int>;

