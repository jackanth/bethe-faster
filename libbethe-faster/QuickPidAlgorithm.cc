/**
 *  @file   bethe-faster/src/QuickPidAlgorithm.cc
 *
 *  @brief  Implementation of the QuickPidAlgorithm class.
 *
 *  $Log: $
 */

#include "QuickPidAlgorithm.h"
#include "PhysicalConstants.h"

#include <iostream>

namespace bf
{

QuickPidAlgorithm::QuickPidAlgorithm(Detector detector) :
    m_detector{std::move_if_noexcept(detector)},
    m_xiPrimePartial{0.},
    m_chiPartial{0.}
{
    if (2. * m_detector.m_atomicMass < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize quick PID algorithm because the atomic mass was too small");

    m_xiPrimePartial = PhysicalConstants::m_kCoefficient * static_cast<double>(m_detector.m_atomicNumber) * m_detector.m_density /
                  (2. * m_detector.m_atomicMass);

    const double denominator = m_detector.m_avgIonizationEnergy * m_detector.m_avgIonizationEnergy * 1.e-12;

    if (denominator < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize quick PID algorithm because the average ionization energy was too small");

    m_chiPartial = std::log(2. * PhysicalConstants::m_electronMass / denominator) + PhysicalConstants::m_vavilovJ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

double QuickPidAlgorithm::CalculateScaledKineticEnergy(const double residualRange, const double mass, const double deltaX) const
{
    const double xiPrime = this->XiPrime(deltaX);
    const double chi = this->Chi(xiPrime);
    const double prefactor = mass * deltaX / (2. * xiPrime);

    // Newton's method
    double scaledKineticEnergy = 0.2; // initial estimate
    double func = this->NewtonMethodFunc(scaledKineticEnergy, residualRange, prefactor, chi);

    while (std::abs(func) > 0.000001)
    {
        double funcDeriv = this->NewtonMethodFuncDeriv(scaledKineticEnergy, prefactor, chi);
        scaledKineticEnergy -= func / funcDeriv;
        func = this->NewtonMethodFunc(scaledKineticEnergy, residualRange, prefactor, chi);
    }

    if (scaledKineticEnergy < -std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"Newton's method failed to find the positive root"};

    return scaledKineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double QuickPidAlgorithm::NewtonMethodFunc(const double scaledKineticEnergy, const double residualRange, const double prefactor, const double chi) const
{
    const double logTerm = std::log(1. + 2. * scaledKineticEnergy * scaledKineticEnergy / chi) / 3.;
    
    const double arctanPrefactor = std::sqrt(0.5 * chi);
    const double arctanTerm = arctanPrefactor * std::atan(scaledKineticEnergy / arctanPrefactor);

    return 3. * prefactor * (logTerm + arctanTerm - scaledKineticEnergy) - residualRange;
}

} // namespace bf
