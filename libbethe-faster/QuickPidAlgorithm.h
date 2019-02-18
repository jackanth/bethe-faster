/**
 *  @file   bethe-faster/include/QuickPidAlgorithm.h
 *
 *  @brief  Header file for the quick PID algorithm class.
 *
 *  $Log: $
 */

#ifndef BF_QUICK_PID_ALGORITHM_H
#define BF_QUICK_PID_ALGORITHM_H 1

#include "Detector.h"
#include "HitCharge.h"
#include "QuickPidHelper.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace bf
{

/**
 *  @brief  QuickPidAlgorithm class
 */
class QuickPidAlgorithm
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  detector the detector
     */
    explicit QuickPidAlgorithm(Detector detector);

    /**
     *  @brief  Estimate the dE/dx using Newton's method
     *
     *  @param  residualRange the residual range
     *  @param  mass the mass
     *  @param  deltaX the effective detector thickness
     *
     *  @return the scaled kinetic energy
     */
    double EstimatedEdx(const double residualRange, const double mass, const double deltaX) const;

    /**
     *  @brief  Calculate the scaled kinetic energy using Newton's method
     *
     *  @param  residualRange the residual range
     *  @param  mass the mass
     *  @param  deltaX the effective detector thickness
     *
     *  @return the scaled kinetic energy
     */
    double CalculateScaledKineticEnergy(const double residualRange, const double mass, const double deltaX) const;

    /**
     *  @brief  Calculate the Bragg gradient and intercept
     *
     *  @param  hitChargeVector the hit charge vector
     *
     *  @return the Bragg gradient and intercept
     */
    static std::tuple<double, double> CalculateBraggGradient(std::vector<HitCharge> hitChargeVector);

private:
    Detector m_detector;       ///< The detector parameters
    double   m_xiPrimePartial; ///< The value of partial xi prime
    double   m_chiPartial;     ///< The value of partial chi

    /**
     *  @brief  Get the xi prime value
     *
     *  @param  deltaX the effective thickness
     *
     *  @return the value of xi prime
     */
    double XiPrime(const double deltaX) const;

    /**
     *  @brief  Get the chi value
     *
     *  @param  xiPrime the value of xi prime
     *
     *  @return the value of chi
     */
    double Chi(const double xiPrime) const;

    /**
     *  @brief  Calculate the scaled kinetic energy function for Newton's method
     *
     *  @parma  scaledKineticEnergy the scaled kinetic energy
     *  @param  residualRange the residual range
     *  @param  prefactor the prefactor
     *  @param  chi the chi value
     *
     *  @return the value of the function
     */
    double NewtonMethodFunc(const double scaledKineticEnergy, const double residualRange, const double prefactor, const double chi) const;

    /**
     *  @brief  Calculate the derivative of the scaled kinetic energy function for Newton's method
     *
     *  @parma  scaledKineticEnergy the scaled kinetic energy
     *  @param  prefactor the prefactor
     *  @param  chi the chi value
     *
     *  @return the value of the derivative
     */
    double NewtonMethodFuncDeriv(const double scaledKineticEnergy, const double prefactor, const double chi) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline std::tuple<double, double> QuickPidAlgorithm::CalculateBraggGradient(std::vector<HitCharge> hitChargeVector)
{
    return QuickPidHelper::CalculateBraggGradient(hitChargeVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double QuickPidAlgorithm::EstimatedEdx(const double residualRange, const double mass, const double deltaX) const
{
    const double scaledKineticEnergy = this->CalculateScaledKineticEnergy(residualRange, mass, deltaX);

    const double xiPrime = this->XiPrime(deltaX);
    const double chi     = this->Chi(xiPrime);

    const double gamma = 1. + scaledKineticEnergy;
    const double beta2 = 1. - 1. / (gamma * gamma);

    return xiPrime / (beta2 * deltaX) * (chi + 2. * std::log(gamma) - beta2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double QuickPidAlgorithm::XiPrime(const double deltaX) const
{
    return m_xiPrimePartial * deltaX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double QuickPidAlgorithm::Chi(const double xiPrime) const
{
    return m_chiPartial + std::log(xiPrime);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double QuickPidAlgorithm::NewtonMethodFuncDeriv(const double scaledKineticEnergy, const double prefactor, const double chi) const
{
    const double fracTerm = (3. * chi + 4. * scaledKineticEnergy) / (chi + 2. * scaledKineticEnergy * scaledKineticEnergy);
    return prefactor * (fracTerm - 3.);
}

} // namespace bf

#endif // #ifndef BF_QUICK_PID_ALGORITHM_H
