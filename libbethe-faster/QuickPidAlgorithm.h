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

#include "TGraph.h"

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
    static bool CalculateBraggGradient(std::vector<HitCharge> hitChargeVector, double &gradient, double &intercept);

    bool CalculateSecondOrderBraggGradient(std::vector<HitCharge> hitChargeVector, double &gradient, double &intercept) const;

    TGraph GetSecondOrderBraggGradientGraph(const std::vector<HitCharge> &hitChargeVector) const;

private:
    Detector m_detector;   ///< The detector parameters
    double   m_xiPrime;    ///< The value of xi prime
    double   m_chiPartial; ///< The value of partial chi

    /**
     *  @brief  Get the chi value
     *
     *  @param  deltaX the value of deltaX
     *
     *  @return the value of chi
     */
    double Chi(const double deltaX) const;

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

    double CalculateTPrime(const double dEdx, const double deltaX) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline bool QuickPidAlgorithm::CalculateBraggGradient(std::vector<HitCharge> hitChargeVector, double &gradient, double &intercept)
{
    return QuickPidHelper::CalculateBraggGradient(hitChargeVector, gradient, intercept);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double QuickPidAlgorithm::EstimatedEdx(const double residualRange, const double mass, const double deltaX) const
{
    const double scaledKineticEnergy = this->CalculateScaledKineticEnergy(residualRange, mass, deltaX);
    const double chi     = this->Chi(deltaX);

    const double gamma = 1. + scaledKineticEnergy;
    const double beta2 = 1. - 1. / (gamma * gamma);

    return m_xiPrime / beta2 * (chi + 2. * std::log(gamma) - beta2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double QuickPidAlgorithm::Chi(const double deltaX) const
{
    return m_chiPartial + std::log(deltaX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double QuickPidAlgorithm::NewtonMethodFuncDeriv(const double scaledKineticEnergy, const double prefactor, const double chi) const
{
    const double fracTerm = (3. * chi + 4. * scaledKineticEnergy) / (chi + 2. * scaledKineticEnergy * scaledKineticEnergy);
    return prefactor * (fracTerm - 3.);
}

} // namespace bf

#endif // #ifndef BF_QUICK_PID_ALGORITHM_H
