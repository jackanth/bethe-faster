/**
 *  @file   bethe-faster/include/Propagator.h
 *
 *  @brief  Header file for the propagator class.
 *
 *  $Log: $
 */

#ifndef BF_PROPAGATOR_H
#define BF_PROPAGATOR_H 1

#include "Detector.h"
#include "Particle.h"
#include "ParticleHelper.h"

#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"
#include "Math/VavilovFast.h"

#include <iostream>
#include <map>
#include <memory>
#include <set>

namespace bf
{

/**
 *  @brief  Forward declaration of ParticleFilter class
 */
class ParticleFilter;

/**
 *  @brief  Propagator class
 */
class Propagator
{
public:
    /**
     *  @brief  An enumeration of the different types of propagation
     */
    enum class PROPAGATION_MODE
    {
        STOCHASTIC, ///< Stochastic propagation
        MEAN,       ///< Mean propagation
        MODAL       ///< Modal propagation
    };

    /**
     *  @brief  Constructor
     *
     *  @param  detector the detector
     */
    explicit Propagator(Detector detector);

    /**
     *  @brief  Propagate a single step particle backwards
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  stepSize the step size
     *  @param  mode the propagation mode
     *
     *  @return the dE/dx value
     */
    double PropagateBackwards(const std::shared_ptr<Particle> &spParticle, const double stepSize,
        const PROPAGATION_MODE mode = PROPAGATION_MODE::STOCHASTIC) const;

    /**
     *  @brief  Calculate Vavilov's kappa
     *
     *  @param  mass the particle mass
     *  @param  energy the particle energy
     *  @param  deltaX the effective thickness
     *
     *  @return the kappa value
     */
    double CalculateKappa(const double mass, const double energy, const double deltaX) const;

    /**
     *  @brief  Get the density correction value
     *
     *  @param  mass the particle mass
     *  @param  energy the particle energy
     *
     *  @return the density correction
     */
    double DensityCorrection(const double mass, const double energy) const;

    /**
     *  @brief  Get the xi value
     *
     *  @param  mass the particle mass
     *  @param  energy the particle energy
     *  @param  deltaX the effective thickness
     *
     *  @return the value of xi
     */
    double Xi(const double mass, const double energy, const double deltaX) const;

protected:
    /**
     *  @brief  Calculate the probability of an observation for a given particle
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  observeddEdx the observed dE/dx value
     *  @param  observeddx the observed effective thickness
     *
     *  @return the observation probability
     */
    double CalculateObservationProbability(const std::shared_ptr<Particle> &spParticle, const double observeddEdx, const double observeddx) const;

    /**
     *  @brief  Calculate the probability of an energy transition
     *
     *  @param  deltaE the change in energy
     *  @param  mass the mass
     *  @param  previousEnergy the previous energy
     *  @param  deltaX the effective thickness
     *
     *  @return the transition probability
     */
    double CalculateTransitionProbability(const double deltaE, const double mass, const double previousEnergy, const double deltaX) const;

    /**
     *  @brief  Sample from the Vavilov distribution
     *
     *  @param  kappa the kappa value
     *  @param  beta2 the beta^2 value
     *  @param  maxError the maximum error
     *  @param  maxIterations the maximum number of iterations
     *
     *  @return the sample
     */
    double SampleFromVavilovDistribution(
        const double kappa, const double beta2, const double maxError = 0.01, const std::size_t maxIterations = 10000UL) const;

    /**
     *  @brief  Get the Vavilov probability density
     *
     *  @param  kappa the kappa value
     *  @param  beta2 the beta^2 value
     *  @param  lambda the lambda value
     *
     *  @return the probability density
     */
    double GetVavilovProbabilityDensity(const double kappa, const double beta2, const double lambda) const;

    /**
     *  @brief  Get the Vavilov mode
     *
     *  @param  kappa the kappa value
     *  @param  beta2 the beta^2 value
     *
     *  @return the mode
     */
    double GetVavilovMode(const double kappa, const double beta2) const;

    friend class ParticleFilter;

private:
    Detector                                         m_detector;  ///< The detector parameters
    double                                           m_cBar;      ///< The value of cBar
    double                                           m_xiPartial; ///< The value of partial xi
    mutable ROOT::Math::Random<ROOT::Math::GSLRngMT> m_randomGen; ///< A random number generator

    /**
     *  @brief  Get the xi value
     *
     *  @param  beta the beta^2 value
     *  @param  deltaX the effective thickness
     *
     *  @return the value of xi
     */
    double Xi(const double beta2, const double deltaX) const;

    /**
     *  @brief  Get the density correction value
     *
     *  @param  beta the value of beta
     *
     *  @return the density correction
     */
    double DensityCorrection(const double beta) const;

    /**
     *  @brief  Sample from the delta E distribution
     *
     *  @param  meanEnergyLoss the mean energy loss
     *  @param  kappa the kappa value
     *  @param  beta2 the beta^2 value
     *  @param  xi the xi value
     *
     *  @return the sample
     */
    double SampleDeltaE(const double meanEnergyLoss, const double kappa, const double beta2, const double xi) const;

    /**
     *  @brief  Get the mode of the delta E distribution
     *
     *  @param  meanEnergyLoss the mean energy loss
     *  @param  kappa the kappa value
     *  @param  beta2 the beta^2 value
     *  @param  xi the xi value
     *
     *  @return the mode
     */
    double GetDeltaEMode(const double meanEnergyLoss, const double kappa, const double beta2, const double xi) const;

    /**
     *  @brief  Get the maximum energy transfer
     *
     *  @param  beta the particle beta
     *  @param  mass the particle mass
     *
     *  @return the maximum energy transfer
     */
    double GetMaxEnergyTransfer(const double beta, const double mass) const;

    /**
     *  @brief  Get the expected energy loss
     *
     *  @param  beta the beta value
     *  @param  beta2 the beta^2 value
     *  @param  xi the xi value
     *  @param  maxEnergyLoss the maximum energy loss
     *
     *  @return the expected energy loss
     */
    double GetExpectedEnergyLoss(const double beta, const double beta2, const double xi, const double maxEnergyLoss) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double Propagator::Xi(const double beta2, const double deltaX) const
{
    if (beta2 < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not calculate xi because beta^2 was too small");

    return m_xiPartial * deltaX / beta2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Propagator::CalculateKappa(const double mass, const double energy, const double deltaX) const
{
    const double beta = ParticleHelper::GetParticleBeta(mass, energy);
    return this->Xi(beta * beta, deltaX) / this->GetMaxEnergyTransfer(beta, mass);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Propagator::DensityCorrection(const double mass, const double energy) const
{
    const double beta = ParticleHelper::GetParticleBeta(mass, energy);
    return this->DensityCorrection(beta);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Propagator::Xi(const double mass, const double energy, const double deltaX) const
{
    const double beta = ParticleHelper::GetParticleBeta(mass, energy);
    return this->Xi(beta * beta, deltaX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Propagator::CalculateObservationProbability(const std::shared_ptr<Particle> &spParticle, const double observeddEdx, const double observeddx) const
{
    return this->CalculateTransitionProbability(observeddEdx * observeddx, spParticle->Mass(), spParticle->KineticEnergy(), observeddx);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Propagator::GetVavilovProbabilityDensity(const double kappa, const double beta2, const double lambda) const
{
    return ROOT::Math::VavilovFast{kappa, beta2}.Pdf(lambda);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Propagator::GetVavilovMode(const double kappa, const double beta2) const
{
    return ROOT::Math::VavilovFast{kappa, beta2}.Mode();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Propagator::GetExpectedEnergyLoss(const double beta, const double beta2, const double xi, const double maxEnergyLoss) const
{
    const double lnArg = 2. * PhysicalConstants::m_electronMass * beta2 * maxEnergyLoss /
                         ((1. - beta2) * m_detector.m_avgIonizationEnergy * m_detector.m_avgIonizationEnergy * 1.e-12);

    return xi * (std::log(lnArg) - 2. * beta2 - this->DensityCorrection(beta));
}

} // namespace bf

#endif // #ifndef BF_PROPAGATOR_H
