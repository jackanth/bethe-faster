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

#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"

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
     *  @brief  Constructor
     *
     *  @param  detector the detector
     */
    explicit Propagator(Detector detector);

    /**
     *  @brief  Destructor
     */
    ~Propagator();

    /**
     *  @brief  Propagate a particle
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  deltaX the effective thickness
     *
     *  @return the dE/dx value
     */
    double Propagate(const std::shared_ptr<Particle> &spParticle, const double deltaX) const;

    /**
     *  @brief  Propagate a set of particles until they have all stopped
     *
     *  @param  deltaX the effective thickness
     *  @param  particleSet the set of particles
     */
    void PropagateUntilStopped(const double deltaX, const ParticleHelper::ParticleSet &particleSet) const;

    /**
     *  @brief  Propagate a particle until it stops
     *
     *  @param  deltaX the effective thickness
     *  @param  spParticle shared pointer to the particle
     */
    void PropagateUntilStopped(const double deltaX, const std::shared_ptr<Particle> &spParticle) const;

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
    double SampleFromVavilovDistribution(const double kappa, const double beta2, const double maxError = 0.01, const std::size_t maxIterations = 10000UL) const;

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

    friend class ParticleFilter;

private:
    Detector                       m_detector;               ///< The detector parameters
    double                         m_cBar;                   ///< The value of cBar
    double                         m_xiPartial;              ///< The value of partial xi
    ROOT::Math::Random<ROOT::Math::GSLRngMT> *                      m_pRandom;                ///< Address of a ROOT TRandom object

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

inline void Propagator::PropagateUntilStopped(const double deltaX, const ParticleHelper::ParticleSet &particleSet) const
{
    for (const std::shared_ptr<Particle> &spParticle : particleSet)
        this->PropagateUntilStopped(deltaX, spParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Propagator::PropagateUntilStopped(const double deltaX, const std::shared_ptr<Particle> &spParticle) const
{
    while (spParticle->IsAlive())
        this->Propagate(spParticle, deltaX);
}

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

} // namespace bf

#endif // #ifndef BF_PROPAGATOR_H
