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

#include "gsl/gsl_rng.h"

#include <memory>
#include <set>
#include <iostream>

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
     *  @param  deltaX the distance to propagate
     *
     *  @return the dQ/dx response of the detector
     */
    float Propagate(const std::shared_ptr<Particle> &spParticle, const float deltaX) const;

    /**
     *  @brief  Propagate a set of particle until they have all stopped
     *
     *  @param  deltaX the distance to propagate each iteration
     *  @param  particleSet the set of particles
     */
    void PropagateUntilStopped(const float deltaX, const ParticleHelper::ParticleSet &particleSet) const;

    /**
     *  @brief  Propagate a particle until it stops
     *
     *  @param  deltaX the distance to propagate each iteration
     *  @param  spParticle shared pointer to the particle
     */
    void PropagateUntilStopped(const float deltaX, const std::shared_ptr<Particle> &spParticle) const;

    /**
     *  @brief  Calculate the detector dQ/dx value
     *
     *  @param  dEdx the dE/dx value
     *
     *  @return the dQ/dx value
     */
    float CalculateResponse(const float dEdx) const;

    /**
     *  @brief  Invert the detector response to get dE/dx
     *
     *  @param  dQdx the dQ/dx value
     *
     *  @return the dE/dx value
     */
    float InvertResponse(const float dQdx) const;

    /**
     *  @brief  Calculate Vavilov's kappa
     *
     *  @param  mass the particle mass
     *  @param  energy the particle energy
     *
     *  @return the kappa value
     */
    float CalculateKappa(const float mass, const float energy) const;

protected:
    /**
     *  @brief  Calculate the probability of an observation for a given particle
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  observeddEdx the observed dE/dx value
     *  @param  deltaX the position increment
     * 
     *  @return the observation probability
     */
    float CalculateObservationProbability(const std::shared_ptr<Particle> &spParticle, const float observeddEdx, const float deltaX) const;

    /**
     *  @brief  Calculate the probability of an energy transition
     *
     *  @param  deltaE the change in energy
     *  @param  mass the mass
     *  @param  previousEnergy the previous energy
     * 
     *  @return the transition probability
     */
    float CalculateTransitionProbability(const float deltaE, const float mass, const float previousEnergy) const;

    friend class ParticleFilter;

private:
    Detector  m_detector;    ///< The detector parameters
    float     m_cBar;        ///< The value of cBar
    float     m_chargeMapC1; ///< The value of the first charge map constant
    float     m_chargeMapC2; ///< The value of the second charge map constant
    float     m_xiPartial;   ///< The value of partial xi
    gsl_rng * m_pGenerator;  ///< Address of the GSL random number generator

    /**
     *  @brief  Get the xi value
     *
     *  @param  beta the value of beta
     *
     *  @return the value of xi
     */
    float Xi(const float beta) const;

    /**
     *  @brief  Get the density correction value
     *
     *  @param  beta the value of beta
     *
     *  @return the density correction
     */
    float DensityCorrection(const float beta) const;

    /**
     *  @brief  Sample from the delta E distribution
     *
     *  @param  meanEnergyLoss the mean energy loss
     *  @param  kappa the kappa value
     *  @param  beta2 the beta^2 value
     *  @param  xi the xi value
     *  @param  maxEnergyTransfer the max energy transfer
     *
     *  @return the sample
     */
    float SampleDeltaE(const float meanEnergyLoss, const float kappa, const float beta2, const float xi, const float maxEnergyTransfer) const;

    /**
     *  @brief  Sample from the modified Rutherford cross-section distribution
     *
     *  @param  eMin the min energy loss
     *  @param  eMax the max energy loss
     *  @param  beta2 the beta^2 value
     *
     *  @return the sample
     */
    float SampleFromRutherfordDistribution(const float eMin, const float eMax, const float beta2) const;

    /**
     *  @brief  Get the maximum energy transfer
     *
     *  @param  beta the particle beta
     *  @param  mass the particle mass
     *
     *  @return the maximum energy transfer
     */
    float GetMaxEnergyTransfer(const float beta, const float mass) const;

    /**
     *  @brief  Get the expected energy loss
     *
     *  @param  beta the beta value
     *  @param  beta2 the beta^2 value
     *  @param  xi the xi value
     *
     *  @return the expected energy loss
     */
    float GetExpectedEnergyLoss(const float beta, const float beta2, const float xi) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline void Propagator::PropagateUntilStopped(const float deltaX, const ParticleHelper::ParticleSet &particleSet) const
{
    for (const std::shared_ptr<Particle> &spParticle : particleSet)
        this->PropagateUntilStopped(deltaX, spParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Propagator::PropagateUntilStopped(const float deltaX, const std::shared_ptr<Particle> &spParticle) const
{
    while (spParticle->KineticEnergy() > 0.f)
        this->Propagate(spParticle, deltaX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Propagator::Xi(const float beta) const
{
    const float beta2 = beta * beta;

    if (beta2 < std::numeric_limits<float>::epsilon())
        throw std::runtime_error("Could not calculate xi because beta^2 was too small");

    return m_xiPartial / beta2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Propagator::CalculateResponse(const float dEdx) const
{
    return -m_chargeMapC1 * dEdx / (1.f - m_chargeMapC2 * dEdx);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Propagator::InvertResponse(const float dQdx) const
{
    return std::min(-dQdx / (m_chargeMapC1 - dQdx * m_chargeMapC2), 0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Propagator::CalculateKappa(const float mass, const float energy) const
{
    const float beta = ParticleHelper::GetParticleBeta(mass, energy);
    return this->Xi(beta) / this->GetMaxEnergyTransfer(beta, mass); 
}

} // namespace bf

#endif // #ifndef BF_PROPAGATOR_H
