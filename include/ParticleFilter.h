/**
 *  @file   bethe-faster/include/ParticleFilter.h
 *
 *  @brief  Header file for the particle filter class.
 *
 *  $Log: $
 */

#ifndef BF_PARTICLE_FILTER_H
#define BF_PARTICLE_FILTER_H 1

#include "Particle.h"
#include "Propagator.h"

#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"

#include <memory>
#include <vector>

namespace bf
{

/**
 *  @brief  ObservedParticleState class
 */
class ObservedParticleState
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  residualRange the residual range (cm)
     *  @param  dEdx the dE/dx value (MeV/cm)
     *  @param  dx the change in position (cm)
     */
    ObservedParticleState(const double residualRange, const double dEdx, const double dx) noexcept;

    /**
     *  @brief  Get the residual range (cm)
     *
     *  @return the residual range
     */
    double GetResidualRange() const noexcept;

    /**
     *  @brief  Get the dE/dx value (MeV/cm)
     *
     *  @return the dE/dx value
     */
    double GetdEdx() const noexcept;

    /**
     *  @brief  Get the dx value (cm)
     *
     *  @return the dx value
     */
    double Getdx() const noexcept;

private:
    double m_residualRange; ///< The residual range (cm)
    double m_dEdx;          ///< The dE/dx value (MeV/cm)
    double m_dx;            ///< The dx value (cm)
};

/**
 *  @brief  MassHypothesis class
 */
class MassHypothesis
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  mass the mass
     *  @param  priorFinalEnergyDistribution the prior final energy distribution
     */
    MassHypothesis(const double mass, std::vector<double> priorFinalEnergyDistribution);

    /**
     *  @brief  Constructor
     *
     *  @param  mass the mass
     *  @param  finalEnergy the final energy
     *  @param  numberOfParticles the number of particles
     */
    MassHypothesis(const double mass, const double finalEnergy, const std::size_t numberOfParticles = 1000UL);

    /**
     *  @brief  Get the mass
     *
     *  @return the mass
     */
    double Mass() const noexcept;

    /**
     *  @brief  Get the number of particles
     *
     *  @return the number of particles
     */
    std::size_t NumberOfParticles() const noexcept;

    /**
     *  @brief  Get the prior final energy distribution
     *
     *  @return the prior final energy distribution
     */
    const std::vector<double> &PriorFinalEnergyDistribution() const noexcept;

private:
    std::shared_ptr<Propagator> m_spPropagator;                 ///< The propagator
    double                      m_mass;                         ///< The particle mass
    std::vector<double>         m_priorFinalEnergyDistribution; ///< The prior final energy distribution
    std::size_t                 m_numberOfParticles;            ///< The number of particles
};

/**
 *  @brief  ParticleFilter class
 */
class ParticleFilter
{
public:
    using WeightedParticle     = std::pair<double, std::shared_ptr<Particle>>;     ///< Alias for a particle with a weight
    using ParticleDistribution = std::vector<WeightedParticle>;                    ///< Alias for a particle distribution
    using ObservedStateVector  = std::vector<ObservedParticleState>;               ///< Alias for a vector of observed particle states
    using DistributionRecord   = std::vector<std::pair<double, double>>;           ///< Alias for a distribution record
    using DistributionHistory  = std::vector<std::shared_ptr<DistributionRecord>>; ///< Alias for a distribution history

    /**
     *  @brief  Constructor
     *
     *  @param  detector the detector
     *  @param  options the options
     */
    ParticleFilter(Detector detector);

    /**
     *  @brief  Filter a uniform energy distribution on a particle history
     *
     *  @param  history the particle histroy
     *  @param  massHypothesis the mass hypothesis
     *
     *  @return the distribution history
     */
    DistributionHistory Filter(const Particle::History &history, const MassHypothesis &massHypothesis) const;

    /**
     *  @brief  Filter a uniform energy distribution on a set of observations
     *
     *  @param  observations the ordered vector of observations
     *  @param  massHypothesis the mass hypothesis
     *
     *  @return the distribution history
     */
    DistributionHistory Filter(ObservedStateVector observations, const MassHypothesis &massHypothesis) const;

    /**
     *  @brief  Calculate the marginal likelihood, given a distribution history
     *
     *  @param  distributionHistory the distribution history
     *
     *  @return the marginal likelihood
     */
    static double CalculateMarginalLikelihood(const DistributionHistory &distributionHistory);

    /**
     *  @brief  Calculate the final energy mean and standard deviation, given a distribution history
     *
     *  @param  distributionHistory the distribution history
     *
     *  @return the mean final energy and the standard deviation of the final energy
     */
    static std::tuple<double, double> CalculateFinalEnergy(const DistributionHistory &distributionHistory);

    /**
     *  @brief  Calculate PIDA for a particle history
     *
     *  @param  history the particle history
     *
     *  @return the PIDA value
     */
    static double CalculatePida(const Particle::History &history);

    /**
     *  @brief  Calculate PIDA for a set of observations
     *
     *  @param  observations the observations
     *
     *  @return the PIDA value
     */
    static double CalculatePida(const ObservedStateVector &observations);

private:
    std::shared_ptr<Propagator>                      m_spPropagator;    ///< The propagator
    std::shared_ptr<ParticleDistribution>            m_spDistribution;  ///< The current distribution
    mutable ROOT::Math::Random<ROOT::Math::GSLRngMT> m_randomGenerator; ///< A ROOT TRandom object

    /**
     *  @brief  Filter on an observation
     *
     *  @param  observation the observation
     *  @param  isFirstObservation whether this is the first observation
     *  @param  nSteps the number of propagation steps
     *  @param  spResamplingProbabilityVector shared pointer to the resampling probability vector
     *  @param  spResamplingParticleVector shared pointer to the resampling particle vector
     *  @param  numParticles the number of particles
     *  @param  stepSize the step size
     */
    bool FilterOnObservation(const ObservedParticleState &observation, const bool isFirstObservation, const std::size_t nSteps,
        const std::shared_ptr<std::vector<double>> &spResamplingProbabilityVector,
        const std::shared_ptr<std::vector<unsigned int>> &spResamplingParticleVector, const std::size_t numParticles, const double stepSize) const;

    /**
     *  @brief  Check the distribution weights
     *
     *  @return Whether the weights valid
     */
    bool CheckDistributionWeights() const;

    /**
     *  @brief  Check the mass hypothesis
     *
     *  @param  massHypothesis the mass hypothesis
     */
    void CheckMassHypothesis(const MassHypothesis &massHypothesis) const;

    /**
     *  @brief  Create a uniform prior over the initial energy bounds
     *
     *  @param  massHypothesis the mass hypothesis
     *  @param  initialResidualRange the initial residual range
     */
    void CreateDistribution(const MassHypothesis &massHypothesis, const double initialResidualRange) const;

    /**
     *  @brief  Propagate particles in the distribution
     *
     *  @param  nSteps the number of steps
     *  @param  stepSize the step size
     */
    void PropagateParticles(const std::size_t nSteps, const double stepSize) const;

    /**
     *  @brief  Resample the distribution
     *
     *  @param  spResamplingProbabilityVector shared pointer to the resampling probability vector
     *  @param  spResamplingParticleVector shared pointer to the resampling particle vector
     */
    void ResampleDistribution(const std::shared_ptr<std::vector<double>> &spResamplingProbabilityVector,
        const std::shared_ptr<std::vector<unsigned int>> &                spResamplingParticleVector) const;

    /**
     *  @brief  Get the effective sample size
     *
     *  @return the effective sample size
     */
    double EffectiveSampleSize() const;

    /**
     *  @brief  Get the median step size
     *
     *  @param  observations the observations
     *
     *  @return the median step size
     */
    double GetMedianStepSize(ObservedStateVector observations) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double ObservedParticleState::GetResidualRange() const noexcept
{
    return m_residualRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ObservedParticleState::GetdEdx() const noexcept
{
    return m_dEdx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ObservedParticleState::Getdx() const noexcept
{
    return m_dx;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double MassHypothesis::Mass() const noexcept
{
    return m_mass;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::size_t MassHypothesis::NumberOfParticles() const noexcept
{
    return m_numberOfParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::vector<double> &MassHypothesis::PriorFinalEnergyDistribution() const noexcept
{
    return m_priorFinalEnergyDistribution;
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_FILTER_H
