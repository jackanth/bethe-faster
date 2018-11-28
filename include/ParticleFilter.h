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

#include "gsl/gsl_rng.h"

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
     *  @param  position the position (cm)
     *  @param  dEdx the dE/dx value (MeV/cm)
     *  @param  dx the change in position (cm)
     */
    ObservedParticleState(const double position, const double dEdx, const double dx) noexcept;

    /**
     *  @brief  Get the position (cm)
     *
     *  @return the position
     */
    double GetPosition() const noexcept;

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
    double m_position; ///< The position (cm)
    double m_dEdx;     ///< The dE/dx value (MeV/cm)
    double m_dx;       ///< The dx value (cm)
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
     *  @param  priorEnergyDistribution the prior energy distribution
     */
    MassHypothesis(const double mass, std::vector<double> priorEnergyDistribution);

    /**
     *  @brief  Get the mass
     *
     *  @return the mass
     */
    double Mass() const noexcept;

    /**
     *  @brief  Get the prior energy distribution
     *
     *  @return the prior energy distribution
     */
    const std::vector<double> &PriorEnergyDistribution() const noexcept;

private:
    std::shared_ptr<Propagator> m_spPropagator;            ///< The propagator
    double                      m_mass;                    ///< The particle mass
    std::vector<double>         m_priorEnergyDistribution; ///< The prior energy distribution
};

/**
 *  @brief  FilterOptions struct
 */
struct FilterOptions
{
    /**
     *  @brief  Constructor
     */
    FilterOptions() noexcept;

    double      m_stepSize;            ///< The step size (cm)
    std::size_t m_nParticles;          ///< The number of particles
    std::size_t m_maxUsedObservations; ///< The maximum number of observations to use
};

/**
 *  @brief  ParticleFilter class
 */
class ParticleFilter
{
public:
    using WeightedParticle     = std::pair<double, std::shared_ptr<Particle>>;       ///< Alias for a particle with a weight
    using ParticleDistribution = std::vector<WeightedParticle>;                      ///< Alias for a particle distribution
    using ObservedStateVector  = std::vector<ObservedParticleState>;                 ///< Alias for a vector of observed particle states
    using DistributionHistory  = std::vector<std::shared_ptr<ParticleDistribution>>; ///< Alias for a history of distributions

    /**
     *  @brief  Constructor
     *
     *  @param  spPropagator shared pointer to the propagator
     *  @param  options the options
     */
    ParticleFilter(std::shared_ptr<Propagator> spPropagator, FilterOptions options = FilterOptions());

    /**
     *  @brief  Destructor
     */
    ~ParticleFilter();

    /**
     *  @brief  Filter a uniform energy distribution on a particle history
     *
     *  @param  history the particle histroy
     *  @param  spMassHypothesis shared pointer to the mass hypothesis
     *
     *  @return the distribution history
     */
    DistributionHistory Filter(const Particle::History &history, const std::shared_ptr<MassHypothesis> &spMassHypothesis) const;

    /**
     *  @brief  Filter a uniform energy distribution on a set of observations
     *
     *  @param  observations the ordered vector of observations
     *  @param  spMassHypothesis shared pointer to the mass hypothesis
     *
     *  @return the distribution history
     */
    DistributionHistory Filter(const ObservedStateVector &observations, const std::shared_ptr<MassHypothesis> &spMassHypothesis) const;

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
    std::shared_ptr<Propagator>           m_spPropagator;                ///< The propagator
    std::shared_ptr<ParticleDistribution> m_spDistribution;              ///< The current distribution
    gsl_rng *                             m_pGenerator;                  ///< Address of the GSL random number generator
    FilterOptions                         m_options;                     ///< The options
    double *                              m_pResamplingProbabilityArray; ///< The resampling probability array
    unsigned int *                        m_pResamplingParticleArray;    ///< The resampling particle array

    /**
     *  @brief  Filter on an observation
     *
     *  @param  observation the observation
     *  @param  isFirstObservation whether this is the first observation
     *  @param  nSteps the number of propagation steps
     */
    bool FilterOnObservation(const ObservedParticleState &observation, const bool isFirstObservation, const std::size_t nSteps) const;

    /**
     *  @brief  Check the distribution weights
     *
     *  @return Whether the weights valid
     */
    bool CheckDistributionWeights() const;

    /**
     *  @brief  Check the legality of the filter options
     *
     *  @param  spMassHypothesis shared pointer to the mass hypothesis
     */
    void CheckOptions(const std::shared_ptr<MassHypothesis> &spMassHypothesis) const;

    /**
     *  @brief  Create a uniform prior over the initial energy bounds
     *
     *  @param  spMassHypothesis shared pointer to the mass hypothesis
     *  @param  initialPosition the initial position
     */
    void CreateDistribution(const std::shared_ptr<MassHypothesis> &spMassHypothesis, const double initialPosition) const;

    /**
     *  @brief  Propagate particles in the distribution
     *
     *  @param  nSteps the number of steps
     */
    void PropagateParticles(const std::size_t nSteps) const;

    /**
     *  @brief  Resample the distribution
     */
    void ResampleDistribution() const;

    /**
     *  @brief  Get the effective sample size
     *
     *  @return the effective sample size
     */
    double EffectiveSampleSize() const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double ObservedParticleState::GetPosition() const noexcept
{
    return m_position;
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

inline const std::vector<double> &MassHypothesis::PriorEnergyDistribution() const noexcept
{
    return m_priorEnergyDistribution;
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_FILTER_H
