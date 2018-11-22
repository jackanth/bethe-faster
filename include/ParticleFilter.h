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

#include <vector>
#include <memory>

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
     *  @param  dEdx the observed dE/dx value (MeV/cm)
     */
    ObservedParticleState(const float position, const float dEdx) noexcept;

    /**
     *  @brief  Get the position (cm)
     * 
     *  @return the position
     */
    float GetPosition() const noexcept;

    /**
     *  @brief  Get the dE/dx value (MeV/cm)
     * 
     *  @return the dE/dx value
     */
    float GetdEdx() const noexcept;

private:
    float m_position; ///< The position (cm)
    float m_dEdx;     ///< The dE/dx value (MeV/cm)
};

/**
 *  @brief  ParticleFilter class
 */
class ParticleFilter
{
public:
    using WeightedParticle = std::pair<float, std::shared_ptr<Particle>>; ///< Alias for a particle with a weight
    using ParticleDistribution = std::vector<WeightedParticle>;           ///< Alias for a particle distribution

    /**
     *  @brief  Constructor
     * 
     *  @param  detector the detector
     *  @param  initialDistribution the initial particle distribution
     *  @param  the position step size
     */
    ParticleFilter(const Detector &detector, ParticleDistribution initialDistribution, const float stepSize);

    /**
     *  @brief  Destructor
     */
    ~ParticleFilter();

    /**
     *  @brief  Filter the distribution on an observation
     * 
     *  @param  initialDistribution the initial particle distribution
     */
    const ParticleDistribution &Filter(const ObservedParticleState &observedState);

private:
    std::shared_ptr<Propagator> m_spPropagator;    ///< The propagator
    float                       m_stepSize;        ///< The step size
    ParticleDistribution        m_distribution;    ///< The observed particle state vector
    float                       m_position;        ///< The current position
    float                       m_initialPosition; ///< The initial position
    std::size_t                 m_nParticles;      ///< The number of particles
    gsl_rng                   * m_pGenerator;      ///< Address of the GSL random number generator
    float                       m_shrinkageFactor; ///< The filter shrinkage factor
    bool                        m_firstObservation; ///< Whether this is the first observation

    /**
     *  @brief  Normalize the distribution
     */
    void NormalizeDistribution();
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float ObservedParticleState::GetPosition() const noexcept
{
    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ObservedParticleState::GetdEdx() const noexcept
{
    return m_dEdx;
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_FILTER_H
