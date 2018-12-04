/**
 *  @file   bethe-faster/include/Particle.h
 *
 *  @brief  Header file for the particle class.
 *
 *  $Log: $
 */

#ifndef BF_PARTICLE_H
#define BF_PARTICLE_H 1

#include <memory>
#include <vector>
#include <cmath>

namespace bf
{

/**
 *  @brief  Forward declaration of the Propagator class
 */
class Propagator;

/**
 *  @brief  Forward declaration of the Particle class
 */
class Particle;

/**
 *  @brief  Forward declaration of the ParticleFilter class
 */
class ParticleFilter;

/**
 *  @brief  ParticleState class
 */
class ParticleState
{
public:
    /**
     *  @brief  Get the residual range (cm)
     *
     *  @return the residual range
     */
    double GetResidualRange() const noexcept;

    /**
     *  @brief  Get the kinetic energy (MeV)
     *
     *  @return the kinetic energy
     */
    double GetKineticEnergy() const noexcept;

    /**
     *  @brief  Get the dx value (cm)
     *
     *  @return the dx value
     */
    double Getdx() const noexcept;

    /**
     *  @brief  Get the dE/dx value (MeV/cm)
     *
     *  @return the dE/dx value
     */
    double GetdEdx() const noexcept;

protected:
    /**
     *  @brief  Constructor
     *
     *  @param  residualRange the residual range (cm)
     *  @param  kineticEnergy the kinetic energy (MeV)
     *  @param  dx the delta position (cm)
     *  @param  dEdx the dE/dx value (MeV/cm)
     */
    ParticleState(const double residualRange, const double kineticEnergy, const double dx, const double dEdx) noexcept;

    friend class Particle;

private:
    double m_residualRange; ///< The residualRange (cm)
    double m_kineticEnergy; ///< The kinetic energy (MeV)
    double m_dx;            ///< The delta position
    double m_dEdx;          ///< The dE/dx value (MeV/cm)
};

/**
 *  @brief  Particle class
 */
class Particle
{
public:
    using History = std::vector<std::shared_ptr<ParticleState>>; ///< An alias for the history of particle states

    /**
     *  @brief  Get the mass (MeV)
     *
     *  @return the mass
     */
    double Mass() const noexcept;

    /**
     *  @brief  Get the particle history
     *
     *  @return the particle history
     */
    const History &GetHistory() const;

    /**
     *  @brief  Get the minimum allowed kinetic energy
     *
     *  @return the minimum kinetic energy
     */
    double MinKineticEnergy() const noexcept;

    /**
     *  @brief  Get the maximum allowed kinetic energy
     *
     *  @return the maximum kinetic energy
     */
    double MaxKineticEnergy() const noexcept;

    /**
     *  @brief  Get the current kinetic energy (MeV)
     *
     *  @return the current kinetic energy
     */
    double KineticEnergy() const noexcept;

    /**
     *  @brief  Get the current residual range (cm)
     *
     *  @return the current residual range
     */
    double ResidualRange() const noexcept;

    /**
     *  @brief  Get whether the particle has failed in its propagation
     *
     *  @return whether the particle has failed
     */
    bool HasFailed() const noexcept;

    /**
     *  @brief  Reset the particle to its final state
     */
    void Reset() noexcept;

protected:
    /**
     *  @brief  Constructor
     *
     *  @param  mass the particle mass (MeV)
     *  @param  finalKineticEnergy the final kinetic energy (MeV)
     *  @param  finalResidualRange the final residual range (cm)
     *  @param  recordHistory whether to record the particle history
     */
    Particle(const double mass, const double finalKineticEnergy, const double finalResidualRange = 0., const bool recordHistory = true);

    /**
     *  @brief  Increment the particle
     *
     *  @param  deltaPosition the position value by which to increment (cm)
     *  @param  deltaEnergy the kinetic energy by which to increment (MeV)
     */
    void Increment(const double deltaPosition, const double deltaEnergy);

    /**
     *  @brief  Set the kinetic energy
     *
     *  @param  kineticEnergy the kinetic energy (MeV)
     */
    void SetKineticEnergy(const double kineticEnergy) noexcept;

    /**
     *  @brief  Set whether the particle has failed
     *
     *  @param  hasFailed whether the particle has failed
     */
    void HasFailed(const bool hasFailed) noexcept;

    friend class Propagator;
    friend class ParticleFilter;
    friend class ParticleHelper;

private:
    double  m_minKineticEnergy;   ///< The minimum allowed energy (MeV)
    double  m_maxKineticEnergy;   ///< The maximum allowed energy (MeV)
    double  m_finalKineticEnergy; ///< The final kinetic energy (MeV)
    double  m_finalResidualRange; ///< The final residual range (cm)
    double  m_mass;               ///< The particle mass (MeV)
    double  m_kineticEnergy;      ///< The kinetic energy (MeV)
    double  m_residualRange;      ///< The current residual range (cm)
    History m_history;            ///< The particle history
    bool    m_recordHistory;      ///< Whether to record the particle history
    bool    m_hasFailed;          ///< Whether the particle has failed in its propagation

    /**
     *  @brief  Calculate a particle's kinetic energy from its beta*gamma value
     *
     *  @param  mass the particle's mass
     *  @param  betaGamma the beta*gamma value
     *
     *  @return the kinetic energy
     */
    double CalculateEnergyFromBetaGamma(const double mass, const double betaGamma) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleState::GetResidualRange() const noexcept
{
    return m_residualRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleState::GetKineticEnergy() const noexcept
{
    return m_kineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleState::Getdx() const noexcept
{
    return m_dx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleState::GetdEdx() const noexcept
{
    return m_dEdx;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double Particle::Mass() const noexcept
{
    return m_mass;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const Particle::History &Particle::GetHistory() const
{
    if (!m_recordHistory)
        throw std::runtime_error{"Could not get history for particle because it was not recorded"};

    return m_history;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Particle::MinKineticEnergy() const noexcept
{
    return m_minKineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Particle::MaxKineticEnergy() const noexcept
{
    return m_maxKineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Particle::KineticEnergy() const noexcept
{
    return m_kineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Particle::ResidualRange() const noexcept
{
    return m_residualRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool Particle::HasFailed() const noexcept
{
    return m_hasFailed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Particle::SetKineticEnergy(const double kineticEnergy) noexcept
{
    m_kineticEnergy = kineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Particle::HasFailed(const bool hasFailed) noexcept
{
    m_hasFailed = hasFailed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Particle::CalculateEnergyFromBetaGamma(const double mass, const double betaGamma) const
{
    return mass * (std::sqrt(1. + betaGamma * betaGamma) - 1.);
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_H
