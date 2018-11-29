/**
 *  @file   bethe-faster/include/Particle.h
 *
 *  @brief  Header file for the particle class.
 *
 *  $Log: $
 */
#ifndef BF_PARTICLE_H
#define BF_PARTICLE_H 1

#include <vector>
#include <memory>

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
     *  @brief  Get the position (cm)
     *
     *  @return the position
     */
    double GetPosition() const noexcept;

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
     *  @param  position the position (cm)
     *  @param  kineticEnergy the kinetic energy (MeV)
     *  @param  dx the delta position (cm)
     *  @param  dEdx the dE/dx value (MeV/cm)
     */
    ParticleState(const double position, const double kineticEnergy, const double dx, const double dEdx) noexcept;

    friend class Particle;

private:
    double m_position;      ///< The position (cm)
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
     *  @brief  Constructor
     *
     *  @param  mass the particle mass (MeV)
     *  @param  initialKineticEnergy the initial kinetic energy (MeV)
     *  @param  initialPosition the initial position (cm)
     *  @param  recordHistory whether to record the particle history
     */
    Particle(const double mass, const double initialKineticEnergy, const double initialPosition = 0., const bool recordHistory = true);

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
     *  @brief  Get the current kinetic energy (MeV)
     *
     *  @return the current kinetic energy
     */
    double KineticEnergy() const noexcept;

    /**
     *  @brief  Get the current position (cm)
     *
     *  @return the current position
     */
    double Position() const noexcept;

    /**
     *  @brief  Get the last kinetic energy (MeV)
     *
     *  @return the last kinetic energy
     */
    double LastKineticEnergy() const;

    /**
     *  @brief  Get whether the particle is alive
     *
     *  @return whether the particle is alive
     */
    bool IsAlive() const noexcept;

protected:
    /**
     *  @brief  Increment the particle
     *
     *  @param  deltaPosition the position value by which to increment (cm)
     *  @param  deltaEnergy the kinetic energy by which to increment (MeV)
     *
     *  @return whether beta is still nonzero
     */
    bool Increment(const double deltaPosition, const double deltaEnergy);

    /**
     *  @brief  Set the kinetic energy
     *
     *  @param  kineticEnergy the new kinetic energy
     */
    void SetKineticEnergy(const double kineticEnergy) noexcept;

    /**
     *  @brief  Kill the particle
     */
    void Kill() noexcept;

    friend class Propagator;
    friend class ParticleFilter;

private:
    double  m_mass;          ///< The particle mass (MeV)
    double  m_kineticEnergy; ///< The kinetic energy (MeV)
    double  m_position;      ///< The current position (cm)
    History m_history;       ///< The particle history
    bool    m_recordHistory; ///< Whether to record the particle history
    bool    m_isAlive;       ///< Whether the particle is alive
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleState::GetPosition() const noexcept
{
    return m_position;
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

inline double Particle::KineticEnergy() const noexcept
{
    return m_kineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double Particle::Position() const noexcept
{
    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool Particle::IsAlive() const noexcept
{
    return m_isAlive;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Particle::SetKineticEnergy(const double kineticEnergy) noexcept
{
    m_kineticEnergy = kineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Particle::Kill() noexcept
{
    m_isAlive = false;
    m_kineticEnergy = 0.;
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_H
