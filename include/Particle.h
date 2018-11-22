/**
 *  @file   bethe-faster/include/Particle.h
 *
 *  @brief  Header file for the particle class.
 *
 *  $Log: $
 */
#ifndef BF_PARTICLE_H
#define BF_PARTICLE_H 1

#include <tuple>
#include <utility>
#include <vector>

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
    float GetPosition() const noexcept;

    /**
     *  @brief  Get the kinetic energy (MeV)
     *
     *  @return the kinetic energy
     */
    float GetKineticEnergy() const noexcept;

    /**
     *  @brief  Get the dx value (cm)
     *
     *  @return the dx value
     */
    float Getdx() const noexcept;

    /**
     *  @brief  Get the dE/dx value (MeV/cm)
     *
     *  @return the dE/dx value
     */
    float GetdEdx() const noexcept;

protected:
    /**
     *  @brief  Constructor
     *
     *  @param  position the position (cm)
     *  @param  kineticEnergy the kinetic energy (MeV)
     *  @param  dx the delta position (cm)
     *  @param  dEdx the dE/dx value (MeV/cm)
     */
    ParticleState(const float position, const float kineticEnergy, const float dx, const float dEdx) noexcept;

    friend class Particle;

private:
    float m_position;      ///< The position (cm)
    float m_kineticEnergy; ///< The kinetic energy (MeV)
    float m_dx;            ///< The delta position
    float m_dEdx;          ///< The dE/dx value (MeV/cm)
};

/**
 *  @brief  Particle class
 */
class Particle
{
public:
    using History = std::vector<ParticleState>; ///< An alias for the history of particle states

    /**
     *  @brief  Constructor
     *
     *  @param  mass the particle mass (MeV)
     *  @param  initialKineticEnergy the initial kinetic energy (MeV)
     *  @param  initialPosition the initial position (cm)
     */
    Particle(const float mass, const float initialKineticEnergy, const float initialPosition = 0.f);

    /**
     *  @brief  Get the mass (MeV)
     *
     *  @return the mass
     */
    float Mass() const noexcept;

    /**
     *  @brief  Get the particle history
     *
     *  @return the particle history
     */
    const History &GetHistory() const noexcept;

    /**
     *  @brief  Get the current kinetic energy (MeV)
     *
     *  @return the current kinetic energy
     */
    float KineticEnergy() const noexcept;

    /**
     *  @brief  Get the current position (cm)
     *
     *  @return the current position
     */
    float Position() const noexcept;

    /**
     *  @brief  Get the last kinetic energy (MeV)
     *
     *  @return the last kinetic energy
     */
    float LastKineticEnergy() const;

protected:
    /**
     *  @brief  Increment the particle
     *
     *  @param  deltaPosition the position value by which to increment (cm)
     *  @param  deltaEnergy the kinetic energy by which to increment (MeV)
     *
     *  @return whether beta is still nonzero
     */
    bool Increment(const float deltaPosition, const float deltaEnergy);

    friend class Propagator;

private:
    float   m_mass;          ///< The particle mass (MeV)
    float   m_kineticEnergy; ///< The kinetic energy (MeV)
    float   m_position;      ///< The current position (cm)
    History m_history;       ///< The particle history
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float ParticleState::GetPosition() const noexcept
{
    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ParticleState::GetKineticEnergy() const noexcept
{
    return m_kineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ParticleState::Getdx() const noexcept
{
    return m_dx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ParticleState::GetdEdx() const noexcept
{
    return m_dEdx;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float Particle::Mass() const noexcept
{
    return m_mass;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const Particle::History &Particle::GetHistory() const noexcept
{
    return m_history;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Particle::KineticEnergy() const noexcept
{
    return m_kineticEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Particle::Position() const noexcept
{
    return m_position;
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_H
