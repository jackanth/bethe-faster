/**
 *  @file   bethe-faster/src/Particle.cc
 *
 *  @brief  Implementation of the particle class.
 *
 *  $Log: $
 */

#include "Particle.h"

namespace bf
{

ParticleState::ParticleState(const float position, const float kineticEnergy, const float dx, const float dEdx) noexcept :
    m_position{position},
    m_kineticEnergy{kineticEnergy},
    m_dx{dx},
    m_dEdx{dEdx}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

Particle::Particle(const float mass, const float initialKineticEnergy, const float initialPosition) :
    m_mass(mass),
    m_kineticEnergy(initialKineticEnergy),
    m_position(initialPosition)
{
    if (m_mass < 0.f)
        throw std::runtime_error("Could not initialize particle with negative mass");

    if (m_kineticEnergy < 0.f)
        throw std::runtime_error("Could not initialize particle with negative energy");
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Particle::LastKineticEnergy() const
{
    if (m_history.empty())
        throw std::runtime_error{"Could not get last kinetic energy for particle because it had no history"};

    return m_history.back().GetKineticEnergy();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool Particle::Increment(const float deltaPosition, const float deltaEnergy)
{
    if (m_kineticEnergy == 0.f)
        return false;

    const float dEdx = deltaEnergy / deltaPosition;
    m_history.emplace_back(ParticleState{m_position, m_kineticEnergy, deltaPosition, dEdx}); // add the current state to the history

    const float newKineticEnergy = m_kineticEnergy + deltaEnergy;

    if (newKineticEnergy <= 0.f)
    {
        m_kineticEnergy = 0.f;
        return false;
    }

    m_kineticEnergy = newKineticEnergy;
    m_position += deltaPosition;

    return true;
}

} // namespace bf
