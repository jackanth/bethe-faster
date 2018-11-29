/**
 *  @file   bethe-faster/src/Particle.cc
 *
 *  @brief  Implementation of the particle class.
 *
 *  $Log: $
 */

#include "Particle.h"

#include <iostream>
#include <limits>

namespace bf
{

ParticleState::ParticleState(const double position, const double kineticEnergy, const double dx, const double dEdx) noexcept :
    m_position{position},
    m_kineticEnergy{kineticEnergy},
    m_dx{dx},
    m_dEdx{dEdx}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

Particle::Particle(const double mass, const double initialKineticEnergy, const double initialPosition, const bool recordHistory) :
    m_mass(mass),
    m_kineticEnergy(initialKineticEnergy),
    m_position(initialPosition),
    m_recordHistory(recordHistory),
    m_isAlive{true}
{
    if (m_mass < 0.)
        throw std::runtime_error("Could not initialize particle with negative mass");

    if (m_kineticEnergy < 0.)
        throw std::runtime_error("Could not initialize particle with negative energy");
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Particle::LastKineticEnergy() const
{
    if (m_history.empty())
        throw std::runtime_error{"Could not get last kinetic energy for particle because it had no history"};

    return m_history.back()->GetKineticEnergy();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool Particle::Increment(const double deltaPosition, const double deltaEnergy)
{
    if (m_kineticEnergy <= std::numeric_limits<double>::epsilon())
    {
        m_isAlive = false;
        return false;
    }

    const double newKineticEnergy = m_kineticEnergy + deltaEnergy;

    if (newKineticEnergy <= std::numeric_limits<double>::epsilon())
    {
        m_kineticEnergy = 0.;
        m_isAlive = false;
        return false;
    }

    if (m_recordHistory)
    {
        const double dEdx = deltaEnergy / deltaPosition;
        m_history.emplace_back(std::make_shared<ParticleState>(ParticleState{m_position, m_kineticEnergy, deltaPosition, dEdx})); // add the current state to the history
    }

    m_kineticEnergy = newKineticEnergy;
    m_position += deltaPosition;

    return true;
}

} // namespace bf
