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

ParticleState::ParticleState(const double residualRange, const double kineticEnergy, const double dx, const double dEdx) noexcept :
    m_residualRange{residualRange},
    m_kineticEnergy{kineticEnergy},
    m_dx{dx},
    m_dEdx{dEdx}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

Particle::Particle(const double mass, const double finalKineticEnergy, const double initialResidualRange, const bool recordHistory) :
    m_mass(mass),
    m_kineticEnergy(finalKineticEnergy),
    m_residualRange(initialResidualRange),
    m_recordHistory(recordHistory)
{
    if (m_mass < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize particle with a very small or negative mass");

    if (m_kineticEnergy < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize particle with very small or negative energy");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Particle::Increment(const double deltaPosition, const double deltaEnergy)
{
    if (m_hasFailed)
        return;
        
    if (deltaPosition < 0.)
        throw std::runtime_error{"Cannot decrement position"};

    if (deltaEnergy < 0.)
        throw std::runtime_error{"Cannot decrement energy"};

    if (m_recordHistory)
    {
        const double dEdx = -deltaEnergy / deltaPosition;
        m_history.emplace_back(std::make_shared<ParticleState>(ParticleState{m_residualRange, m_kineticEnergy, deltaPosition, dEdx})); // add the current state to the history
    }

    m_kineticEnergy += deltaEnergy;
    m_residualRange += deltaPosition;
}

} // namespace bf
