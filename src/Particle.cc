/**
 *  @file   bethe-faster/src/Particle.cc
 *
 *  @brief  Implementation of the particle class.
 *
 *  $Log: $
 */

#include "Particle.h"
#include "PhysicalConstants.h"

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

Particle::Particle(const double mass, const double finalKineticEnergy, const double finalResidualRange, const bool recordHistory) :
    m_minKineticEnergy{0.},
    m_maxKineticEnergy{0.},
    m_finalKineticEnergy{finalKineticEnergy},
    m_finalResidualRange{finalResidualRange},
    m_mass(mass),
    m_kineticEnergy(finalKineticEnergy),
    m_residualRange(finalResidualRange),
    m_recordHistory(recordHistory),
    m_hasFailed(false)
{
    if (m_mass < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize particle with a very small or negative mass");

    m_minKineticEnergy = this->CalculateEnergyFromBetaGamma(m_mass, PhysicalConstants::m_betheMinBetaGamma);
    m_maxKineticEnergy = this->CalculateEnergyFromBetaGamma(m_mass, PhysicalConstants::m_betheMaxBetaGamma);

    if (m_kineticEnergy < m_minKineticEnergy)
    {
        m_kineticEnergy = m_minKineticEnergy;
        m_finalKineticEnergy = m_minKineticEnergy;
    }

    if (m_kineticEnergy > m_maxKineticEnergy)
    {
        m_kineticEnergy = m_maxKineticEnergy;
        m_finalKineticEnergy = m_maxKineticEnergy;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Particle::Increment(const double deltaPosition, const double deltaEnergy)
{
    if (m_hasFailed)
        return;

    if (deltaPosition < -std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"Cannot decrement residual range"};

    if (deltaEnergy < -std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"Cannot decrement kinetic energy"};

    const double newEnergy = m_kineticEnergy + deltaEnergy;

    if (newEnergy > m_maxKineticEnergy)
    {
        m_hasFailed = true;
        return;
    }

    if (m_recordHistory)
    {
        const double dEdx = -deltaEnergy / deltaPosition;
        m_history.emplace_back(std::make_shared<ParticleState>(ParticleState{m_residualRange, m_kineticEnergy, deltaPosition, dEdx})); // add the current state to the history
    }

    m_kineticEnergy = newEnergy;
    m_residualRange += deltaPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Particle::Reset() noexcept
{
    m_residualRange = m_finalResidualRange;
    m_kineticEnergy = m_finalKineticEnergy;
    m_hasFailed     = false;
    m_history.clear();
}

} // namespace bf
