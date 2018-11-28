/**
 *  @file   bethe-faster/src/ParticleHelper.cc
 *
 *  @brief  Implementation of the particle helper class.
 *
 *  $Log: $
 */

#include "ParticleHelper.h"
#include "Particle.h"

#include <cmath>

namespace bf
{

std::shared_ptr<Particle> ParticleHelper::GetParticle(const double mass, const double initialKineticEnergy, const double initialPosition)
{
    return std::shared_ptr<Particle>{new Particle{mass, initialKineticEnergy, initialPosition}};
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ParticleHelper::GetParticleBeta(const std::shared_ptr<Particle> &spParticle)
{
    return ParticleHelper::GetParticleBeta(spParticle->Mass(), spParticle->KineticEnergy());
}

} // namespace bf
