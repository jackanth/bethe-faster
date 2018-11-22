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

std::shared_ptr<Particle> ParticleHelper::GetParticle(const float mass, const float initialKineticEnergy, const float initialPosition)
{
    return std::shared_ptr<Particle>{new Particle{mass, initialKineticEnergy, initialPosition}};
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ParticleHelper::GetParticleBeta(const std::shared_ptr<Particle> &spParticle)
{
    return ParticleHelper::GetParticleBeta(spParticle->Mass(), spParticle->KineticEnergy());
}

} // namespace bf
