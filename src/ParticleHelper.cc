/**
 *  @file   bethe-faster/src/ParticleHelper.cc
 *
 *  @brief  Implementation of the particle helper class.
 *
 *  $Log: $
 */

#include "ParticleHelper.h"

namespace bf
{

std::string ParticleHelper::ToString(const PARTICLE_TYPE particleType)
{
    switch (particleType)
    {
        case PARTICLE_TYPE::MUON:
            return "MUON";

        case PARTICLE_TYPE::CHARGED_PION:
            return "CHARGED_PION";

        case PARTICLE_TYPE::CHARGED_KAON:
            return "CHARGED_KAON";

        case PARTICLE_TYPE::PROTON:
            return "PROTON";

        case PARTICLE_TYPE::OTHER:
            return "OTHER";

        default:
            break;
    }

    throw std::runtime_error{"Unknown particle type"};
}

} // namespace bf
