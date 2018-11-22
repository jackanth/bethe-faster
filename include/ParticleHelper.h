/**
 *  @file   bethe-faster/include/ParticleHelper.h
 *
 *  @brief  Header file for the particle helper class.
 *
 *  $Log: $
 */
#ifndef BF_PARTICLE_HELPER_H
#define BF_PARTICLE_HELPER_H 1

#include "Particle.h"
#include "PhysicalConstants.h"

#include <memory>
#include <set>

namespace bf
{

/**
 *  @brief  Forward declaration of Particle class
 */
class Particle;

/**
 *  @brief  ParticleHelper class
 */
class ParticleHelper
{
public:
    using ParticleSet = std::set<std::shared_ptr<Particle>, std::owner_less<std::shared_ptr<Particle>>>; ///< Alias for a set of particle shared pointers

    /**
     *  @brief  Deleted constructor
     */
    ParticleHelper() = delete;

    /**
     *  @brief  Get a particle
     * 
     *  @param  mass the particle mass
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the particle
     */
    static std::shared_ptr<Particle> GetParticle(const float mass, const float initialKineticEnergy, const float initialPosition = 0.f);

    /**
     *  @brief  Get a muon
     * 
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetMuon(const float initialKineticEnergy, const float initialPosition = 0.f);

    /**
     *  @brief  Get a proton
     * 
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetProton(const float initialKineticEnergy, const float initialPosition = 0.f);

    /**
     *  @brief  Get a charged pion
     * 
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetChargedPion(const float initialKineticEnergy, const float initialPosition = 0.f);

    /**
     *  @brief  Get a charged kaon
     * 
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetChargedKaon(const float initialKineticEnergy, const float initialPosition = 0.f);

    /**
     *  @brief  Propagate a particle until it stops
     * 
     *  @param  spParticle shared pointer to the particle
     */
    static void PropagateUntilStopped(const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Get the beta value for a particle
     * 
     *  @param  spParticle shared pointer to the particle
     * 
     *  @return the particle beta value
     */
    static float GetParticleBeta(const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Get the beta value for a particle
     * 
     *  @param  mass the particle mass
     *  @param  energy the particle energy
     * 
     *  @return the particle beta value
     */
    static float GetParticleBeta(const float mass, const float energy);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetMuon(const float initialKineticEnergy, const float initialPosition)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_muonMass, initialKineticEnergy, initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetProton(const float initialKineticEnergy, const float initialPosition)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_protonMass, initialKineticEnergy, initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetChargedPion(const float initialKineticEnergy, const float initialPosition)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedPionMass, initialKineticEnergy, initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetChargedKaon(const float initialKineticEnergy, const float initialPosition)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedKaonMass, initialKineticEnergy, initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ParticleHelper::GetParticleBeta(const float mass, const float energy)
{
    const float gamma = energy / mass + 1.f;
    return 1.f - 1.f / (gamma * gamma);
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_HELPER_H
