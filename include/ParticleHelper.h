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
#include <cmath>

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
    static std::shared_ptr<Particle> GetParticle(const double mass, const double initialKineticEnergy, const double initialPosition = 0.);

    /**
     *  @brief  Get a muon
     * 
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetMuon(const double initialKineticEnergy, const double initialPosition = 0.);

    /**
     *  @brief  Get a proton
     * 
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetProton(const double initialKineticEnergy, const double initialPosition = 0.);

    /**
     *  @brief  Get a charged pion
     * 
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetChargedPion(const double initialKineticEnergy, const double initialPosition = 0.);

    /**
     *  @brief  Get a charged kaon
     * 
     *  @param  initialKineticEnergy the initial kinetic energy
     *  @param  initialPosition the initial position
     * 
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetChargedKaon(const double initialKineticEnergy, const double initialPosition = 0.);

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
    static double GetParticleBeta(const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Get the beta value for a particle
     * 
     *  @param  mass the particle mass
     *  @param  energy the particle energy
     * 
     *  @return the particle beta value
     */
    static double GetParticleBeta(const double mass, const double energy);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetMuon(const double initialKineticEnergy, const double initialPosition)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_muonMass, initialKineticEnergy, initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetProton(const double initialKineticEnergy, const double initialPosition)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_protonMass, initialKineticEnergy, initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetChargedPion(const double initialKineticEnergy, const double initialPosition)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedPionMass, initialKineticEnergy, initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetChargedKaon(const double initialKineticEnergy, const double initialPosition)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedKaonMass, initialKineticEnergy, initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetParticleBeta(const double mass, const double energy)
{
    const double gamma = energy / mass + 1.;
    return std::sqrt(1. - 1. / (gamma * gamma));
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_HELPER_H
