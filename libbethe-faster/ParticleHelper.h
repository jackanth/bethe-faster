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

#include <cmath>
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
     *  @brief  Enumeration of the types of particle
     */
    enum class PARTICLE_TYPE
    {
        MUON,         ///< A muon
        PROTON,       ///< A proton
        CHARGED_PION, ///< A charged pion
        CHARGED_KAON, ///< A charged kaon
        OTHER         ///< Another particle type
    };

    /**
     *  @brief  Deleted constructor
     */
    ParticleHelper() = delete;

    /**
     *  @brief  Get a particle
     *
     *  @param  mass the particle mass
     *  @param  kineticEnergy the kinetic energy
     *  @param  residualRange the residual range
     *
     *  @return shared pointer to the particle
     */
    static std::shared_ptr<Particle> GetParticle(const double mass, const double kineticEnergy, const double residualRange = 0.);

    /**
     *  @brief  Get a muon
     *
     *  @param  kineticEnergy the kinetic energy
     *  @param  residualRange the residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetMuon(const double kineticEnergy = 0., const double residualRange = 0.);

    /**
     *  @brief  Get a proton
     *
     *  @param  kineticEnergy the kinetic energy
     *  @param  residualRange the residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetProton(const double kineticEnergy = 0., const double residualRange = 0.);

    /**
     *  @brief  Get a charged pion
     *
     *  @param  kineticEnergy the kinetic energy
     *  @param  residualRange the residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetChargedPion(const double kineticEnergy = 0., const double residualRange = 0.);

    /**
     *  @brief  Get a charged kaon
     *
     *  @param  kineticEnergy the kinetic energy
     *  @param  residualRange the residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetChargedKaon(const double kineticEnergy = 0., const double residualRange = 0.);

    /**
     *  @brief  Calculate gamma for a particle
     *
     *  @param  spParticle shared pointer to the particle
     *
     *  @return the particle gamma
     */
    static double GetParticleGamma(const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Calculate gamma for a particle
     *
     *  @param  mass the particle mass
     *  @param  energy the particle energy
     *
     *  @return the particle gamma
     */
    static double GetParticleGamma(const double mass, const double energy);

    /**
     *  @brief  Calculate beta for a particle
     *
     *  @param  spParticle shared pointer to the particle
     *
     *  @return the particle beta
     */
    static double GetParticleBeta(const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Calculate beta for a particle
     *
     *  @param  mass the particle mass
     *  @param  energy the particle energy
     *
     *  @return the particle beta
     */
    static double GetParticleBeta(const double mass, const double energy);

    /**
     *  @brief  Turn a particle type enum into a string
     *
     *  @param  particleType the particle type
     *
     *  @return the particle type as a string
     */
    static std::string ToString(const PARTICLE_TYPE particleType);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetMuon(const double kineticEnergy, const double residualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_muonMass, kineticEnergy, residualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetProton(const double kineticEnergy, const double residualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_protonMass, kineticEnergy, residualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetChargedPion(const double kineticEnergy, const double residualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedPionMass, kineticEnergy, residualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetChargedKaon(const double kineticEnergy, const double residualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedKaonMass, kineticEnergy, residualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetParticleGamma(const double mass, const double energy)
{
    return energy / mass + 1.;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetParticleGamma(const std::shared_ptr<Particle> &spParticle)
{
    return ParticleHelper::GetParticleGamma(spParticle->Mass(), spParticle->KineticEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetParticleBeta(const double mass, const double energy)
{
    const double gamma = ParticleHelper::GetParticleGamma(mass, energy);
    return std::sqrt(1. - 1. / (gamma * gamma));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetParticle(const double mass, const double kineticEnergy, const double residualRange)
{
    return std::shared_ptr<Particle>{new Particle{mass, kineticEnergy, residualRange}};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetParticleBeta(const std::shared_ptr<Particle> &spParticle)
{
    return ParticleHelper::GetParticleBeta(spParticle->Mass(), spParticle->KineticEnergy());
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_HELPER_H
