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
     *  @brief  Get a particle to propagate backwards
     *
     *  @param  mass the particle mass
     *  @param  finalKineticEnergy the final kinetic energy
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the particle
     */
    static std::shared_ptr<Particle> GetParticle(const double mass, const double finalKineticEnergy, const double finalResidualRange = 0.);

    /**
     *  @brief  Get a muon to propagate backwards
     *
     *  @param  finalKineticEnergy the final kinetic energy
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetMuon(const double finalKineticEnergy, const double finalResidualRange = 0.);

    /**
     *  @brief  Get a stopped muon to propagate backwards
     *
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetStoppedMuon(const double finalResidualRange = 0.);

    /**
     *  @brief  Get a proton to propagate backwards
     *
     *  @param  finalKineticEnergy the final kinetic energy
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetProton(const double finalKineticEnergy, const double finalResidualRange = 0.);

    /**
     *  @brief  Get a stopped proton to propagate backwards
     *
     *  @param  finalKineticEnergy the final kinetic energy
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetStoppedProton(const double finalResidualRange = 0.);

    /**
     *  @brief  Get a charged pion to propagate backwards
     *
     *  @param  finalKineticEnergy the final kinetic energy
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetChargedPion(const double finalKineticEnergy, const double finalResidualRange = 0.);

    /**
     *  @brief  Get a stopped charged pion to propagate backwards
     *
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetStoppedChargedPion(const double finalResidualRange = 0.);

    /**
     *  @brief  Get a charged kaon to propagate backwards
     *
     *  @param  finalKineticEnergy the final kinetic energy
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetChargedKaon(const double finalKineticEnergy, const double finalResidualRange = 0.);

    /**
     *  @brief  Get a stopped charged kaon to propagate backwards
     *
     *  @param  finalResidualRange the final residual range
     *
     *  @return shared pointer to the muon
     */
    static std::shared_ptr<Particle> GetStoppedChargedKaon(const double finalResidualRange = 0.);

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
     *  @brief  Get the minimum muon energy
     *
     *  @return the minimum muon energy
     */
    static double GetMinimumMuonEnergy();

    /**
     *  @brief  Get the minimum charged pion energy
     *
     *  @return the minimum charged pion energy
     */
    static double GetMinimumChargedPionEnergy();

    /**
     *  @brief  Get the minimum charged kaon energy
     *
     *  @return the minimum charged kaon energy
     */
    static double GetMinimumChargedKaonEnergy();

    /**
     *  @brief  Get the minimum proton energy
     *
     *  @return the minimum proton energy
     */
    static double GetMinimumProtonEnergy();

    /**
     *  @brief  Turn a particle type enum into a string
     *
     *  @param  particleType the particle type
     *
     *  @return the particle type as a string
     */
    static std::string ToString(const PARTICLE_TYPE particleType);

private:
    static constexpr double m_minimumMuonEnergy        = 0.01;  ///< The minimum calculable muon energy (MeV)
    static constexpr double m_minimumChargedPionEnergy = 0.015; ///< The minimum calculable charged pion energy (MeV)
    static constexpr double m_minimumChargedKaonEnergy = 0.05;  ///< The minimum calculable charged kaon energy (MeV)
    static constexpr double m_minimumProtonEnergy      = 0.095; ///< The minimum calculable proton energy (MeV)
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetMuon(const double finalKineticEnergy, const double finalResidualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_muonMass, finalKineticEnergy, finalResidualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetStoppedMuon(const double finalResidualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_muonMass, m_minimumMuonEnergy, finalResidualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetProton(const double finalKineticEnergy, const double finalResidualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_protonMass, finalKineticEnergy, finalResidualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetStoppedProton(const double finalResidualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_protonMass, m_minimumProtonEnergy, finalResidualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetChargedPion(const double finalKineticEnergy, const double finalResidualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedPionMass, finalKineticEnergy, finalResidualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetStoppedChargedPion(const double finalResidualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedPionMass, m_minimumChargedPionEnergy, finalResidualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetChargedKaon(const double finalKineticEnergy, const double finalResidualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedKaonMass, finalKineticEnergy, finalResidualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetStoppedChargedKaon(const double finalResidualRange)
{
    return ParticleHelper::GetParticle(PhysicalConstants::m_chargedKaonMass, m_minimumChargedKaonEnergy, finalResidualRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetParticleBeta(const double mass, const double energy)
{
    const double gamma = energy / mass + 1.;
    return std::sqrt(1. - 1. / (gamma * gamma));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::shared_ptr<Particle> ParticleHelper::GetParticle(const double mass, const double finalKineticEnergy, const double finalResidualRange)
{
    return std::shared_ptr<Particle>{new Particle{mass, finalKineticEnergy, finalResidualRange}};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetParticleBeta(const std::shared_ptr<Particle> &spParticle)
{
    return ParticleHelper::GetParticleBeta(spParticle->Mass(), spParticle->KineticEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetMinimumMuonEnergy()
{
    return m_minimumMuonEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetMinimumChargedPionEnergy()
{
    return m_minimumChargedPionEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetMinimumChargedKaonEnergy()
{
    return m_minimumChargedKaonEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ParticleHelper::GetMinimumProtonEnergy()
{
    return m_minimumProtonEnergy;
}

} // namespace bf

#endif // #ifndef BF_PARTICLE_HELPER_H
