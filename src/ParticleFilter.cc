/**
 *  @file   bethe-faster/src/ParticleFilter.cc
 *
 *  @brief  Implementation of the particle filter class.
 *
 *  $Log: $
 */

#include "ParticleFilter.h"

#include "gsl/gsl_randist.h"
#include "gsl/gsl_linalg.h"

#include <iostream>
#include <ctime>

namespace bf
{

ObservedParticleState::ObservedParticleState(const float position, const float dEdx) noexcept :
    m_position{position},
    m_dEdx{dEdx}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::ParticleFilter(const Detector &detector, ParticleDistribution initialDistribution, const float stepSize) :
    m_spPropagator{new Propagator{detector}},
    m_distribution{std::move_if_noexcept(initialDistribution)},
    m_stepSize{stepSize},
    m_position{0.f},
    m_initialPosition{0.f},
    m_nParticles{m_distribution.size()},
    m_pGenerator{nullptr},
    m_shrinkageFactor{0.01f},
    m_firstObservation{true}
{
    if (m_nParticles == 0UL)
        throw std::runtime_error{"Cannot filter with an empty prior distribution"};

    bool firstParticle = true;

    for (const auto &pair : m_distribution)
    {
        if (firstParticle)
        {
            m_initialPosition = pair.second->Position();
            m_position = m_initialPosition;
            firstParticle = false;
        }

        else if (m_position != pair.second->Position())
            throw std::runtime_error{"The distribution particles had different starting positions"};
    }

    this->NormalizeDistribution();

    // Set up the random number generator
    const gsl_rng_type * pType;
    gsl_rng_env_setup();
    pType = gsl_rng_default;
    m_pGenerator = gsl_rng_alloc(pType);
    gsl_rng_set(m_pGenerator, std::time(nullptr));
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::~ParticleFilter()
{
    gsl_rng_free(m_pGenerator);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFilter::ParticleDistribution &ParticleFilter::Filter(const ObservedParticleState &observedState)
{
    if (!m_firstObservation)
    {
        double *probabilityArray = new double[m_nParticles];
        unsigned int *particleArray = new unsigned int[m_nParticles];

        for (std::size_t i = 0UL; i < m_nParticles; ++i)
        {
            const auto &pair = m_distribution.at(i);
            probabilityArray[i] = pair.first;
            particleArray[i] = 0UL;
        }

        gsl_ran_multinomial(m_pGenerator, m_nParticles, m_nParticles, probabilityArray, particleArray);
        ParticleDistribution newDistribution;

        for (std::size_t i = 0UL; i < m_nParticles; ++i)
        {
            const std::shared_ptr<Particle> &spParticle = m_distribution.at(i).second;

            for (std::size_t j = 0UL; j < particleArray[i]; ++j)
            {
                const auto spNewParticle = std::shared_ptr<Particle>{new Particle{*spParticle}};
                newDistribution.emplace_back(1.f, spNewParticle);
            }
        }

        delete[] probabilityArray;
        delete[] particleArray;

        m_distribution = std::move(newDistribution);
        this->NormalizeDistribution();
    }

    if (observedState.GetPosition() < m_position)
        throw std::runtime_error{"Cannot filter on out-of-order observations"};

    m_position += m_stepSize;
    float weightSum = 0.f;

    for (auto &pair : m_distribution)
    {
        if (pair.first == 0.f)
            continue;

        try
        {   
            while (pair.second->Position() < observedState.GetPosition() && pair.second->KineticEnergy() > 0.f)
                m_spPropagator->Propagate(pair.second, m_stepSize);
        }

        catch (std::runtime_error &err)
        {
            std::cerr << "Filter was forced to reject particle: " << err.what() << std::endl;
            pair.first = 0.f;
            continue;
        }

        if (pair.second->KineticEnergy() == 0.f)
        {
            pair.first = 0.f;
            continue;
        }

        if (m_firstObservation)
            pair.first = m_spPropagator->CalculateObservationProbability(pair.second, observedState.GetdEdx(), m_stepSize);

        else
        {
            const float lastKineticEnergy = pair.second->LastKineticEnergy();
            const float currentKineticEnergy = pair.second->KineticEnergy();
            const float observedDeltaE = observedState.GetdEdx() * m_stepSize;
            const float energyResidual = currentKineticEnergy - lastKineticEnergy - observedDeltaE;

            const float pObservation = m_spPropagator->CalculateObservationProbability(pair.second, observedState.GetdEdx(), m_stepSize);
            const float pTransition = m_spPropagator->CalculateTransitionProbability(currentKineticEnergy - lastKineticEnergy,  pair.second->Mass(), lastKineticEnergy);
            const float pTransitionGivenObs = gsl_ran_gaussian_pdf(energyResidual, observedDeltaE / 2.f);
            
            pair.first = pObservation * pTransition * pTransitionGivenObs;
        }

        weightSum += pair.first;
    }

    if (weightSum < std::numeric_limits<float>::epsilon())
        throw std::runtime_error{"All particles stopped before observations completed"};

    this->NormalizeDistribution();

    if (m_firstObservation)
        m_firstObservation = false;

    return m_distribution;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleFilter::NormalizeDistribution()
{
    float totalWeight(0.f);

    for (const auto &pair : m_distribution)
        totalWeight += pair.first;

    if (totalWeight == 1.f)
        return;

    for (auto &pair : m_distribution)
        pair.first /= totalWeight;
}

} // namespace bf
