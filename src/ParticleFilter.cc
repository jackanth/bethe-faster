/**
 *  @file   bethe-faster/src/ParticleFilter.cc
 *
 *  @brief  Implementation of the particle filter class.
 *
 *  $Log: $
 */

#include "ParticleFilter.h"

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_randist.h"

#include <ctime>
#include <iostream>

namespace bf
{

ObservedParticleState::ObservedParticleState(const double position, const double dEdx, const double dx) noexcept :
    m_position{position},
    m_dEdx{dEdx},
    m_dx{dx}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MassHypothesis::MassHypothesis(const double mass, std::vector<double> priorEnergyDistribution) :
    m_mass{mass},
    m_priorEnergyDistribution{std::move_if_noexcept(priorEnergyDistribution)}
{
    if (m_priorEnergyDistribution.empty())
        throw std::runtime_error{"Must provide a prior energy distribution"};
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FilterOptions::FilterOptions() noexcept :
    m_stepSize{0.03},
    m_nParticles{1000UL},
    m_maxUsedObservations{std::numeric_limits<std::size_t>::max()}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::ParticleFilter(std::shared_ptr<Propagator> spPropagator, FilterOptions options) :
    m_spPropagator{std::move(spPropagator)},
    m_spDistribution{new ParticleDistribution{}},
    m_pGenerator{nullptr},
    m_options{std::move(options)},
    m_pResamplingProbabilityArray{nullptr},
    m_pResamplingParticleArray{nullptr}
{
    // Set up the random number generator
    const gsl_rng_type *pType;
    gsl_rng_env_setup();
    pType        = gsl_rng_default;
    m_pGenerator = gsl_rng_alloc(pType);
    gsl_rng_set(m_pGenerator, static_cast<std::uint32_t>(std::time(nullptr)));

    // Check the options
    if (m_options.m_nParticles < 2UL)
        throw std::runtime_error{"Number of filter particles must be 2 or greater"};

    if (1. / static_cast<double>(m_options.m_nParticles) < std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"The inverse number of filter particles must be greater than epsilon"};

    if (m_options.m_stepSize < std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"Step size must be positive and greater than epsilon"};

    if (m_options.m_maxUsedObservations == 0UL)
        throw std::runtime_error{"Max number of observations to use must be greater than 0"};

    // Prepare the resampling arrays
    m_pResamplingProbabilityArray = new double[m_options.m_nParticles];
    m_pResamplingParticleArray    = new unsigned int[m_options.m_nParticles];
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::~ParticleFilter()
{
    gsl_rng_free(m_pGenerator);

    if (m_pResamplingProbabilityArray)
        delete[] m_pResamplingProbabilityArray;

    if (m_pResamplingParticleArray)
        delete[] m_pResamplingParticleArray;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::DistributionHistory ParticleFilter::Filter(const Particle::History &history, const std::shared_ptr<MassHypothesis> &spMassHypothesis) const
{
    ObservedStateVector observations;

    for (const std::shared_ptr<bf::ParticleState> &spState : history)
        observations.emplace_back(spState->GetPosition(), spState->GetdEdx(), spState->Getdx());

    return this->Filter(observations, spMassHypothesis);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::DistributionHistory ParticleFilter::Filter(
    const ObservedStateVector &observations, const std::shared_ptr<MassHypothesis> &spMassHypothesis) const
{
    // Perform some checks and create the prior
    if (observations.empty())
        throw std::runtime_error{"Cannot filter on an empty observation set"};

    double currentPosition = observations.front().GetPosition() - m_options.m_stepSize;
    this->CheckOptions(spMassHypothesis);
    this->CreateDistribution(spMassHypothesis, currentPosition);
    
    const std::size_t maxObservations = std::min(observations.size(), m_options.m_maxUsedObservations);
    const double      sampleRate      = static_cast<double>(observations.size()) / static_cast<double>(maxObservations);

    // Filter on the distribution
    DistributionHistory history;
    bool                isFirstObservation = true;

    for (std::size_t i = 0UL; i < maxObservations; ++i)
    {
        const std::size_t            index       = static_cast<std::size_t>(std::round(sampleRate * static_cast<double>(i)));
        const ObservedParticleState &observation = observations.at(index);

        // Ensure observations are in order and there are no repeated positions
        const double observedPosition = observation.GetPosition();

        if (!isFirstObservation && (observedPosition - currentPosition <= std::numeric_limits<double>::epsilon()))
            throw std::runtime_error{"Cannot filter on out-of-order or position-degenerate observations"};

        // Calculate the number of propagation steps required
        std::size_t nSteps = 0UL;
    
        while (observedPosition - currentPosition > std::numeric_limits<double>::epsilon())
        {
            currentPosition += m_options.m_stepSize;
            ++nSteps;
        }

        if (!this->FilterOnObservation(observation, isFirstObservation, nSteps)) // all the particles finished or were too unlikely
        {
            // It's not terrible if there's a few extra hits at the end of a distribution - more than 10% would signify a real problem
            if (static_cast<float>(maxObservations - i) / static_cast<float>(maxObservations) > 0.1f)
                history.emplace_back(std::make_shared<ParticleDistribution>(*m_spDistribution)); // deep copy of null distribution will mean a likelihood of 0

            m_spDistribution->clear();
            return history;
        }

        history.emplace_back(std::make_shared<ParticleDistribution>(*m_spDistribution)); // deep copy
        isFirstObservation = false;
    }

    m_spDistribution->clear();
    return history;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParticleFilter::FilterOnObservation(const ObservedParticleState &observation, const bool isFirstObservation, const std::size_t nSteps) const
{
    if (!isFirstObservation)
    {
        if (this->EffectiveSampleSize() < static_cast<double>(m_options.m_nParticles) / 2.)
            this->ResampleDistribution();

        this->PropagateParticles(nSteps);
    }

    // Calculate the new particle weights
    #pragma omp parallel for
    for (auto it = m_spDistribution->begin(); it < m_spDistribution->end(); it++)
    {
        const std::shared_ptr<Particle> spParticle = it->second;

        if (spParticle->KineticEnergy() == 0.) // either the particle has truly stopped or there was an issue with its propagation
        {
            it->first = 0.;
            continue;
        }
        
        it->first = m_spPropagator->CalculateObservationProbability(spParticle, observation.GetdEdx(), observation.Getdx());
    }

    if (!this->CheckDistributionWeights()) // all the particles finished or were too unlikely; we're done
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParticleFilter::CheckDistributionWeights() const
{
    double totalWeight = 0.;

    for (const auto &pair : *m_spDistribution)
        totalWeight += pair.first;

    if (totalWeight < std::numeric_limits<double>::epsilon())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleFilter::CheckOptions(const std::shared_ptr<MassHypothesis> &spMassHypothesis) const
{
    if (spMassHypothesis->Mass() < std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"Particle mass must be positive and greater than epsilon"};

    if (spMassHypothesis->PriorEnergyDistribution().size() != m_options.m_nParticles)
        throw std::runtime_error{"The number of entries in the prior energy distirbution must match the number of filter particles"};
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleFilter::CreateDistribution(const std::shared_ptr<MassHypothesis> &spMassHypothesis, const double initialPosition) const
{
    m_spDistribution->clear();
    const double initialWeight = 1. / static_cast<double>(m_options.m_nParticles); // guaranteed to be ok by options check

    for (const double energy : spMassHypothesis->PriorEnergyDistribution())
        m_spDistribution->emplace_back(initialWeight, std::make_shared<Particle>(spMassHypothesis->Mass(), energy, initialPosition, false));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleFilter::PropagateParticles(const std::size_t nSteps) const
{
    #pragma omp parallel for
    for (auto it = m_spDistribution->begin(); it < m_spDistribution->end(); it++)
    {
        const std::shared_ptr<Particle> &spParticle = it->second;

        if (spParticle->KineticEnergy() == 0.)
            continue;

        try
        {
            for (std::size_t i = 0UL; i < nSteps && spParticle->KineticEnergy() > 0.; ++i)
                m_spPropagator->Propagate(spParticle, m_options.m_stepSize);
        }

        catch (std::runtime_error &err)
        {
            std::cerr << "Filter was forced to reject particle due to error: " << err.what() << std::endl;
            spParticle->SetKineticEnergy(0.);
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleFilter::ResampleDistribution() const
{
    // Set up probability and particle count arrays
    for (std::size_t i = 0UL; i < m_options.m_nParticles; ++i)
    {
        const auto &weightParticlePair = m_spDistribution->at(i);
        m_pResamplingProbabilityArray[i] = weightParticlePair.first;
        m_pResamplingParticleArray[i]    = 0UL;
    }

    // Sample from multinomial distribution and repopulate distribution with samples
    gsl_ran_multinomial(m_pGenerator, m_options.m_nParticles, static_cast<unsigned int>(m_options.m_nParticles), m_pResamplingProbabilityArray, m_pResamplingParticleArray);
    auto spNewDistribution = std::shared_ptr<ParticleDistribution>{new ParticleDistribution()};

    for (std::size_t i = 0UL; i < m_options.m_nParticles; ++i)
    {
        if (m_pResamplingParticleArray[i] == 0UL)
            continue;

        const std::shared_ptr<Particle> &spParticle = m_spDistribution->at(i).second;
        spNewDistribution->emplace_back(1., spParticle); // we can reuse the original particle once

        for (std::size_t j = 1UL; j < m_pResamplingParticleArray[i]; ++j)
        {
            const auto spNewParticle = std::shared_ptr<Particle>{new Particle{*spParticle}}; // a deep copy
            spNewDistribution->emplace_back(1., spNewParticle);
        }
    }

    *m_spDistribution = std::move(*spNewDistribution);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ParticleFilter::CalculateMarginalLikelihood(const DistributionHistory &distributionHistory)
{
    double marginalLogLikelihood = 0.;

    for (const std::shared_ptr<ParticleDistribution> &spDistribution : distributionHistory)
    {
        if (spDistribution->empty())
            throw std::runtime_error{"Distribution cannot be empty"};

        double avgWeight = 0.;

        for (const auto &weightParticlePair : *spDistribution)
            avgWeight += weightParticlePair.first;

        avgWeight /= static_cast<double>(spDistribution->size());

        if (avgWeight < std::numeric_limits<double>::epsilon())
        {
            marginalLogLikelihood = std::numeric_limits<double>::lowest();
            break;
        }

        marginalLogLikelihood += std::log(avgWeight);
    }

    return marginalLogLikelihood;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<double, double> ParticleFilter::CalculateFinalEnergy(const DistributionHistory &distributionHistory)
{
    // Calculate the mean
    double weightedEnergy = 0.;
    double totalWeight = 0.;

    for (const auto &weightParticlePair : *distributionHistory.back())
    {
        weightedEnergy += weightParticlePair.first * weightParticlePair.second->KineticEnergy();
        totalWeight += weightParticlePair.first;
    }

    const double mean = weightedEnergy / totalWeight;

    // Calculate the standard deviation
    double weightedSquaredEnergy = 0.;

    for (const auto &weightParticlePair : *distributionHistory.back())
        weightedSquaredEnergy += weightParticlePair.first * weightParticlePair.second->KineticEnergy() * weightParticlePair.second->KineticEnergy();

    const double denominator = (totalWeight - 1.f > std::numeric_limits<double>::epsilon()) ? totalWeight * (totalWeight - 1.f) : totalWeight * totalWeight;
    const double stdev = std::sqrt((totalWeight * weightedSquaredEnergy - weightedEnergy * weightedEnergy) / denominator);

    return {mean, stdev};
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ParticleFilter::CalculatePida(const Particle::History &history)
{
    ObservedStateVector observations;

    for (const std::shared_ptr<bf::ParticleState> &spState : history)
        observations.emplace_back(spState->GetPosition(), spState->GetdEdx(), spState->Getdx());

    return ParticleFilter::CalculatePida(observations);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ParticleFilter::CalculatePida(const ObservedStateVector &observations)
{
    double pida = 0.;

    const double maxPosition = observations.back().GetPosition();
    ObservedStateVector selectedObservations;

    for (const ObservedParticleState &observation : observations)
    {
        if (observation.GetPosition() > maxPosition - 20.f)
            selectedObservations.push_back(observation);
    }

    for (const ObservedParticleState &observation : selectedObservations)
        pida += -observation.GetdEdx() * std::pow(maxPosition - observation.GetPosition(), 0.42);
    
    return pida / static_cast<double>(selectedObservations.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ParticleFilter::EffectiveSampleSize() const
{
    // Resample the distribution if the effective sample size is too low
    double summedSquaredWeights = 0.;
    double summedWeights = 0.;

    for (const auto &weightParticlePair : *m_spDistribution)
    {
        summedSquaredWeights += weightParticlePair.first * weightParticlePair.first;
        summedWeights += weightParticlePair.first;
    }

    return summedWeights * summedWeights / summedSquaredWeights;
}

} // namespace bf
