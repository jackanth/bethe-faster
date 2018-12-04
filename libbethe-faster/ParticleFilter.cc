/**
 *  @file   bethe-faster/src/ParticleFilter.cc
 *
 *  @brief  Implementation of the particle filter class.
 *
 *  $Log: $
 */

#include "ParticleFilter.h"
#include "Math/DistFunc.h"

#include <cassert>
#include <ctime>
#include <iostream>

namespace bf
{

ObservedParticleState::ObservedParticleState(const double residualRange, const double dEdx, const double dx) noexcept :
    m_residualRange{residualRange},
    m_dEdx{dEdx},
    m_dx{dx}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MassHypothesis::MassHypothesis(const double mass, std::vector<double> priorFinalEnergyDistribution) :
    m_mass{mass},
    m_priorFinalEnergyDistribution{std::move_if_noexcept(priorFinalEnergyDistribution)},
    m_numberOfParticles{m_priorFinalEnergyDistribution.size()}
{
    if (m_numberOfParticles == 0UL)
        throw std::runtime_error{"Must provide a prior final energy distribution"};
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::ParticleFilter(Detector detector) :
    m_spPropagator{new Propagator{std::move(detector)}},
    m_spDistribution{new ParticleDistribution{}},
    m_randomGenerator{static_cast<std::uint32_t>(std::time(nullptr))}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::DistributionHistory ParticleFilter::Filter(const Particle::History &history, const MassHypothesis &massHypothesis) const
{
    ObservedStateVector observations;

    for (const std::shared_ptr<bf::ParticleState> &spState : history)
        observations.emplace_back(spState->GetResidualRange(), spState->GetdEdx(), spState->Getdx());

    return this->Filter(observations, massHypothesis);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFilter::DistributionHistory ParticleFilter::Filter(ObservedStateVector observations, const MassHypothesis &massHypothesis) const
{
    // Perform some checks and create the prior
    if (observations.empty())
        throw std::runtime_error{"Cannot filter on an empty observation set"};

    this->CheckMassHypothesis(massHypothesis);
    std::sort(observations.begin(), observations.end(),
        [](const auto lhs, const auto rhs) { return lhs.GetResidualRange() < rhs.GetResidualRange(); });

    double       currentPosition = observations.front().GetResidualRange() - observations.front().Getdx();
    this->CreateDistribution(massHypothesis, currentPosition);

    // Set up the sampling vectors
    const std::size_t numParticles                  = massHypothesis.NumberOfParticles();
    const auto        spResamplingProbabilityVector = std::shared_ptr<std::vector<double>>{new std::vector<double>(numParticles, 0.)};
    const auto spResamplingParticleVector = std::shared_ptr<std::vector<unsigned int>>{new std::vector<unsigned int>(numParticles, 0UL)};

    // Filter on the distribution
    DistributionHistory history;
    bool                isFirstObservation = true;

    for (std::size_t i = 0UL, maxObservations = observations.size(); i < maxObservations; ++i)
    {
        const ObservedParticleState &observation = observations.at(i);

        // Ensure observations are in order and there are no repeated positions
        const double observedPosition = observation.GetResidualRange();

        if (!isFirstObservation && (observedPosition - currentPosition <= std::numeric_limits<double>::epsilon()))
            throw std::runtime_error{"Cannot filter on out-of-order or position-degenerate observations"};

        // Calculate the number of propagation steps required
        std::size_t nSteps = 0UL;

        while (observedPosition - currentPosition > std::numeric_limits<double>::epsilon())
        {
            currentPosition += observation.Getdx();
            ++nSteps;
        }

        if (!this->FilterOnObservation(observation, isFirstObservation, nSteps, spResamplingProbabilityVector, spResamplingParticleVector,
                numParticles)) // all the particles finished or were too unlikely
        {
            // It's not terrible if there's a few extra hits at the end of a distribution - more than 10% would signify a real problem
            if (static_cast<float>(maxObservations - i) / static_cast<float>(maxObservations) > 0.1f)
            {
                DistributionRecord distributionRecord;

                for (const auto &weightParticlePair : *m_spDistribution)
                    distributionRecord.emplace_back(weightParticlePair.first, weightParticlePair.second->KineticEnergy());

                history.emplace_back(std::make_shared<DistributionRecord>(std::move(distributionRecord)));
            }

            m_spDistribution->clear();
            return history;
        }

        DistributionRecord distributionRecord;

        for (const auto &weightParticlePair : *m_spDistribution)
            distributionRecord.emplace_back(weightParticlePair.first, weightParticlePair.second->KineticEnergy());

        history.emplace_back(std::make_shared<DistributionRecord>(std::move(distributionRecord)));
        isFirstObservation = false;
    }

    m_spDistribution->clear();
    return history;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParticleFilter::FilterOnObservation(const ObservedParticleState &observation, const bool isFirstObservation, const std::size_t nSteps,
    const std::shared_ptr<std::vector<double>> &      spResamplingProbabilityVector,
    const std::shared_ptr<std::vector<unsigned int>> &spResamplingParticleVector, const std::size_t numParticles) const
{
    if (!isFirstObservation)
    {
        if (this->EffectiveSampleSize() < static_cast<double>(numParticles) / 2.)
            this->ResampleDistribution(spResamplingProbabilityVector, spResamplingParticleVector);

        this->PropagateParticles(nSteps, observation.Getdx());
    }

// Calculate the new particle weights
#pragma omp parallel for
    for (auto it = m_spDistribution->begin(); it < m_spDistribution->end(); it++)
    {
        const std::shared_ptr<Particle> spParticle = it->second;

        if (spParticle->HasFailed()) // there was an issue with its propagation
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

void ParticleFilter::CheckMassHypothesis(const MassHypothesis &massHypothesis) const
{
    const std::size_t numParticles = massHypothesis.NumberOfParticles();

    if (numParticles < 2UL)
        throw std::runtime_error{"Number of filter particles must be 2 or greater"};

    if (1. / static_cast<double>(numParticles) < std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"The inverse number of filter particles must be greater than epsilon"};

    if (massHypothesis.Mass() < std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"Particle mass must be positive and greater than epsilon"};
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleFilter::CreateDistribution(const MassHypothesis &massHypothesis, const double initialPosition) const
{
    m_spDistribution->clear();
    const double initialWeight = 1. / static_cast<double>(massHypothesis.NumberOfParticles()); // guaranteed to be ok by options check

    for (const double energy : massHypothesis.PriorFinalEnergyDistribution())
    {
        auto spParticle = std::shared_ptr<Particle>{new Particle{massHypothesis.Mass(), energy, initialPosition, false}};
        m_spDistribution->emplace_back(initialWeight, std::move(spParticle));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleFilter::PropagateParticles(const std::size_t nSteps, const double stepSize) const
{
#pragma omp parallel for
    for (auto it = m_spDistribution->begin(); it < m_spDistribution->end(); it++)
    {
        const std::shared_ptr<Particle> &spParticle = it->second;

        try
        {
            for (std::size_t i = 0UL; i < nSteps; ++i)
                m_spPropagator->PropagateBackwards(spParticle, stepSize);
        }

        catch (std::runtime_error &err)
        {
            std::cerr << "Filter was forced to reject particle due to error: " << err.what() << std::endl;
            spParticle->HasFailed(true);
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleFilter::ResampleDistribution(const std::shared_ptr<std::vector<double>> &spResamplingProbabilityVector,
    const std::shared_ptr<std::vector<unsigned int>> &                                spResamplingParticleVector) const
{
    // Set up probability vector
    for (std::size_t i = 0UL; i < spResamplingProbabilityVector->size(); ++i)
        spResamplingProbabilityVector->at(i) = m_spDistribution->at(i).first;

    // Sample from multinomial distribution and repopulate distribution with samples
    *spResamplingParticleVector =
        m_randomGenerator.Multinomial(static_cast<unsigned int>(spResamplingProbabilityVector->size()), *spResamplingProbabilityVector);
    auto spNewDistribution = std::shared_ptr<ParticleDistribution>{new ParticleDistribution()};

    for (std::size_t i = 0UL; i < spResamplingProbabilityVector->size(); ++i)
    {
        if (spResamplingParticleVector->at(i) == 0UL)
            continue;

        const std::shared_ptr<Particle> &spParticle = m_spDistribution->at(i).second;
        spNewDistribution->emplace_back(1., spParticle); // we can reuse the original particle once

        for (std::size_t j = 1UL; j < spResamplingParticleVector->at(i); ++j)
        {
            const auto spNewParticle = std::shared_ptr<Particle>{new Particle{*spParticle}}; // a deep copy
            spNewDistribution->emplace_back(1., spNewParticle);
        }
    }

    *m_spDistribution = std::move(*spNewDistribution);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<double> ParticleFilter::GetMarginalLikelihoodHistory(const DistributionHistory &distributionHistory)
{
    if (distributionHistory.empty())
        throw std::runtime_error{"Distribution history may not be empty"};

    std::vector<double> likelihoodHistory;
    double marginalLogLikelihood = 0.;

    for (const std::shared_ptr<DistributionRecord> &spRecord : distributionHistory)
    {
        if (spRecord->empty())
            throw std::runtime_error{"Distribution may not be empty"};

        double avgWeight = 0.;

        for (const auto &weightParticlePair : *spRecord)
            avgWeight += weightParticlePair.first;

        avgWeight /= static_cast<double>(spRecord->size());

        if (avgWeight < std::numeric_limits<double>::epsilon())
        {
            likelihoodHistory.push_back(std::numeric_limits<double>::lowest());
            break;
        }

        marginalLogLikelihood += std::log(avgWeight);
        likelihoodHistory.push_back(marginalLogLikelihood);
    }

    return likelihoodHistory;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<double, double> ParticleFilter::CalculateFinalEnergy(const DistributionHistory &distributionHistory)
{
    // Calculate the mean
    double weightedEnergy = 0.;
    double totalWeight    = 0.;

    for (const auto &weightParticlePair : *distributionHistory.back())
    {
        weightedEnergy += weightParticlePair.first * weightParticlePair.second;
        totalWeight += weightParticlePair.first;
    }

    const double mean = weightedEnergy / totalWeight;

    // Calculate the standard deviation
    double weightedSquaredEnergy = 0.;

    for (const auto &weightParticlePair : *distributionHistory.back())
        weightedSquaredEnergy += weightParticlePair.first * weightParticlePair.second * weightParticlePair.second;

    const double denominator = (totalWeight - 1.f > std::numeric_limits<double>::epsilon()) ? totalWeight * (totalWeight - 1.f) : totalWeight * totalWeight;
    const double stdev = std::sqrt((totalWeight * weightedSquaredEnergy - weightedEnergy * weightedEnergy) / denominator);

    return {mean, stdev};
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ParticleFilter::CalculatePida(const Particle::History &history)
{
    ObservedStateVector observations;

    for (const std::shared_ptr<bf::ParticleState> &spState : history)
        observations.emplace_back(spState->GetResidualRange(), spState->GetdEdx(), spState->Getdx());

    return ParticleFilter::CalculatePida(observations);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ParticleFilter::CalculatePida(const ObservedStateVector &observations)
{
    double      pida          = 0.;
    std::size_t nObservations = 0UL;

    for (const ObservedParticleState &observation : observations)
    {
        if (observation.GetResidualRange() > 20.f)
            continue;

        pida += -observation.GetdEdx() * std::pow(observation.GetResidualRange(), 0.42);
        ++nObservations;
    }

    if (nObservations == 0UL)
        throw std::runtime_error{"Could not calculate PIDA because there were no observations near the end of the track"};

    return pida / static_cast<double>(nObservations);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ParticleFilter::EffectiveSampleSize() const
{
    // Resample the distribution if the effective sample size is too low
    double summedSquaredWeights = 0.;
    double summedWeights        = 0.;

    for (const auto &weightParticlePair : *m_spDistribution)
    {
        summedSquaredWeights += weightParticlePair.first * weightParticlePair.first;
        summedWeights += weightParticlePair.first;
    }

    return summedWeights * summedWeights / summedSquaredWeights;
}

} // namespace bf
