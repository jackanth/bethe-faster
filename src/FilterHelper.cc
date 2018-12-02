/**
 *  @file   bethe-faster/src/FilterHelper.cc
 *
 *  @brief  Implementation of the filter helper class.
 *
 *  $Log: $
 */

#include "FilterHelper.h"

#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"

#include <ctime>

namespace bf
{

std::vector<double> FilterHelper::GetGaussianEnergyDistribution(const std::size_t nParticles, const double mean, const double sigma)
{
    std::vector<double> distribution;
    auto generator = ROOT::Math::Random<ROOT::Math::GSLRngMT>{static_cast<std::uint32_t>(std::time(nullptr))};

    while (distribution.size() < nParticles)
    {
        const double data = generator.Gaus(mean, sigma);

        if (data > std::numeric_limits<double>::epsilon())
            distribution.push_back(data);
    }

    return distribution;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<double> FilterHelper::GetUniformEnergyDistribution(const std::size_t nParticles, const double min, const double max)
{
    std::vector<double> distribution;
    auto generator = ROOT::Math::Random<ROOT::Math::GSLRngMT>{static_cast<std::uint32_t>(std::time(nullptr))};

    while (distribution.size() < nParticles)
    {
        const double data = generator.Uniform(min, max);

        if (data > std::numeric_limits<double>::epsilon())
            distribution.push_back(data);
    }

    return distribution;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<double> FilterHelper::GetPerfectlyUniformEnergyDistribution(const std::size_t nParticles, const double min, const double max)
{
    if (min < std::numeric_limits<double>::epsilon())
        throw std::runtime_error{"Lowest energy bound must be larger than epsilon"};

    if (nParticles < 2UL)
        throw std::runtime_error{"There must be at least two particles"};

    std::vector<double> distribution;
    const double deltaE = (max - min) / static_cast<double>(nParticles - 1UL);

    for (std::size_t i = 0; i < nParticles; ++i)
        distribution.push_back(min + static_cast<double>(i) * deltaE);

    return distribution;
}

} // namespace bf
