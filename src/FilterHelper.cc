/**
 *  @file   bethe-faster/src/FilterHelper.cc
 *
 *  @brief  Implementation of the filter helper class.
 *
 *  $Log: $
 */

#include "FilterHelper.h"
#include "TRandom.h"

namespace bf
{

std::vector<double> FilterHelper::GetGaussianEnergyPrior(const std::size_t nParticles, const double mean, const double sigma) noexcept
{
    std::vector<double> prior;
    TRandom generator = TRandom{static_cast<std::uint32_t>(std::time(nullptr))};

    while (prior.size() < nParticles)
    {
        const double data = generator.Gaus(mean, sigma);

        if (data > std::numeric_limits<double>::epsilon())
            prior.push_back(data);
    }

    return prior;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<double> FilterHelper::GetUniformEnergyPrior(const std::size_t nParticles, const double min, const double max) noexcept
{
    std::vector<double> prior;
    TRandom generator = TRandom{static_cast<std::uint32_t>(std::time(nullptr))};

    while (prior.size() < nParticles)
    {
        const double data = generator.Uniform(min, max);

        if (data > std::numeric_limits<double>::epsilon())
            prior.push_back(data);
    }

    return prior;
}

} // namespace bf
