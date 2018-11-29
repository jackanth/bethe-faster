/**
 *  @file   bethe-faster/include/FilterHelper.h
 *
 *  @brief  Header file for the filter helper class.
 *
 *  $Log: $
 */
#ifndef BF_FILTER_HELPER_H
#define BF_FILTER_HELPER_H 1

#include "Detector.h"
#include <vector>

namespace bf
{

/**
 *  @brief  FilterHelper class
 */
class FilterHelper
{
public:
    /**
     *  @brief  Deleted constructor
     */
    FilterHelper() = delete;

    /**
     *  @brief  Get a Gaussian energy prior in MeV
     * 
     *  @param nParticles the number of particles
     *  @param mean the mean energy (MeV)
     *  @param sigma the standard deviation of the energy (MeV)
     * 
     *  @return the energy prior
     */
    static std::vector<double> GetGaussianEnergyPrior(const std::size_t nParticles, const double mean, const double sigma);

    /**
     *  @brief  Get a uniform energy prior in MeV
     * 
     *  @param nParticles the number of particles
     *  @param min the minimum energy (MeV)
     *  @param max the maximum energy (MeV)
     * 
     *  @return the energy prior
     */
    static std::vector<double> GetUniformEnergyPrior(const std::size_t nParticles, const double min, const double max);

    /**
     *  @brief  Get an evenly-distributed uniform energy prior in MeV
     * 
     *  @param nParticles the number of particles
     *  @param min the minimum energy (MeV)
     *  @param max the maximum energy (MeV)
     * 
     *  @return the energy prior
     */
    static std::vector<double> GetPerfectlyUniformEnergyPrior(const std::size_t nParticles, const double min, const double max);
};

} // namespace bf

#endif // #ifndef BF_FILTER_HELPER_H
