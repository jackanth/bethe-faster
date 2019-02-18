/**
 *  @file   bethe-faster/test/Test.h
 *
 *  @brief  Header file for the test script
 *
 *  $Log: $
 */

#ifndef BF_TEST_H
#define BF_TEST_H 1

#include "BetheFaster.h"

/**
 *  @brief  Test the generation step
 *
 *  @param  detector the detector
 */
void TestGeneration(const bf::Detector &detector);

/**
 *  @brief  Generate a set of particles and save the plots
 *
 *  @param  propagator the propagator
 *  @param  mode the mode
 *  @param  modeName the mode name
 *  @param  deltaX the effective thickness (cm)
 *  @param  maxEnergy the maximum energy (MeV)
 */
void GenerateParticles(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode, const std::string &modeName,
    const double deltaX, const double maxEnergy);

/**
 *  @brief  Test the filtering step
 *
 *  @param  detector the detector
 */
void TestFiltering(const bf::Detector &detector);

/**
 *  @brief  Run the filter on a particle
 *
 *  @param  detector the detector
 *  @param  propagator the propagator
 *  @param  type the particle type
 *  @param  spParticle shared pointer to the particle
 *  @param  deltaX the step size
 *  @param  maxEnergy the maximum energy
 */
void FilterOnParticle(const bf::Detector &detector, const bf::Propagator &propagator, const bf::ParticleHelper::PARTICLE_TYPE type,
    const std::shared_ptr<bf::Particle> &spParticle, const double deltaX, const double maxEnergy);

#endif // #ifndef BF_TEST_H
