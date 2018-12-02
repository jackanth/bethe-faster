/**
 *  @file   bethe-faster/include/PlotHelper.h
 *
 *  @brief  Header file for the plot helper class.
 *
 *  $Log: $
 */

#ifndef BF_PLOT_HELPER_H
#define BF_PLOT_HELPER_H 1

#include "Particle.h"
#include "Propagator.h"

#include "TGraph.h"
#include "TApplication.h"

namespace bf
{

/**
 *  @brief  PlotHelper class
 */
class PlotHelper
{
public:
    /**
     *  @brief  Deleted constructor
     */
    PlotHelper() = delete;

    /**
     *  @brief  Get or initialize the TApplication object
     * 
     *  @return address of the TApplication object
     */
    static TApplication *InitializeApplication();

    /**
     *  @brief  Set the global plot style
     */
    static void SetGlobalPlotStyle();

    /**
     *  @brief  Plot a particle's dEdx
     * 
     *  @param  spParticle shared pointer to the particle
     * 
     *  @return the graph
     */
    static TGraph PlotParticledEdx(const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Plot a particle's kappa
     * 
     *  @param  propagator the propagator
     *  @param  spParticle shared pointer to the particle
     * 
     *  @return the graph
     */
    static TGraph PlotParticleKappa(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Pause to process plotting events
     */
    static void Pause();
};

} // namespace bf

#endif // #ifndef BF_PLOT_HELPER_H
