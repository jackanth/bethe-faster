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
#include "ParticleFilter.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

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
     *  @brief  Get particle energy energy
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return the graph
     */
    static TGraph GetParticleEnergyGraph(const std::shared_ptr<Particle> &spParticle, const bool useResidualRange = false);

    /**
     *  @brief  Plot a particle's energy using a line with some default settings
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  lineWidth the line width
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticleEnergyLine(const std::shared_ptr<Particle> &spParticle, const unsigned int colour = 0UL,
        const std::int16_t lineWidth = 2, const bool useResidualRange = false);

    /**
     *  @brief  Plot a particle's energy using markers with some default settings
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  markerStyle the marker style
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticleEnergyMarkers(const std::shared_ptr<Particle> &spParticle, const unsigned int colour = 0UL,
        const std::int16_t markerStyle = 6, const bool useResidualRange = false);

    /**
     *  @brief  Get particle dEdx versus x graph
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return the graph
     */
    static TGraph GetParticledEdxVersusXGraph(const std::shared_ptr<Particle> &spParticle, const bool useResidualRange = false);

    /**
     *  @brief  Plot a particle's dEdx versus x using a line with some default settings
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  lineWidth the line width
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticledEdxVersusXLine(const std::shared_ptr<Particle> &spParticle, const unsigned int colour = 0UL,
        const std::int16_t lineWidth = 2, const bool useResidualRange = false);

    /**
     *  @brief  Plot a particle's dEdx versus x using markers with some default settings
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  markerStyle the marker style
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticledEdxVersusXMarkers(const std::shared_ptr<Particle> &spParticle, const unsigned int colour = 0UL,
        const std::int16_t markerStyle = 6, const bool useResidualRange = false);

    /**
     *  @brief  Get particle dEdx versus T graph
     *
     *  @param  spParticle shared pointer to the particle
     *
     *  @return the graph
     */
    static TGraph GetParticledEdxVersusTGraph(const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Plot a particle's dEdx versus T using a line with some default settings
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  lineWidth the line width
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticledEdxVersusTLine(
        const std::shared_ptr<Particle> &spParticle, const unsigned int colour = 0UL, const std::int16_t lineWidth = 2);

    /**
     *  @brief  Plot a particle's dEdx versus T using markers with some default settings
     *
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  markerStyle the marker style
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticledEdxVersusTMarkers(
        const std::shared_ptr<Particle> &spParticle, const unsigned int colour = 0UL, const std::int16_t markerStyle = 6);

    /**
     *  @brief  Get particle kappa versus x graph
     *
     *  @param  propagator the propagator
     *  @param  spParticle shared pointer to the particle
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return the graph
     */
    static TGraph GetParticleKappaVersusXGraph(
        const Propagator &propagator, const std::shared_ptr<Particle> &spParticle, const bool useResidualRange = false);

    /**
     *  @brief  Plot a particle's kappa versus x using a line with some default settings
     *
     *  @param  propagator the propagator
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  lineWidth the line width
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticleKappaVersusXLine(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle,
        const unsigned int colour = 0UL, const std::int16_t lineWidth = 2, const bool useResidualRange = false);

    /**
     *  @brief  Plot a particle's kappa versus x using markers with some default settings
     *
     *  @param  propagator the propagator
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  markerStyle the marker style
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticleKappaVersusXMarkers(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle,
        const unsigned int colour = 0UL, const std::int16_t markerStyle = 6, const bool useResidualRange = false);

    /**
     *  @brief  Get particle kappa versus T graph
     *
     *  @param  propagator the propagator
     *  @param  spParticle shared pointer to the particle
     *
     *  @return the graph
     */
    static TGraph GetParticleKappaVersusTGraph(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle);

    /**
     *  @brief  Plot a particle's kappa versus T using a line with some default settings
     *
     *  @param  propagator the propagator
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  lineWidth the line width
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticleKappaVersusTLine(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle,
        const unsigned int colour = 0UL, const std::int16_t lineWidth = 2);

    /**
     *  @brief  Plot a particle's kappa versus T using markers with some default settings
     *
     *  @param  propagator the propagator
     *  @param  spParticle shared pointer to the particle
     *  @param  colour the line colour
     *  @param  markerStyle the marker style
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticleKappaVersusTMarkers(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle,
        const unsigned int colour = 0UL, const std::int16_t markerStyle = 6);

    /**
     *  @brief  Get a particle likelihood history graph
     *
     *  @param  distributionHistory the distribution history
     *
     *  @return the graph
     */
    static TGraph GetParticleLikelihoodHistoryGraph(const ParticleFilter::DistributionHistory &distributionHistory);

    /**
     *  @brief  Plot a particle's likelihood history with some default settings
     *
     *  @param  distributionHistory the distribution history
     *  @param  colour the line colour
     *  @param  lineWidth the line width
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotParticleLikelihoodHistory(const ParticleFilter::DistributionHistory &distributionHistory,
        const unsigned int colour = 0UL, const std::int16_t lineWidth = 2);

    /**
     *  @brief  Pause to process plotting events
     */
    static void Pause();

    /**
     *  @brief  Get the indexed scheme colour
     *
     *  @param  number the number
     *
     *  @return the scheme colour
     */
    static std::int16_t GetSchemeColour(const unsigned int number);
};

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

inline std::int16_t PlotHelper::GetSchemeColour(const unsigned int number)
{
    return static_cast<std::int16_t>(gStyle->GetColorPalette(number % 8UL));
}

} // namespace bf

#endif // #ifndef BF_PLOT_HELPER_H
