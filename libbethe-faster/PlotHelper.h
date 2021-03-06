/**
 *  @file   bethe-faster/include/PlotHelper.h
 *
 *  @brief  Header file for the plot helper class.
 *
 *  $Log: $
 */

#ifndef BF_PLOT_HELPER_H
#define BF_PLOT_HELPER_H 1

#include "HitCharge.h"
#include "Particle.h"
#include "ParticleFilter.h"
#include "Propagator.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TStyle.h"

namespace bf
{

/**
 *  @brief  A multigraph entry
 */
class MultiGraphEntry
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  graph the graph
     *  @param  legendText the legend text
     *  @param  colour the colour
     *  @param  drawLine whether to draw a line instead of points
     *  @param  style the line width or marker style
     */
    MultiGraphEntry(TGraph graph, std::string legendText, const unsigned int colour, bool drawLine, const std::int16_t style) noexcept;

    /**
     *  @brief  Get the graph
     *
     *  @return the graph
     */
    TGraph &Graph() noexcept;

    /**
     *  @brief  Get the legend text
     *
     *  @return the legend text
     */
    const std::string &LegendText() const noexcept;

    /**
     *  @brief  Get the colour
     *
     *  @return the colour
     */
    unsigned int Colour() const noexcept;

    /**
     *  @brief  Get whether to draw a line instead of markers
     *
     *  @return whether to draw a line
     */
    bool DrawLine() const noexcept;

    /**
     *  @brief  Get the line width or marker style
     *
     *  @return the line width or marker style
     */
    std::int16_t Style() const noexcept;

private:
    TGraph       m_graph;      ///< The graph
    std::string  m_legendText; ///< The legend text
    unsigned int m_colour;     ///< Entry colour
    bool         m_drawLine;   ///< Whether to draw a line
    std::int16_t m_style;      ///< The line width or marker style
};

/**
 *  @brief  Some plotting options
 */
struct PlotOptions
{
    /**
     *  @brief  Constructor
     */
    PlotOptions() noexcept;

    std::string m_xAxisTitle; ///< The x-axis title
    std::string m_yAxisTitle; ///< The y-axis title
    bool        m_drawLegend; ///< Whether to draw the legend
    bool        m_yLogScale;  ///< Whether y has a log scale
    bool        m_xLogScale;  ///< Whether y has a log scale
    double      m_legendX1;   ///< The legend x1 coordinate
    double      m_legendX2;   ///< The legend x2 coordinate
    double      m_legendY1;   ///< The legend y1 coordinate
    double      m_legendY2;   ///< The legend y2 coordinate
};

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
     *  @param  offset the offset to use
     *
     *  @return the graph
     */
    static TGraph GetParticledEdxVersusXGraph(const std::shared_ptr<Particle> &spParticle, const bool useResidualRange = false, const double offset = 0.);

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
    static TCanvas *PlotParticleLikelihoodHistory(
        const ParticleFilter::DistributionHistory &distributionHistory, const unsigned int colour = 0UL, const std::int16_t lineWidth = 2);

    /**
     *  @brief  Get a hit charge graph
     *
     *  @param  hitChargeVector the hit charge vector
     *  @param  useResidualRange whether the use residual range instead of position
     *
     *  @return the graph
     */
    static TGraph GetHitChargeGraph(const std::vector<HitCharge> &hitChargeVector, const bool useResidualRange = true);

    /**
     *  @brief  Plot a charge graph
     *
     *  @param  hitChargeVector the hit charge vector
     *  @param  useResidualRange whether the use residual range instead of position
     *  @param  colour the marker colour
     *  @param  markerStyle the marker style
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotHitChargeGraph(const std::vector<HitCharge> &hitChargeVector, const bool useResidualRange,
        const unsigned int colour, const std::int16_t markerStyle);

    /**
     *  @brief  Get a Bragg gradient graph
     *
     *  @param  hitChargeVector the hit charge vector
     *
     *  @return the graph
     */
    static TGraph GetBraggGradientGraph(const std::vector<HitCharge> &hitChargeVector);

    /**
     *  @brief  Plot a Bragg gradient graph
     *
     *  @param  hitChargeVector the hit charge vector
     *  @param  colour the marker colour
     *  @param  markerStyle the marker style
     *
     *  @return address of the TCanvas object
     */
    static TCanvas *PlotBraggGradientGraph(const std::vector<HitCharge> &hitChargeVector, const unsigned int colour, const std::int16_t markerStyle);

    /**
     *  @brief  Get a multigraph
     *
     *  @param  graphEntries the graph entries
     *  @param  options the plotting options
     *  @param  pMultiGraph address of the new TMultiGraph
     *  @param  pLegend address of the new TLegend
     */
    static void GetMultiGraph(const std::vector<std::reference_wrapper<MultiGraphEntry>> &graphEntries, const PlotOptions &options,
        TMultiGraph *&pMultiGraph, TLegend *&pLegend);

    /**
     *  @brief  Draw a multigraph
     *
     *  @param  graphEntries the graph entries
     *  @param  options the plotting options
     *
     *  @return address of the canvas
     */
    static TCanvas *DrawMultiGraph(const std::vector<std::reference_wrapper<MultiGraphEntry>> &graphEntries, const PlotOptions &options);

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

    /**
     *  @brief  Get the light version of the indexed scheme colour
     *
     *  @param  number the number
     *
     *  @return the scheme colour
     */
    static std::int16_t GetSchemeColourLight(const unsigned int number);
};

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

inline MultiGraphEntry::MultiGraphEntry(
    TGraph graph, std::string legendText, const unsigned int colour, const bool drawLine, const std::int16_t style) noexcept :
    m_graph{std::move_if_noexcept(graph)},
    m_legendText{std::move_if_noexcept(legendText)},
    m_colour{colour},
    m_drawLine{drawLine},
    m_style{style}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline TGraph &MultiGraphEntry::Graph() noexcept
{
    return m_graph;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &MultiGraphEntry::LegendText() const noexcept
{
    return m_legendText;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int MultiGraphEntry::Colour() const noexcept
{
    return m_colour;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline bool MultiGraphEntry::DrawLine() const noexcept
{
    return m_drawLine;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline std::int16_t MultiGraphEntry::Style() const noexcept
{
    return m_style;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

inline PlotOptions::PlotOptions() noexcept :
    m_xAxisTitle{},
    m_yAxisTitle{},
    m_drawLegend{true},
    m_yLogScale{false},
    m_xLogScale{false},
    m_legendX1{0.8},
    m_legendX2{0.88},
    m_legendY1{0.72},
    m_legendY2{0.88}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

inline std::int16_t PlotHelper::GetSchemeColour(const unsigned int number)
{
    return static_cast<std::int16_t>(gStyle->GetColorPalette(number % 8UL));
}

} // namespace bf

#endif // #ifndef BF_PLOT_HELPER_H
