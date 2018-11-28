#include "Propagator.h"
#include "ParticleHelper.h"
#include "ParticleFilter.h"

#ifdef USE_ROOT
    #include "TApplication.h"
    #include "TCanvas.h"
    #include "TGraph.h"
#endif // #ifdef USE_ROOT

#include <functional>

#ifdef USE_ROOT
void SetPlotStyle();
TApplication *Initialize();
TGraph *PlotParticledQdx(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle);
TGraph *PlotParticledEdx(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle);
TGraph *PlotParticleBetas(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle);
TGraph *PlotParticleKappas(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle);
TGraph *PlotParticleEnergies(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle);
void Pause();
TCanvas * PlotForParticles(const std::string &name, const std::function<TGraph *(const bf::Propagator &, const std::shared_ptr<bf::Particle> &)> &graphGetter,
    const std::string &yaxisLabel, const bool logY, const std::shared_ptr<bf::Particle> &spMuon, const std::shared_ptr<bf::Particle> &spProton,
    const std::shared_ptr<bf::Particle> &spPion, const std::shared_ptr<bf::Particle> &spKaon, std::int16_t *palette, const std::int16_t markerStyle, const bf::Propagator &propagator);
TCanvas * PlotForParticle(const std::string &name, const std::function<TGraph *(const bf::Propagator &, const std::shared_ptr<bf::Particle> &)> &graphGetter,
    const std::string &yaxisLabel, const bool logY, const std::shared_ptr<bf::Particle> &spParticle, std::int16_t color, const std::int16_t markerStyle, const bf::Propagator &propagator);
TCanvas * PlotDistributionError(const std::string &name, const std::function<TGraph *(const bf::Propagator &, const std::shared_ptr<bf::Particle> &)> &graphGetter,
    const std::string &yaxisLabel, const bool logY, const std::shared_ptr<bf::Particle> &spParticle, const std::int16_t color, const std::int16_t markerStyle, const bf::Propagator &propagator,
    const bf::ParticleFilter::DistributionHistory &distributionHistory, const bf::Particle::History &particleHistory, const std::int16_t errorColor);
#endif // #ifdef USE_ROOT

std::shared_ptr<bf::Particle> PropagateParticle(const bf::Propagator &propagator, const double mass, const double energy, const double deltaX);
void TestGeneration();
void TestFiltering();
