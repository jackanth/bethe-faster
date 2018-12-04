/**
 *  @file   bethe-faster/test/Test.cc
 *
 *  @brief  Implementation of the tests.
 *
 *  $Log: $
 */

#include "Test.h"

#include "TAxis.h"
#include "TLegend.h"
#include "TMultiGraph.h"

#include <iostream>

int main()
{
    try
    {
        bf::PlotHelper::SetGlobalPlotStyle();
        const auto detector = bf::DetectorHelper::GetMicroBooNEDetector();

        std::cout << "Testing generation step" << std::endl;
        TestGeneration(detector);

        std::cout << "Testing filtering step" << std::endl;
        TestFiltering(detector);
    }

    catch (const std::runtime_error &err)
    {
        std::cerr << "Testing failed: " << err.what() << std::endl;
        return 1;
    }

    catch (...)
    {
        std::cerr << "Testing failed but exception was not caught" << std::endl;
        return 1;
    }

    return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TestGeneration(const bf::Detector &detector)
{
    const auto propagator = bf::Propagator{detector};

    // Generate at effective thickness 0.03cm
    GenerateParticles(propagator, bf::Propagator::PROPAGATION_MODE::STOCHASTIC, "stochastic", 0.03, 1000.f);
    GenerateParticles(propagator, bf::Propagator::PROPAGATION_MODE::MEAN, "mean", 0.03, 1000.f);
    GenerateParticles(propagator, bf::Propagator::PROPAGATION_MODE::MODAL, "modal", 0.03, 1000.f);

    // Generate at effective thickness 0.3cm
    GenerateParticles(propagator, bf::Propagator::PROPAGATION_MODE::STOCHASTIC, "stochastic", 0.3, 1000.f);
    GenerateParticles(propagator, bf::Propagator::PROPAGATION_MODE::MEAN, "mean", 0.3, 1000.f);
    GenerateParticles(propagator, bf::Propagator::PROPAGATION_MODE::MODAL, "modal", 0.3, 1000.f);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void GenerateParticles(const bf::Propagator &propagator, const bf::Propagator::PROPAGATION_MODE mode, const std::string &modeName,
    const double deltaX, const double maxEnergy)
{
    // Define the particles
    auto spMuon        = bf::ParticleHelper::GetMuon();
    auto spProton      = bf::ParticleHelper::GetProton();
    auto spChargedPion = bf::ParticleHelper::GetChargedPion();
    auto spChargedKaon = bf::ParticleHelper::GetChargedKaon();

    // Propagate particles
    while ((spMuon->KineticEnergy() < maxEnergy) && !spMuon->HasFailed())
        propagator.PropagateBackwards(spMuon, deltaX, mode);

    while ((spProton->KineticEnergy() < maxEnergy) && !spProton->HasFailed())
        propagator.PropagateBackwards(spProton, deltaX, mode);

    while ((spChargedPion->KineticEnergy() < maxEnergy) && !spChargedPion->HasFailed())
        propagator.PropagateBackwards(spChargedPion, deltaX, mode);

    while ((spChargedKaon->KineticEnergy() < maxEnergy) && !spChargedKaon->HasFailed())
        propagator.PropagateBackwards(spChargedKaon, deltaX, mode);

    // Draw plots
    const std::string nameSuffix = "_" + modeName + "_" + std::to_string(deltaX) + "cm_" + std::to_string(maxEnergy) + "MeV.eps";
    bf::PlotOptions options;
    const bool drawLine = (mode != bf::Propagator::PROPAGATION_MODE::STOCHASTIC);
    const std::int16_t style = drawLine ? 2 : 6;
    const std::int16_t altStyle = (!drawLine && (deltaX < 0.1)) ? 1 : style;

    // Plot the energy graph
    auto muonEnergyGraph        = bf::MultiGraphEntry{bf::PlotHelper::GetParticleEnergyGraph(spMuon), "\\mu", 0UL, drawLine, style};
    auto chargedPionEnergyGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleEnergyGraph(spChargedPion), "\\pi^\\pm", 1UL, drawLine, style};
    auto chargedKaonEnergyGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleEnergyGraph(spChargedKaon), "K^\\pm", 2UL, drawLine, style};
    auto protonEnergyGraph      = bf::MultiGraphEntry{bf::PlotHelper::GetParticleEnergyGraph(spProton), "p\\", 3UL, drawLine, style};

    options.m_xAxisTitle = "x\\text{ (cm)}";
    options.m_yAxisTitle = "T \\text{ (MeV)}";
    options.m_yLogScale  = false;
    options.m_xLogScale  = false;
    options.m_drawLegend = true;

    bf::PlotHelper::DrawMultiGraph({muonEnergyGraph, chargedPionEnergyGraph, chargedKaonEnergyGraph, protonEnergyGraph}, options)->SaveAs(("energy" + nameSuffix).c_str());

    // Plot the dE/dx versus x graph
    auto muondEdxVersusXGraph        = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spMuon), "\\mu", 0UL, drawLine, altStyle};
    auto chargedPiondEdxVersusXGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spChargedPion), "\\pi^\\pm", 1UL, drawLine, altStyle};
    auto chargedKaondEdxVersusXGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spChargedKaon), "K^\\pm", 2UL, drawLine, altStyle};
    auto protondEdxVersusXGraph      = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusXGraph(spProton), "p\\", 3UL, drawLine, altStyle};

    options.m_xAxisTitle = "x\\text{ (cm)}";
    options.m_yAxisTitle = "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}";
    options.m_yLogScale  = true;
    options.m_xLogScale  = false;
    options.m_drawLegend = true;

    bf::PlotHelper::DrawMultiGraph({muondEdxVersusXGraph, chargedPiondEdxVersusXGraph, chargedKaondEdxVersusXGraph, protondEdxVersusXGraph}, options)->SaveAs(("dEdx_versus_x" + nameSuffix).c_str());

    // Plot the dE/dx versus T graph
    auto muondEdxVersusTGraph        = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusTGraph(spMuon), "\\mu", 0UL, drawLine, altStyle};
    auto chargedPiondEdxVersusTGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusTGraph(spChargedPion), "\\pi^\\pm", 1UL, drawLine, altStyle};
    auto chargedKaondEdxVersusTGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusTGraph(spChargedKaon), "K^\\pm", 2UL, drawLine, altStyle};
    auto protondEdxVersusTGraph      = bf::MultiGraphEntry{bf::PlotHelper::GetParticledEdxVersusTGraph(spProton), "p\\", 3UL, drawLine, altStyle};

    options.m_xAxisTitle = "T\\text{ (MeV)}";
    options.m_yAxisTitle = "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}";
    options.m_yLogScale  = true;
    options.m_xLogScale  = false;
    options.m_drawLegend = true;

    bf::PlotHelper::DrawMultiGraph({muondEdxVersusTGraph, chargedPiondEdxVersusTGraph, chargedKaondEdxVersusTGraph, protondEdxVersusTGraph}, options)->SaveAs(("dEdx_versus_T" + nameSuffix).c_str());

    // Plot the kappa versus x graph
    auto muonKappaVersusXGraph        = bf::MultiGraphEntry{bf::PlotHelper::GetParticleKappaVersusXGraph(propagator, spMuon), "\\mu", 0UL, drawLine, style};
    auto chargedPionKappaVersusXGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleKappaVersusXGraph(propagator, spChargedPion), "\\pi^\\pm", 1UL, drawLine, style};
    auto chargedKaonKappaVersusXGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleKappaVersusXGraph(propagator, spChargedKaon), "K^\\pm", 2UL, drawLine, style};
    auto protonKappaVersusXGraph      = bf::MultiGraphEntry{bf::PlotHelper::GetParticleKappaVersusXGraph(propagator, spProton), "p\\", 3UL, drawLine, style};

    options.m_xAxisTitle = "x\\text{ (cm)}";
    options.m_yAxisTitle = "\\kappa";
    options.m_yLogScale  = true;
    options.m_xLogScale  = false;
    options.m_drawLegend = true;

    bf::PlotHelper::DrawMultiGraph({muonKappaVersusXGraph, chargedPionKappaVersusXGraph, chargedKaonKappaVersusXGraph, protonKappaVersusXGraph}, options)->SaveAs(("kappa_versus_x" + nameSuffix).c_str());

    // Plot the kappa versus T graph
    auto muonKappaVersusTGraph        = bf::MultiGraphEntry{bf::PlotHelper::GetParticleKappaVersusTGraph(propagator, spMuon), "\\mu", 0UL, drawLine, style};
    auto chargedPionKappaVersusTGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleKappaVersusTGraph(propagator, spChargedPion), "\\pi^\\pm", 1UL, drawLine, style};
    auto chargedKaonKappaVersusTGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleKappaVersusTGraph(propagator, spChargedKaon), "K^\\pm", 2UL, drawLine, style};
    auto protonKappaVersusTGraph      = bf::MultiGraphEntry{bf::PlotHelper::GetParticleKappaVersusTGraph(propagator, spProton), "p\\", 3UL, drawLine, style};

    options.m_xAxisTitle = "T\\text{ (MeV)}";
    options.m_yAxisTitle = "\\kappa";
    options.m_yLogScale  = true;
    options.m_xLogScale  = false;
    options.m_drawLegend = true;

    bf::PlotHelper::DrawMultiGraph({muonKappaVersusTGraph, chargedPionKappaVersusTGraph, chargedKaonKappaVersusTGraph, protonKappaVersusTGraph}, options)->SaveAs(("kappa_versus_T" + nameSuffix).c_str());
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TestFiltering(const bf::Detector &detector)
{
    // Define the particles
    auto spMuon        = bf::ParticleHelper::GetMuon();
    auto spChargedPion = bf::ParticleHelper::GetChargedPion();
    auto spChargedKaon = bf::ParticleHelper::GetChargedKaon();
    auto spProton      = bf::ParticleHelper::GetProton();

    const auto propagator = bf::Propagator{detector};
    FilterOnParticle(detector, propagator, bf::ParticleHelper::PARTICLE_TYPE::MUON, spMuon, 0.3, 500.);
    FilterOnParticle(detector, propagator, bf::ParticleHelper::PARTICLE_TYPE::CHARGED_PION, spChargedPion, 0.3, 500.);
    FilterOnParticle(detector, propagator, bf::ParticleHelper::PARTICLE_TYPE::CHARGED_KAON, spChargedKaon, 0.3, 500.);
    FilterOnParticle(detector, propagator, bf::ParticleHelper::PARTICLE_TYPE::PROTON, spProton, 0.3, 500.);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void FilterOnParticle(const bf::Detector &detector, const bf::Propagator &propagator, const bf::ParticleHelper::PARTICLE_TYPE type,
    const std::shared_ptr<bf::Particle> &spParticle, const double deltaX, const double maxEnergy)
{
    // Prepare the test particle
    spParticle->Reset();

    while (spParticle->KineticEnergy() < maxEnergy)
        propagator.PropagateBackwards(spParticle, deltaX);

    // Prepare the hypotheses
    const auto finalEnergyDistribution = bf::FilterHelper::GetGaussianEnergyDistribution(1000UL, 5., 10.);
    const auto muonHypothesis          = bf::MassHypothesis{bf::PhysicalConstants::m_muonMass, finalEnergyDistribution};
    const auto chargedPionHypothesis   = bf::MassHypothesis{bf::PhysicalConstants::m_chargedPionMass, finalEnergyDistribution};
    const auto chargedKaonHypothesis   = bf::MassHypothesis{bf::PhysicalConstants::m_chargedKaonMass, finalEnergyDistribution};
    const auto protonHypothesis        = bf::MassHypothesis{bf::PhysicalConstants::m_protonMass, finalEnergyDistribution};

    // Run the filter
    const auto particleFilter = bf::ParticleFilter{detector};

    const auto   muonDistributionHistory = particleFilter.Filter(spParticle->GetHistory(), muonHypothesis);
    const double muonLogLikelihood       = bf::ParticleFilter::CalculateMarginalLikelihood(muonDistributionHistory);

    const auto   chargedPionDistributionHistory = particleFilter.Filter(spParticle->GetHistory(), chargedPionHypothesis);
    const double chargedPionLogLikelihood       = bf::ParticleFilter::CalculateMarginalLikelihood(chargedPionDistributionHistory);

    const auto   chargedKaonDistributionHistory = particleFilter.Filter(spParticle->GetHistory(), chargedKaonHypothesis);
    const double chargedKaonLogLikelihood       = bf::ParticleFilter::CalculateMarginalLikelihood(chargedKaonDistributionHistory);

    const auto   protonDistributionHistory = particleFilter.Filter(spParticle->GetHistory(), protonHypothesis);
    const double protonLogLikelihood       = bf::ParticleFilter::CalculateMarginalLikelihood(protonDistributionHistory);

    // Check the results
    auto predictedType = bf::ParticleHelper::PARTICLE_TYPE::OTHER;

    if ((muonLogLikelihood > chargedPionLogLikelihood) && (muonLogLikelihood > chargedKaonLogLikelihood) && (muonLogLikelihood > protonLogLikelihood))
        predictedType = bf::ParticleHelper::PARTICLE_TYPE::MUON;

    else if ((chargedPionLogLikelihood > muonLogLikelihood) && (chargedPionLogLikelihood > chargedKaonLogLikelihood) &&
             (chargedPionLogLikelihood > protonLogLikelihood))
        predictedType = bf::ParticleHelper::PARTICLE_TYPE::CHARGED_PION;

    else if ((chargedKaonLogLikelihood > chargedPionLogLikelihood) && (chargedKaonLogLikelihood > muonLogLikelihood) &&
             (chargedKaonLogLikelihood > protonLogLikelihood))
        predictedType = bf::ParticleHelper::PARTICLE_TYPE::CHARGED_KAON;

    else
        predictedType = bf::ParticleHelper::PARTICLE_TYPE::PROTON;

    // State the results
    if (type == predictedType)
        std::cout << "Correctly predicted " << bf::ParticleHelper::ToString(type) << std::endl;

    else
        std::cout << "Mispredicted " << bf::ParticleHelper::ToString(type) << " as " << bf::ParticleHelper::ToString(predictedType) << std::endl;

    // Plot the likelihood graphs
    auto muonGraph        = bf::MultiGraphEntry{bf::PlotHelper::GetParticleLikelihoodHistoryGraph(muonDistributionHistory), "\\mu", 0UL, true, 2};
    auto chargedPionGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleLikelihoodHistoryGraph(chargedPionDistributionHistory), "\\pi^\\pm", 1UL, true, 2};
    auto chargedKaonGraph = bf::MultiGraphEntry{bf::PlotHelper::GetParticleLikelihoodHistoryGraph(chargedKaonDistributionHistory), "K^\\pm", 2UL, true, 2};
    auto protonGraph      = bf::MultiGraphEntry{bf::PlotHelper::GetParticleLikelihoodHistoryGraph(protonDistributionHistory), "p\\", 3UL, true, 2};

    bf::PlotOptions options;
    options.m_xAxisTitle = "\\text{Observation number}";
    options.m_yAxisTitle = "\\text{Log-likelihood}";
    options.m_yLogScale  = false;
    options.m_xLogScale  = false;
    options.m_drawLegend = true;

    const auto fileName = "likelihood_" + bf::ParticleHelper::ToString(type) + "_" + std::to_string(deltaX) + "cm_" + std::to_string(maxEnergy) + "MeV.eps";
    bf::PlotHelper::DrawMultiGraph({muonGraph, chargedPionGraph, chargedKaonGraph, protonGraph}, options)->SaveAs(fileName.c_str());
}