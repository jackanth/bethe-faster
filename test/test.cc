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

    PlotOptions options;

    // Plot the energy graph
    auto muonEnergyGraph        = bf::PlotHelper::GetParticleEnergyGraph(spMuon);
    auto chargedPionEnergyGraph = bf::PlotHelper::GetParticleEnergyGraph(spChargedPion);
    auto chargedKaonEnergyGraph = bf::PlotHelper::GetParticleEnergyGraph(spChargedKaon);
    auto protonEnergyGraph      = bf::PlotHelper::GetParticleEnergyGraph(spProton);

    options.m_isLine     = (mode != bf::Propagator::PROPAGATION_MODE::STOCHASTIC);
    options.m_style      = options.m_isLine ? 2 : 6;
    options.m_fileName   = "energy" + nameSuffix;
    options.m_xAxisTitle = "x\\text{ (cm)}";
    options.m_yAxisTitle = "T \\text{ (MeV)}";
    options.m_yLogScale  = false;

    WriteMultiGraph(muonEnergyGraph, chargedPionEnergyGraph, chargedKaonEnergyGraph, protonEnergyGraph, options);

    // Plot the dE/dx versus x graph
    auto muondEdxVersusXGraph        = bf::PlotHelper::GetParticledEdxVersusXGraph(spMuon);
    auto chargedPiondEdxVersusXGraph = bf::PlotHelper::GetParticledEdxVersusXGraph(spChargedPion);
    auto chargedKaondEdxVersusXGraph = bf::PlotHelper::GetParticledEdxVersusXGraph(spChargedKaon);
    auto protondEdxVersusXGraph      = bf::PlotHelper::GetParticledEdxVersusXGraph(spProton);

    options.m_isLine     = (mode != bf::Propagator::PROPAGATION_MODE::STOCHASTIC);
    options.m_style      = options.m_isLine ? 2 : 6;
    options.m_fileName   = "dEdx_versus_x" + nameSuffix;
    options.m_xAxisTitle = "x\\text{ (cm)}";
    options.m_yAxisTitle = "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}";
    options.m_yLogScale  = true;

    if (!options.m_isLine && (deltaX < 0.1))
        options.m_style = 1;

    WriteMultiGraph(muondEdxVersusXGraph, chargedPiondEdxVersusXGraph, chargedKaondEdxVersusXGraph, protondEdxVersusXGraph, options);

    // Plot the dE/dx versus T graph
    auto muondEdxVersusTGraph        = bf::PlotHelper::GetParticledEdxVersusTGraph(spMuon);
    auto chargedPiondEdxVersusTGraph = bf::PlotHelper::GetParticledEdxVersusTGraph(spChargedPion);
    auto chargedKaondEdxVersusTGraph = bf::PlotHelper::GetParticledEdxVersusTGraph(spChargedKaon);
    auto protondEdxVersusTGraph      = bf::PlotHelper::GetParticledEdxVersusTGraph(spProton);

    options.m_isLine     = (mode != bf::Propagator::PROPAGATION_MODE::STOCHASTIC);
    options.m_style      = options.m_isLine ? 2 : 6;
    options.m_fileName   = "dEdx_versus_T" + nameSuffix;
    options.m_xAxisTitle = "T\\text{ (MeV)}";
    options.m_yAxisTitle = "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}";
    options.m_yLogScale  = true;

    if (!options.m_isLine && (deltaX < 0.1))
        options.m_style = 1;

    WriteMultiGraph(muondEdxVersusTGraph, chargedPiondEdxVersusTGraph, chargedKaondEdxVersusTGraph, protondEdxVersusTGraph, options);

    // Plot the kappa versus x graph
    auto muonKappaVersusXGraph        = bf::PlotHelper::GetParticleKappaVersusXGraph(propagator, spMuon);
    auto chargedPionKappaVersusXGraph = bf::PlotHelper::GetParticleKappaVersusXGraph(propagator, spChargedPion);
    auto chargedKaonKappaVersusXGraph = bf::PlotHelper::GetParticleKappaVersusXGraph(propagator, spChargedKaon);
    auto protonKappaVersusXGraph      = bf::PlotHelper::GetParticleKappaVersusXGraph(propagator, spProton);

    options.m_isLine     = (mode != bf::Propagator::PROPAGATION_MODE::STOCHASTIC);
    options.m_style      = options.m_isLine ? 2 : 6;
    options.m_fileName   = "kappa_versus_x" + nameSuffix;
    options.m_xAxisTitle = "x\\text{ (cm)}";
    options.m_yAxisTitle = "\\kappa";
    options.m_yLogScale  = true;

    WriteMultiGraph(muonKappaVersusXGraph, chargedPionKappaVersusXGraph, chargedKaonKappaVersusXGraph, protonKappaVersusXGraph, options);

    // Plot the kappa versus T graph
    auto muonKappaVersusTGraph        = bf::PlotHelper::GetParticleKappaVersusTGraph(propagator, spMuon);
    auto chargedPionKappaVersusTGraph = bf::PlotHelper::GetParticleKappaVersusTGraph(propagator, spChargedPion);
    auto chargedKaonKappaVersusTGraph = bf::PlotHelper::GetParticleKappaVersusTGraph(propagator, spChargedKaon);
    auto protonKappaVersusTGraph      = bf::PlotHelper::GetParticleKappaVersusTGraph(propagator, spProton);

    options.m_isLine     = (mode != bf::Propagator::PROPAGATION_MODE::STOCHASTIC);
    options.m_style      = options.m_isLine ? 2 : 6;
    options.m_fileName   = "kappa_versus_T" + nameSuffix;
    options.m_xAxisTitle = "T\\text{ (MeV)}";
    options.m_yAxisTitle = "\\kappa";
    options.m_yLogScale  = true;

    WriteMultiGraph(muonKappaVersusTGraph, chargedPionKappaVersusTGraph, chargedKaonKappaVersusTGraph, protonKappaVersusTGraph, options);
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
    auto muonGraph        = bf::PlotHelper::GetParticleLikelihoodHistoryGraph(muonDistributionHistory);
    auto chargedPionGraph = bf::PlotHelper::GetParticleLikelihoodHistoryGraph(chargedPionDistributionHistory);
    auto chargedKaonGraph = bf::PlotHelper::GetParticleLikelihoodHistoryGraph(chargedKaonDistributionHistory);
    auto protonGraph      = bf::PlotHelper::GetParticleLikelihoodHistoryGraph(protonDistributionHistory);

    PlotOptions options;
    options.m_isLine = true;
    options.m_style  = 2;
    options.m_fileName = "likelihood_" + bf::ParticleHelper::ToString(type) + "_" + std::to_string(deltaX) + "cm_" + std::to_string(maxEnergy) + "MeV.eps";
    options.m_xAxisTitle = "\\text{Observation number}";
    options.m_yAxisTitle = "\\text{Log-likelihood}";
    options.m_yLogScale  = false;

    WriteMultiGraph(muonGraph, chargedPionGraph, chargedKaonGraph, protonGraph, options);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void WriteMultiGraph(TGraph &muonGraph, TGraph &chargedPionGraph, TGraph &chargedKaonGraph, TGraph &protonGraph, const PlotOptions &options)
{
    // Create the canvas and get the graphs
    TCanvas canvas{"bfCanvas", "bfCanvas", 800, 600};
    canvas.SetGrid();

    if (options.m_yLogScale)
        canvas.SetLogy();

    // Format the graphs
    const auto muonColour = bf::PlotHelper::GetSchemeColour(0UL);
    muonGraph.SetFillColor(muonColour);

    const auto chargedPionColour = bf::PlotHelper::GetSchemeColour(1UL);
    chargedPionGraph.SetFillColor(chargedPionColour);

    const auto chargedKaonColour = bf::PlotHelper::GetSchemeColour(2UL);
    chargedKaonGraph.SetFillColor(chargedKaonColour);

    const auto protonColour = bf::PlotHelper::GetSchemeColour(3UL);
    protonGraph.SetFillColor(protonColour);

    if (options.m_isLine)
    {
        muonGraph.SetLineWidth(options.m_style);
        chargedPionGraph.SetLineWidth(options.m_style);
        chargedKaonGraph.SetLineWidth(options.m_style);
        protonGraph.SetLineWidth(options.m_style);

        muonGraph.SetLineColor(muonColour);
        chargedPionGraph.SetLineColor(chargedPionColour);
        chargedKaonGraph.SetLineColor(chargedKaonColour);
        protonGraph.SetLineColor(protonColour);
    }

    else
    {
        muonGraph.SetMarkerStyle(options.m_style);
        muonGraph.SetMarkerStyle(options.m_style);
        chargedPionGraph.SetMarkerStyle(options.m_style);
        chargedKaonGraph.SetMarkerStyle(options.m_style);
        protonGraph.SetMarkerStyle(options.m_style);

        muonGraph.SetMarkerColor(muonColour);
        chargedPionGraph.SetMarkerColor(chargedPionColour);
        chargedKaonGraph.SetMarkerColor(chargedKaonColour);
        protonGraph.SetMarkerColor(protonColour);
    }

    // Create and draw the multi-graph
    TMultiGraph multiGraph{};
    multiGraph.Add(&muonGraph);
    multiGraph.Add(&chargedPionGraph);
    multiGraph.Add(&chargedKaonGraph);
    multiGraph.Add(&protonGraph);

    multiGraph.GetXaxis()->SetTitle(options.m_xAxisTitle.c_str());
    multiGraph.GetYaxis()->SetTitle(options.m_yAxisTitle.c_str());

    multiGraph.Draw(options.m_isLine ? "AL" : "AP");

    // Add legend
    TLegend legend{0.8, 0.72, 0.88, 0.88};
    legend.AddEntry(&muonGraph, "\\mu^-", "f");
    legend.AddEntry(&chargedPionGraph, "\\pi^\\pm", "f");
    legend.AddEntry(&chargedKaonGraph, "K^\\pm", "f");
    legend.AddEntry(&protonGraph, "p\\", "f");
    legend.SetBorderSize(1);
    legend.Draw();

    canvas.SaveAs(options.m_fileName.c_str());
}