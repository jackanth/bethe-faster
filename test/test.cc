#include "test.h"
#include "DetectorHelper.h"
#include "FilterHelper.h"

#ifdef USE_ROOT
    #include "TAxis.h"
    #include "TColor.h"
    #include "TLegend.h"
    #include "TMultiGraph.h"
    #include "TROOT.h"
    #include "TStyle.h"
    #include "TSystem.h"
    #include "TVector.h"
    #include "TLine.h"
    #include "TText.h"
    #include "TH1F.h"
    #include "TGaxis.h"
    #include "TF1.h"
#endif // #ifdef USE_ROOT

#include <iostream>
#include <chrono>

#ifdef USE_ROOT
void SetPlotStyle()
{
    TStyle *pStyle = new TStyle("NewStyle", "NewStyle");
    pStyle->SetOptStat(0);
    pStyle->SetOptTitle(0);
    pStyle->SetOptDate(0);
    pStyle->SetLabelFont(42, "xyz");
    pStyle->SetLabelSize(0.03f, "xyz");
    pStyle->SetTitleSize(0.035f, "xyz");
    pStyle->SetTitleFont(42, "xyz");
    pStyle->SetCanvasDefW(500);
    pStyle->SetCanvasDefH(500);
    pStyle->SetCanvasColor(0);
    pStyle->SetPadColor(0);
    pStyle->SetCanvasBorderMode(0);
    pStyle->SetCanvasBorderSize(0);
    pStyle->SetPadBorderMode(0);
    pStyle->SetPadBorderSize(0);
    pStyle->SetPadBottomMargin(0.1f);
    pStyle->SetPadTopMargin(0.1f);
    pStyle->SetPadRightMargin(0.1f);
    pStyle->SetPadGridX(0);
    pStyle->SetPadGridY(0);
    pStyle->SetPadTickX(1);
    pStyle->SetPadTickY(1);
    pStyle->SetHistFillColor(0);
    pStyle->SetFrameBorderMode(0);
    gROOT->SetStyle("NewStyle");
}
#endif // #ifdef USE_ROOT

#ifdef USE_ROOT
TApplication *Initialize()
{
    TApplication *pApplication = nullptr;

    if (gApplication)
        pApplication = gApplication;

    else
    {
        int   argc   = 0;
        char *argv[1UL];
        argv[0] = nullptr;
        pApplication = new TApplication("PandoraMonitoring", &argc, const_cast<char **>(argv));
        pApplication->SetReturnFromRun(kTRUE);
    }

    SetPlotStyle();
    return pApplication;
}
#endif // #ifdef USE_ROOT

#ifdef USE_ROOT
TGraph *PlotParticledEdx(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    (void) propagator;

    std::vector<double> positions, responses;
    const auto &       history = spParticle->GetHistory();

    for (const auto &spState : history)
    {
        positions.push_back(spState->GetPosition());
        responses.push_back(-spState->GetdEdx());
    }

    return new TGraph(static_cast<Int_t>(history.size()), positions.data(), responses.data());
}
#endif // #ifdef USE_ROOT

#ifdef USE_ROOT
TGraph *PlotParticleBetas(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    (void) propagator;

    std::vector<double> positions, betas;
    const auto &       history = spParticle->GetHistory();

    for (const auto &spState : history)
    {
        positions.push_back(spState->GetPosition());
        betas.push_back(bf::ParticleHelper::GetParticleBeta(spParticle->Mass(), spState->GetKineticEnergy()));
    }

    return new TGraph(static_cast<Int_t>(history.size()), positions.data(), betas.data());
}
#endif // #ifdef USE_ROOT

#ifdef USE_ROOT
TGraph *PlotParticleKappas(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<double> positions, kappas;
    const auto &       history = spParticle->GetHistory();

    for (const auto &spState : history)
    {
        positions.push_back(spState->GetPosition());
        kappas.push_back(propagator.CalculateKappa(spParticle->Mass(), spState->GetKineticEnergy(), spState->Getdx()));
    }

    return new TGraph(static_cast<Int_t>(history.size()), positions.data(), kappas.data());
}
#endif // #ifdef USE_ROOT

#ifdef USE_ROOT
TGraph *PlotParticleEnergies(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    (void) propagator;
    
    std::vector<double> positions, responses;
    const auto &       history = spParticle->GetHistory();

    for (const auto &spState : history)
    {
        positions.push_back(spState->GetPosition());
        responses.push_back(spState->GetKineticEnergy());
    }

    return new TGraph(static_cast<Int_t>(history.size()), positions.data(), responses.data());
}
#endif // #ifdef USE_ROOT

std::shared_ptr<bf::Particle> PropagateParticle(const bf::Propagator &propagator, const double mass, const double energy, const double deltaX)
{
    const auto  spParticle = std::make_shared<bf::Particle>(mass, energy);

    while (spParticle->IsAlive())
        propagator.Propagate(spParticle, deltaX);

    return spParticle;
}

#ifdef USE_ROOT
void Pause()
{
    std::cout << "Press return to continue ..." << std::endl;
    int flag = fcntl(1, F_GETFL, 0);

    int key = 0;
    while (true)
    {
        gSystem->ProcessEvents();
        (void)fcntl(1, F_SETFL, flag | O_NONBLOCK);
        key = getchar();

        if ((key == '\n') || (key == '\r'))
            break;

        usleep(1000);
    }

    fcntl(1, F_SETFL, flag);
}
#endif // #ifdef USE_ROOT

#ifdef USE_ROOT
TCanvas * PlotForParticles(const std::string &name, const std::function<TGraph *(const bf::Propagator &, const std::shared_ptr<bf::Particle> &)> &graphGetter,
    const std::string &yaxisLabel, const bool logY, const std::shared_ptr<bf::Particle> &spMuon, const std::shared_ptr<bf::Particle> &spProton,
    const std::shared_ptr<bf::Particle> &spPion, const std::shared_ptr<bf::Particle> &spKaon, std::int16_t *palette, const std::int16_t markerStyle, const bf::Propagator &propagator)
{
    TCanvas *pCanvas = new TCanvas(name.c_str(), name.c_str(), 800, 600);

    if (logY)
        pCanvas->SetLogy();

    TGraph *pMuonResponse = graphGetter(propagator, spMuon);
    pMuonResponse->SetMarkerStyle(markerStyle);
    pMuonResponse->SetMarkerColor(palette[0]);
    pMuonResponse->SetFillColor(palette[0]);

    TGraph *pProtonResponse = graphGetter(propagator, spProton);
    pProtonResponse->SetMarkerStyle(markerStyle);
    pProtonResponse->SetMarkerColor(palette[1]);
    pProtonResponse->SetFillColor(palette[1]);

    TGraph *pPionResponse = graphGetter(propagator, spPion);
    pPionResponse->SetMarkerStyle(markerStyle);
    pPionResponse->SetMarkerColor(palette[2]);
    pPionResponse->SetFillColor(palette[2]);

    TGraph *pKaonResponse = graphGetter(propagator, spKaon);
    pKaonResponse->SetMarkerStyle(markerStyle);
    pKaonResponse->SetMarkerColor(palette[3]);
    pKaonResponse->SetFillColor(palette[3]);

    TMultiGraph *pMultiGraph = new TMultiGraph();
    pMultiGraph->Add(pMuonResponse);
    pMultiGraph->Add(pProtonResponse);
    pMultiGraph->Add(pPionResponse);
    pMultiGraph->Add(pKaonResponse);
    pMultiGraph->Draw("ap");

    // Extra formatting.
    pCanvas->SetGrid();
    pMultiGraph->GetYaxis()->SetTitle(yaxisLabel.c_str());
    pMultiGraph->GetXaxis()->SetTitle("x\\text{ (cm)}");
    pMultiGraph->GetXaxis()->SetLimits(0., 250.);

    // Redraw.
    pMultiGraph->Draw("ap");
    pCanvas->Update();

    // Add legend.
    TLegend *pLegend = new TLegend(0.8, 0.72, 0.88, 0.88);
    pLegend->AddEntry(pMuonResponse, "\\mu^-", "f");
    pLegend->AddEntry(pProtonResponse, "p\\", "f");
    pLegend->AddEntry(pPionResponse, "\\pi^\\pm", "f");
    pLegend->AddEntry(pKaonResponse, "K^\\pm", "f");
    pLegend->SetBorderSize(1);
    pLegend->Draw();

    return pCanvas;
}
#endif // #ifdef USE_ROOT

#ifdef USE_ROOT
TCanvas * PlotForParticle(const std::string &name, const std::function<TGraph *(const bf::Propagator &, const std::shared_ptr<bf::Particle> &)> &graphGetter,
    const std::string &yaxisLabel, const bool logY, const std::shared_ptr<bf::Particle> &spParticle, const std::int16_t color, const std::int16_t markerStyle, const bf::Propagator &propagator)
{
     TCanvas *pCanvas = new TCanvas(name.c_str(), name.c_str(), 800, 600);

    if (logY)
        pCanvas->SetLogy();

    TGraph *pParticleResponse = graphGetter(propagator, spParticle);
    pParticleResponse->SetMarkerStyle(markerStyle);
    pParticleResponse->SetMarkerColor(color);

    pParticleResponse->Draw("ap");

    // Extra formatting.
    pCanvas->SetGrid();
    pParticleResponse->GetYaxis()->SetTitle(yaxisLabel.c_str());
    pParticleResponse->GetXaxis()->SetTitle("x\\text{ (cm)}");
    pParticleResponse->GetXaxis()->SetLimits(0, 250.);

    // Redraw.
    pParticleResponse->Draw("ap");
    pCanvas->Update();

    return pCanvas;
}
#endif // #ifdef USE_ROOT

#ifdef USE_ROOT
TCanvas * PlotDistributionError(const std::string &name, const std::function<TGraph *(const bf::Propagator &, const std::shared_ptr<bf::Particle> &)> &graphGetter,
    const std::string &yaxisLabel, const bool logY, const std::shared_ptr<bf::Particle> &spParticle, const std::int16_t color, const std::int16_t markerStyle, const bf::Propagator &propagator,
    const bf::ParticleFilter::DistributionHistory &distributionHistory, const bf::Particle::History &particleHistory, const std::int16_t errorColor)
{
    TCanvas *pCanvas = new TCanvas(name.c_str(), name.c_str(), 800, 600);

    // Get the response graph
    TGraph *pParticleResponse = graphGetter(propagator, spParticle);
    pParticleResponse->SetMarkerStyle(markerStyle);
    pParticleResponse->SetMarkerColor(color);
    pParticleResponse->SetFillColor(color);

    // Get the error graph
    std::vector<double> positions, errors;
   // const std::size_t nExamples = std::min(distributionHistory.size(), particleHistory.size());
   // const double      sampleRate      = static_cast<double>(particleHistory.size()) / static_cast<double>(distributionHistory.size());

    for (std::size_t i = 0UL; i < distributionHistory.size(); ++i)
    {
        const auto &spDistribution = distributionHistory.at(i);
        const auto &spState        = particleHistory.at(i);

        double meanEnergy = 0.f;
        double totalWeight = 0.f;

        for (const auto &pair : *spDistribution)
        {
            meanEnergy += pair.first * pair.second->KineticEnergy();
            totalWeight += pair.first;
        }

        meanEnergy /= totalWeight;

        const double error = meanEnergy - spState->GetKineticEnergy();
        positions.push_back(spState->GetPosition());
        errors.push_back(error);
    }

    TGraph *pErrorGraph = new TGraph(static_cast<Int_t>(positions.size()), positions.data(), errors.data());
    pErrorGraph->SetMarkerColor(errorColor);
    pErrorGraph->SetFillColor(errorColor);
    pErrorGraph->SetMarkerStyle(markerStyle);

    TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
    p1->SetFillStyle(0); // transparent

    TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
    p2->SetFillStyle(0); // transparent

    p1->SetTicky(0);
     if (logY)
        {}//p1->SetLogy();

    p2->SetTicky(0);
    p1->Draw();
    p1->cd();

    pParticleResponse->GetYaxis()->SetTitle(yaxisLabel.c_str());
    pParticleResponse->GetXaxis()->SetTitle("x\\text{ (cm)}");
    pParticleResponse->GetXaxis()->SetLimits(0., 250.);
    pParticleResponse->Draw("AP");
    
    // Add legend.
    TLegend *pLegend = new TLegend(0.68, 0.72, 0.88, 0.88);
    pLegend->AddEntry(pParticleResponse, "-dQ/dx\\", "f");
    pLegend->AddEntry(pErrorGraph, "\\text{Mean KE estimate error}", "f");
    pLegend->SetBorderSize(1);
    pLegend->Draw();

    gPad->Update();

    Double_t xmin = p1->GetUxmin();
    Double_t xmax = p1->GetUxmax();
    Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
    Double_t ymin = pErrorGraph->GetHistogram()->GetMinimum();
    Double_t ymax = pErrorGraph->GetHistogram()->GetMaximum();
    Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
    p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
    p2->Draw();
    p2->cd();
    pErrorGraph->Draw("P");
    gPad->Update();

    TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
    axis->SetLabelSize(0.03f);
    axis->SetTitleSize(0.035f);
    axis->SetLabelFont(42);
    axis->SetTitleFont(42);
    axis->SetTitle("Mean KE estimate error (MeV)");

    axis->Draw();
    gPad->Update();
    pCanvas->cd();

    return pCanvas;
}
#endif // #ifdef USE_ROOT

void TestGeneration()
{
    const double deltaX = 0.03;

    // Propagate some particles
    const auto propagator = bf::Propagator(bf::DetectorHelper::GetMicroBooNEDetector());

    const auto spMuon   = bf::ParticleHelper::GetMuon(500.);
    const auto spProton = bf::ParticleHelper::GetProton(500.);
    const auto spPion   = bf::ParticleHelper::GetChargedPion(500.);
    const auto spKaon   = bf::ParticleHelper::GetChargedKaon(500.);

    propagator.PropagateUntilStopped(deltaX, {spMuon, spProton, spPion, spKaon});

#ifdef USE_ROOT
    const std::int16_t markerStyle = 1; // 1, 6 are 7 are dots of increasing size

    // Plot the results
    std::int16_t palette[8];
    palette[0] = static_cast<std::int16_t>(TColor::GetColor("#1B9E77"));
    palette[1] = static_cast<std::int16_t>(TColor::GetColor("#D95F02"));
    palette[2] = static_cast<std::int16_t>(TColor::GetColor("#7570B3"));
    palette[3] = static_cast<std::int16_t>(TColor::GetColor("#E7298A"));
    palette[4] = static_cast<std::int16_t>(TColor::GetColor("#66A61E"));
    palette[5] = static_cast<std::int16_t>(TColor::GetColor("#E6AB02"));
    palette[6] = static_cast<std::int16_t>(TColor::GetColor("#A6761D"));
    palette[7] = static_cast<std::int16_t>(TColor::GetColor("#666666"));

    PlotForParticle("dEdx_muon", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spMuon, palette[0], markerStyle, propagator);
    PlotForParticle("dEdx_proton", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spProton, palette[1], markerStyle, propagator);
    PlotForParticle("dEdx_pion", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spPion, palette[2], markerStyle, propagator);
    PlotForParticle("dEdx_kaon", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spKaon, palette[3], markerStyle, propagator);

    PlotForParticles("beta", PlotParticleBetas, "\\beta", false, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);
    PlotForParticles("energy", PlotParticleEnergies, "T \\text{ (MeV)}", false, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);
    PlotForParticles("dEdx", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);

    TCanvas * pKappaCanvas = PlotForParticles("kappa", PlotParticleKappas, "\\kappa", true, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);

    TText *pTextLandau = new TText(11., 4.e-4, "Truncated Landau");
    pTextLandau->SetTextFont(42);
    pTextLandau->SetTextSize(0.04f);
    pTextLandau->Draw();

    TLine *pLineLandau = new TLine(0., 0.01, pKappaCanvas->GetUxmax(), 0.01);
    pLineLandau->SetLineColor(kBlack);
    pLineLandau->SetLineWidth(2);
    pLineLandau->SetLineStyle(2);
    pLineLandau->Draw();

    TText *pTextConvolution = new TText(11., 0.5, "Vavilov");
    pTextConvolution->SetTextFont(42);
    pTextConvolution->SetTextSize(0.04f);
    pTextConvolution->Draw();

    TLine *pLineLogNormal = new TLine(0., 10., pKappaCanvas->GetUxmax(), 10.);
    pLineLogNormal->SetLineColor(kBlack);
    pLineLogNormal->SetLineWidth(2);
    pLineLogNormal->SetLineStyle(2);
    pLineLogNormal->Draw();

    TText *pTextNormal = new TText(11., 20., "Normal");
    pTextNormal->SetTextFont(42);
    pTextNormal->SetTextSize(0.04f);
    pTextNormal->Draw();
    
    Pause();
#endif // #ifdef USE_ROOT
}

void TestFiltering()
{
    std::size_t pidaCorrect = 0UL;
    std::size_t meCorrect = 0UL;
    std::size_t total = 0UL;

    const auto spPropagator = std::make_shared<bf::Propagator>(bf::DetectorHelper::GetMicroBooNEDetector());

    // Set up the particle filter
    bf::FilterOptions filterOptions;
    filterOptions.m_stepSize            = 0.03;
    filterOptions.m_nParticles          = 1000UL;
    filterOptions.m_maxUsedObservations = 10000UL;
    const auto particleFilter = bf::ParticleFilter{spPropagator, filterOptions};

    for (int i = 0; i < 100; ++i)
    {
        std::cout << i << std::endl;
        // Propagate a particle
        const double particleMass   = bf::PhysicalConstants::m_protonMass;
        const double particleEnergy = 500.; // MeV
        const double deltaX         = 0.03; // cm

        const auto spParticle = bf::ParticleHelper::GetParticle(particleMass, particleEnergy);
        spPropagator->PropagateUntilStopped(deltaX, spParticle);

        const auto &particleHistory     = spParticle->GetHistory();
        //const std::size_t nObservations = particleHistory.size();

        //Some plotting options
        const std::int16_t markerStyle = 6; // 1, 6 are 7 are dots of increasing size
        
        std::int16_t palette[8];
        palette[0] = static_cast<std::int16_t>(TColor::GetColor("#1B9E77"));
        palette[1] = static_cast<std::int16_t>(TColor::GetColor("#D95F02"));
        palette[2] = static_cast<std::int16_t>(TColor::GetColor("#7570B3"));
        palette[3] = static_cast<std::int16_t>(TColor::GetColor("#E7298A"));
        palette[4] = static_cast<std::int16_t>(TColor::GetColor("#66A61E"));
        palette[5] = static_cast<std::int16_t>(TColor::GetColor("#E6AB02"));
        palette[6] = static_cast<std::int16_t>(TColor::GetColor("#A6761D"));
        palette[7] = static_cast<std::int16_t>(TColor::GetColor("#666666"));

        // Prepare the hypotheses
        auto muonPrior = bf::FilterHelper::GetPerfectlyUniformEnergyPrior(filterOptions.m_nParticles, 50., 10000.);
        const auto spMuonHypothesis = std::shared_ptr<bf::MassHypothesis>{new bf::MassHypothesis{bf::PhysicalConstants::m_muonMass, std::move(muonPrior)}};

        auto pionPrior = bf::FilterHelper::GetPerfectlyUniformEnergyPrior(filterOptions.m_nParticles, 50., 10000.);
        const auto spPionHypothesis = std::shared_ptr<bf::MassHypothesis>{new bf::MassHypothesis{bf::PhysicalConstants::m_chargedPionMass, std::move(pionPrior)}};

        auto kaonPrior = bf::FilterHelper::GetPerfectlyUniformEnergyPrior(filterOptions.m_nParticles, 50., 10000.);
        const auto spKaonHypothesis = std::shared_ptr<bf::MassHypothesis>{new bf::MassHypothesis{bf::PhysicalConstants::m_chargedKaonMass, std::move(kaonPrior)}};

        auto protonPrior = bf::FilterHelper::GetPerfectlyUniformEnergyPrior(filterOptions.m_nParticles, 50., 1000.);
        const auto spProtonHypothesis = std::shared_ptr<bf::MassHypothesis>{new bf::MassHypothesis{bf::PhysicalConstants::m_protonMass, std::move(protonPrior)}};

        // Test the muon hypothesis
        const auto muonDistributionHistory = particleFilter.Filter(particleHistory, spMuonHypothesis);
        const double muonLogLikelihood = bf::ParticleFilter::CalculateMarginalLikelihood(muonDistributionHistory);

        PlotDistributionError("distribution_error_muon", PlotParticledEdx, "-\\mathrm{d}Q/\\mathrm{d}x \\text{ (MeV/cm)}", true, spParticle, palette[0], markerStyle, *spPropagator, muonDistributionHistory, particleHistory, palette[5]);
                const auto protonDistributionHistory = particleFilter.Filter(particleHistory, spProtonHypothesis);


        auto distributionHistory = protonDistributionHistory;
        const std::size_t numViews = std::min(20UL, distributionHistory.size());
        const std::size_t observationFrequency = static_cast<std::size_t>(distributionHistory.size() / numViews);
        TCanvas *pCanvas = new TCanvas("filter_progress", "filter_progress", 800, 600);

        for (std::size_t j = 0UL; j < numViews; ++j)
        {
            TH1F *pHistogram = new TH1F("particle_filter", "particle_filter", 50, 0., 1000.);

            const std::size_t index    = j * observationFrequency;
            const auto &spDistribution = distributionHistory.at(index);
            const auto &spState        = particleHistory.at(index);

            double totalWeight = 0.;
            for (const auto &pair : *spDistribution)
                totalWeight += pair.first;

            for (const auto &pair : *spDistribution)
                pHistogram->Fill(pair.second->KineticEnergy(), pair.first / totalWeight);

            pHistogram->GetXaxis()->SetTitle("Kinetic energy (MeV)");
            pHistogram->GetYaxis()->SetTitle("Probability density");

            pHistogram->Draw("hist");
            gPad->Update();

            // Add the truth line
            TLine *pTrueLine = new TLine(spState->GetKineticEnergy(), 0., spState->GetKineticEnergy(), pCanvas->GetUymax());
            pTrueLine->SetLineColor(kBlack);
            pTrueLine->SetLineWidth(2);
            pTrueLine->SetLineStyle(2);
            pTrueLine->Draw();

            // Add a label
            TText *pLabel = new TText(30., pCanvas->GetUymax() * 0.9, ("Observation " + std::to_string(index)).c_str());
            pLabel->SetTextFont(42);
            pLabel->SetTextSize(0.04f);
            pLabel->Draw();

            gPad->Update();
            Pause();

            delete pHistogram;
        }

        double muonMeanFinalEnergy = 0.f, muonStdevFinalEnergy = 0.f;
        std::tie(muonMeanFinalEnergy, muonStdevFinalEnergy) = bf::ParticleFilter::CalculateFinalEnergy(muonDistributionHistory);

        

        // Test the pion hypothesis
        const auto pionDistributionHistory = particleFilter.Filter(particleHistory, spPionHypothesis);
        const double pionLogLikelihood = bf::ParticleFilter::CalculateMarginalLikelihood(pionDistributionHistory);

        double pionMeanFinalEnergy = 0.f, pionStdevFinalEnergy = 0.f;
        std::tie(pionMeanFinalEnergy, pionStdevFinalEnergy) = bf::ParticleFilter::CalculateFinalEnergy(pionDistributionHistory);

        // Test the kaon hypothesis
        const auto kaonDistributionHistory = particleFilter.Filter(particleHistory, spKaonHypothesis);
        const double kaonLogLikelihood = bf::ParticleFilter::CalculateMarginalLikelihood(kaonDistributionHistory);

        double kaonMeanFinalEnergy = 0.f, kaonStdevFinalEnergy = 0.f;
        std::tie(kaonMeanFinalEnergy, kaonStdevFinalEnergy) = bf::ParticleFilter::CalculateFinalEnergy(kaonDistributionHistory);

        std::cout << "Kaon log-likelihood was " << kaonLogLikelihood << ", final T = " << kaonMeanFinalEnergy << " +- " << kaonStdevFinalEnergy << " MeV"<< std::endl;
        PlotDistributionError("distribution_error_kaon", PlotParticledEdx, "-\\mathrm{d}Q/\\mathrm{d}x \\text{ (MeV/cm)}", true, spParticle, palette[0], markerStyle, *spPropagator, kaonDistributionHistory, particleHistory, palette[5]);
        Pause();

        // Test the proton hypothesis
        const double protonLogLikelihood = bf::ParticleFilter::CalculateMarginalLikelihood(protonDistributionHistory);

        double protonMeanFinalEnergy = 0.f, protonStdevFinalEnergy = 0.f;
        std::tie(protonMeanFinalEnergy, protonStdevFinalEnergy) = bf::ParticleFilter::CalculateFinalEnergy(protonDistributionHistory);

        std::cout << "Proton log-likelihood was " << protonLogLikelihood << ", final T = " << protonMeanFinalEnergy << " +- " << protonStdevFinalEnergy << " MeV"<< std::endl;
        PlotDistributionError("distribution_error_proton", PlotParticledEdx, "-\\mathrm{d}Q/\\mathrm{d}x \\text{ (MeV/cm)}", true, spParticle, palette[0], markerStyle, *spPropagator, protonDistributionHistory, particleHistory, palette[5]);
        Pause();

        // if (muonLogLikelihood > protonLogLikelihood && muonLogLikelihood > pionLogLikelihood && muonLogLikelihood > kaonLogLikelihood)
        //     std::cout << "This was a muon" << std::endl;

        // else if (protonLogLikelihood > muonLogLikelihood && protonLogLikelihood > pionLogLikelihood && protonLogLikelihood > kaonLogLikelihood)
        //     std::cout << "This was a proton" << std::endl;

        if (muonLogLikelihood > pionLogLikelihood && muonLogLikelihood > protonLogLikelihood && muonLogLikelihood > kaonLogLikelihood)
        {
            ++meCorrect;
            std::cout << "Me correct" << std::endl;
        }

        else
        {
            std::cout << "Me incorrect" << std::endl;
        }
        // else
        //     std::cout << "This was a kaon" << std::endl;

        const double pida = bf::ParticleFilter::CalculatePida(particleHistory);

        const double pidaPionError = std::abs(pida - 8.52747);
        const double pidaMuonError = std::abs(pida - 7.82418);
        const double pidaProtonError = std::abs(pida - 17.4296);
        const double pidaKaonError = std::abs(pida - 13.4559);

        if (pidaMuonError < pidaPionError && pidaMuonError < pidaProtonError && pidaMuonError < pidaKaonError)
        {
            ++pidaCorrect;
            std::cout << "PIDA correct" << std::endl;
        }

        else
            std::cout << "PIDA incorrect" << std::endl;

        
        std::cout << "Pion log-likelihood was " << pionLogLikelihood << ", final T = " << pionMeanFinalEnergy << " +- " << pionStdevFinalEnergy << " MeV"<< std::endl;
        PlotDistributionError("distribution_error_pion", PlotParticledEdx, "-\\mathrm{d}Q/\\mathrm{d}x \\text{ (MeV/cm)}", true, spParticle, palette[0], markerStyle, *spPropagator, pionDistributionHistory, particleHistory, palette[5]);
        Pause();

        // else if (pidaMuonError < pidaPionError && pidaMuonError < pidaProtonError && pidaMuonError < pidaKaonError)
        //     std::cout << "PIDA says this was a muon" << std::endl;

        // else if (pidaProtonError < pidaMuonError && pidaProtonError < pidaPionError && pidaProtonError < pidaKaonError)
        //     std::cout << "PIDA says this was a proton" << std::endl;

        // else
        //     std::cout << "PIDA says this was a kaon" << std::endl;

        ++total;
    }

    std::cout << "Pions pida correct = " << pidaCorrect << std::endl;
    std::cout << "Pions me correct   = " << meCorrect << std::endl;
    std::cout << "Pions total        = " << total << std::endl;

    pidaCorrect = 0UL;
    meCorrect = 0UL;
    total = 0UL;

    for (int i = 0; i < 100; ++i)
    {
        std::cout << i << std::endl;
        // Propagate a particle
        const double particleMass   = bf::PhysicalConstants::m_chargedPionMass;
        const double particleEnergy = 500.; // MeV
        const double deltaX         = 0.3; // cm
        
        const auto spParticle = bf::ParticleHelper::GetParticle(particleMass, particleEnergy);
        spPropagator->PropagateUntilStopped(deltaX, spParticle);

        const auto &particleHistory     = spParticle->GetHistory();
        //const std::size_t nObservations = particleHistory.size();

        // Some plotting options
        //const std::int16_t markerStyle = 6; // 1, 6 are 7 are dots of increasing size
        
        // std::int16_t palette[8];
        // palette[0] = static_cast<std::int16_t>(TColor::GetColor("#1B9E77"));
        // palette[1] = static_cast<std::int16_t>(TColor::GetColor("#D95F02"));
        // palette[2] = static_cast<std::int16_t>(TColor::GetColor("#7570B3"));
        // palette[3] = static_cast<std::int16_t>(TColor::GetColor("#E7298A"));
        // palette[4] = static_cast<std::int16_t>(TColor::GetColor("#66A61E"));
        // palette[5] = static_cast<std::int16_t>(TColor::GetColor("#E6AB02"));
        // palette[6] = static_cast<std::int16_t>(TColor::GetColor("#A6761D"));
        // palette[7] = static_cast<std::int16_t>(TColor::GetColor("#666666"));

        // Prepare the hypotheses
        auto muonPrior = bf::FilterHelper::GetUniformEnergyPrior(filterOptions.m_nParticles, 50., 10000.);
        const auto spMuonHypothesis = std::shared_ptr<bf::MassHypothesis>{new bf::MassHypothesis{bf::PhysicalConstants::m_muonMass, std::move(muonPrior)}};

        auto pionPrior = bf::FilterHelper::GetUniformEnergyPrior(filterOptions.m_nParticles, 50., 10000.);
        const auto spPionHypothesis = std::shared_ptr<bf::MassHypothesis>{new bf::MassHypothesis{bf::PhysicalConstants::m_chargedPionMass, std::move(pionPrior)}};

        auto kaonPrior = bf::FilterHelper::GetUniformEnergyPrior(filterOptions.m_nParticles, 50., 10000.);
        const auto spKaonHypothesis = std::shared_ptr<bf::MassHypothesis>{new bf::MassHypothesis{bf::PhysicalConstants::m_chargedKaonMass, std::move(kaonPrior)}};

        auto protonPrior = bf::FilterHelper::GetUniformEnergyPrior(filterOptions.m_nParticles, 50., 1000.);
        const auto spProtonHypothesis = std::shared_ptr<bf::MassHypothesis>{new bf::MassHypothesis{bf::PhysicalConstants::m_protonMass, std::move(protonPrior)}};

        // Test the muon hypothesis
        const auto muonDistributionHistory = particleFilter.Filter(particleHistory, spMuonHypothesis);
        const double muonLogLikelihood = bf::ParticleFilter::CalculateMarginalLikelihood(muonDistributionHistory);

        double muonMeanFinalEnergy = 0.f, muonStdevFinalEnergy = 0.f;
        std::tie(muonMeanFinalEnergy, muonStdevFinalEnergy) = bf::ParticleFilter::CalculateFinalEnergy(muonDistributionHistory);

        //std::cout << "Muon log-likelihood was " << muonLogLikelihood << ", final T = " << muonMeanFinalEnergy << " +- " << muonStdevFinalEnergy << " MeV"<< std::endl;
        // PlotDistributionError("distribution_error_muon", PlotParticledEdx, "-\\mathrm{d}Q/\\mathrm{d}x \\text{ (MeV/cm)}", false, spParticle, palette[0], markerStyle, propagator, muonDistributionHistory, particleHistory, palette[5]);
        // Pause();

        // Test the pion hypothesis
        const auto pionDistributionHistory = particleFilter.Filter(particleHistory, spPionHypothesis);
        const double pionLogLikelihood = bf::ParticleFilter::CalculateMarginalLikelihood(pionDistributionHistory);

        double pionMeanFinalEnergy = 0.f, pionStdevFinalEnergy = 0.f;
        std::tie(pionMeanFinalEnergy, pionStdevFinalEnergy) = bf::ParticleFilter::CalculateFinalEnergy(pionDistributionHistory);

        //std::cout << "Pion log-likelihood was " << pionLogLikelihood << ", final T = " << pionMeanFinalEnergy << " +- " << pionStdevFinalEnergy << " MeV"<< std::endl;
        // PlotDistributionError("distribution_error_pion", PlotParticledEdx, "-\\mathrm{d}Q/\\mathrm{d}x \\text{ (MeV/cm)}", false, spParticle, palette[0], markerStyle, propagator, pionDistributionHistory, particleHistory, palette[5]);
        // Pause();

        // Test the kaon hypothesis
        const auto kaonDistributionHistory = particleFilter.Filter(particleHistory, spKaonHypothesis);
        const double kaonLogLikelihood = bf::ParticleFilter::CalculateMarginalLikelihood(kaonDistributionHistory);

        double kaonMeanFinalEnergy = 0.f, kaonStdevFinalEnergy = 0.f;
        std::tie(kaonMeanFinalEnergy, kaonStdevFinalEnergy) = bf::ParticleFilter::CalculateFinalEnergy(kaonDistributionHistory);

       // std::cout << "Kaon log-likelihood was " << kaonLogLikelihood << ", final T = " << kaonMeanFinalEnergy << " +- " << kaonStdevFinalEnergy << " MeV"<< std::endl;
        // PlotDistributionError("distribution_error_kaon", PlotParticledEdx, "-\\mathrm{d}Q/\\mathrm{d}x \\text{ (MeV/cm)}", false, spParticle, palette[0], markerStyle, propagator, kaonDistributionHistory, particleHistory, palette[5]);
        // Pause();

        // Test the proton hypothesis
        const auto protonDistributionHistory = particleFilter.Filter(particleHistory, spProtonHypothesis);
        const double protonLogLikelihood = bf::ParticleFilter::CalculateMarginalLikelihood(protonDistributionHistory);

        double protonMeanFinalEnergy = 0.f, protonStdevFinalEnergy = 0.f;
        std::tie(protonMeanFinalEnergy, protonStdevFinalEnergy) = bf::ParticleFilter::CalculateFinalEnergy(protonDistributionHistory);

       // std::cout << "Proton log-likelihood was " << protonLogLikelihood << ", final T = " << protonMeanFinalEnergy << " +- " << protonStdevFinalEnergy << " MeV"<< std::endl;
        // PlotDistributionError("distribution_error_proton", PlotParticledEdx, "-\\mathrm{d}Q/\\mathrm{d}x \\text{ (MeV/cm)}", false, spParticle, palette[0], markerStyle, propagator, protonDistributionHistory, particleHistory, palette[5]);
        // Pause();

        // if (muonLogLikelihood > protonLogLikelihood && muonLogLikelihood > pionLogLikelihood && muonLogLikelihood > kaonLogLikelihood)
        //     std::cout << "This was a muon" << std::endl;

        // else if (protonLogLikelihood > muonLogLikelihood && protonLogLikelihood > pionLogLikelihood && protonLogLikelihood > kaonLogLikelihood)
        //     std::cout << "This was a proton" << std::endl;

        if (pionLogLikelihood > muonLogLikelihood && pionLogLikelihood > protonLogLikelihood && pionLogLikelihood > kaonLogLikelihood)
            ++meCorrect;

        // else
        //     std::cout << "This was a kaon" << std::endl;

        const double pida = bf::ParticleFilter::CalculatePida(particleHistory);

        const double pidaPionError = std::abs(pida - 8.52747);
        const double pidaMuonError = std::abs(pida - 7.82418);
        const double pidaProtonError = std::abs(pida - 17.4296);
        const double pidaKaonError = std::abs(pida - 13.4559);

        if (pidaPionError < pidaMuonError && pidaPionError < pidaProtonError && pidaPionError < pidaKaonError)
            ++pidaCorrect;

        // else if (pidaMuonError < pidaPionError && pidaMuonError < pidaProtonError && pidaMuonError < pidaKaonError)
        //     std::cout << "PIDA says this was a muon" << std::endl;

        // else if (pidaProtonError < pidaMuonError && pidaProtonError < pidaPionError && pidaProtonError < pidaKaonError)
        //     std::cout << "PIDA says this was a proton" << std::endl;

        // else
        //     std::cout << "PIDA says this was a kaon" << std::endl;

        ++total;
    }

    std::cout << "Muons pida correct = " << pidaCorrect << std::endl;
    std::cout << "Muons me correct   = " << meCorrect << std::endl;
    std::cout << "Muons total        = " << total << std::endl;
}

int main()
{  
    Initialize();
    const bool testGeneration = true;
    const bool testFiltering = false;

    try
    {
#ifdef USE_ROOT
        Initialize();
#endif // #ifdef USE_ROOT

        if (testGeneration)
        {
            std::cout << "Testing generation step" << std::endl;
            TestGeneration();
        }

        if (testFiltering)
        {
            std::cout << "Testing filtering step" << std::endl;
            TestFiltering();
        }
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
