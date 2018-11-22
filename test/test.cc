#include "DetectorHelper.h"
#include "ParticleHelper.h"
#include "Propagator.h"
#include "ParticleFilter.h"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVector.h"
#include "TLine.h"
#include "TText.h"
#include "TH1F.h"

#include <iostream>

void SetPlotStyle()
{
    TStyle *pStyle = new TStyle("NewStyle", "NewStyle");
    pStyle->SetOptStat(0);
    pStyle->SetOptTitle(0);
    pStyle->SetOptDate(0);
    pStyle->SetLabelSize(0.03, "xyz");
    pStyle->SetTitleSize(0.035, "xyz");
    pStyle->SetCanvasDefW(500);
    pStyle->SetCanvasDefH(500);
    pStyle->SetCanvasColor(0);
    pStyle->SetCanvasBorderMode(0);
    pStyle->SetCanvasBorderSize(0);
    pStyle->SetPadBottomMargin(0.1);
    pStyle->SetPadTopMargin(0.1);
    pStyle->SetPadRightMargin(0.1);
    pStyle->SetPadGridX(0);
    pStyle->SetPadGridY(0);
    pStyle->SetPadTickX(1);
    pStyle->SetPadTickY(1);
    pStyle->SetFrameBorderMode(0);
    gROOT->SetStyle("NewStyle");
}

TApplication *Initialize()
{
    TApplication *pApplication = nullptr;

    if (gApplication)
        pApplication = gApplication;

    else
    {
        int   argc   = 0;
        char *argv   = (char *)"";
        pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
        pApplication->SetReturnFromRun(kTRUE);
    }

    SetPlotStyle();
    return pApplication;
}

TGraph *PlotParticledQdx(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<float> positions, responses;
    const auto &       history = spParticle->GetHistory();

    for (const auto &state : history)
    {
        positions.push_back(state.GetPosition());
        responses.push_back(propagator.CalculateResponse(state.GetdEdx()));
    }

    return new TGraph(history.size(), positions.data(), responses.data());
}

TGraph *PlotParticledEdx(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<float> positions, responses;
    const auto &       history = spParticle->GetHistory();

    for (const auto &state : history)
    {
        positions.push_back(state.GetPosition());
        responses.push_back(-state.GetdEdx());
    }

    return new TGraph(history.size(), positions.data(), responses.data());
}

TGraph *PlotParticleBetas(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<float> positions, betas;
    const auto &       history = spParticle->GetHistory();

    for (const auto &state : history)
    {
        positions.push_back(state.GetPosition());
        betas.push_back(bf::ParticleHelper::GetParticleBeta(spParticle->Mass(), state.GetKineticEnergy()));
    }

    return new TGraph(history.size(), positions.data(), betas.data());
}

TGraph *PlotParticleKappas(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<float> positions, kappas;
    const auto &       history = spParticle->GetHistory();

    for (const auto &state : history)
    {
        positions.push_back(state.GetPosition());
        kappas.push_back(propagator.CalculateKappa(spParticle->Mass(), state.GetKineticEnergy()));
    }

    return new TGraph(history.size(), positions.data(), kappas.data());
}

TGraph *PlotParticleEnergies(const bf::Propagator &propagator, const std::shared_ptr<bf::Particle> &spParticle)
{
    std::vector<float> positions, responses;
    const auto &       history = spParticle->GetHistory();

    for (const auto &state : history)
    {
        positions.push_back(state.GetPosition());
        responses.push_back(state.GetKineticEnergy());
    }

    return new TGraph(history.size(), positions.data(), responses.data());
}

std::shared_ptr<bf::Particle> PropagateParticle(const bf::Propagator &propagator, const float mass, const float energy, const float deltaX)
{
    const auto  spParticle = std::make_shared<bf::Particle>(mass, energy);

    while (spParticle->KineticEnergy() > 0.f)
        propagator.Propagate(spParticle, deltaX);

    return spParticle;
}

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

TCanvas * PlotForParticles(const std::string &name, const std::function<TGraph *(const bf::Propagator &, const std::shared_ptr<bf::Particle> &)> &graphGetter,
    const std::string &yaxisLabel, const bool logY, const std::shared_ptr<bf::Particle> &spMuon, const std::shared_ptr<bf::Particle> &spProton,
    const std::shared_ptr<bf::Particle> &spPion, const std::shared_ptr<bf::Particle> &spKaon, Int_t *palette, const int markerStyle, const bf::Propagator &propagator)
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
    pMultiGraph->GetXaxis()->SetLimits(0.f, 200.);

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

TCanvas * PlotForParticle(const std::string &name, const std::function<TGraph *(const bf::Propagator &, const std::shared_ptr<bf::Particle> &)> &graphGetter,
    const std::string &yaxisLabel, const bool logY, const std::shared_ptr<bf::Particle> &spParticle, Int_t color, const int markerStyle, const bf::Propagator &propagator)
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
    pParticleResponse->GetXaxis()->SetLimits(0.f, 200.);

    // Redraw.
    pParticleResponse->Draw("ap");
    pCanvas->Update();

    return pCanvas;
}

void TestGeneration()
{
    const float deltaX = 0.3f;
    const int markerStyle = 6; // 1, 6 are 7 are dots of increasing size

    // Propagate some particles
    const auto propagator = bf::Propagator(bf::DetectorHelper::GetMicroBooNEDetector());

    const auto spMuon   = bf::ParticleHelper::GetMuon(500.f);
    const auto spProton = bf::ParticleHelper::GetProton(500.f);
    const auto spPion   = bf::ParticleHelper::GetChargedPion(500.f);
    const auto spKaon   = bf::ParticleHelper::GetChargedKaon(500.f);

    propagator.PropagateUntilStopped(deltaX, {spMuon, spProton, spPion, spKaon});

    // Plot the results
    Int_t palette[8];
    palette[0] = TColor::GetColor("#1B9E77");
    palette[1] = TColor::GetColor("#D95F02");
    palette[2] = TColor::GetColor("#7570B3");
    palette[3] = TColor::GetColor("#E7298A");
    palette[4] = TColor::GetColor("#66A61E");
    palette[5] = TColor::GetColor("#E6AB02");
    palette[6] = TColor::GetColor("#A6761D");
    palette[7] = TColor::GetColor("#666666");

    PlotForParticle("dEdx_muon", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spMuon, palette[0], markerStyle, propagator);
    PlotForParticle("dQdx_muon", PlotParticledQdx, "\\mathrm{d}Q/\\mathrm{d}x \\text{ (ADC/cm)}", false, spMuon, palette[0], markerStyle, propagator);

    PlotForParticle("dEdx_proton", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spProton, palette[1], markerStyle, propagator);
    PlotForParticle("dQdx_proton", PlotParticledQdx, "\\mathrm{d}Q/\\mathrm{d}x \\text{ (ADC/cm)}", false, spProton, palette[1], markerStyle, propagator);

    PlotForParticle("dEdx_pion", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spPion, palette[2], markerStyle, propagator);
    PlotForParticle("dQdx_pion", PlotParticledQdx, "\\mathrm{d}Q/\\mathrm{d}x \\text{ (ADC/cm)}", false, spPion, palette[2], markerStyle, propagator);

    PlotForParticle("dEdx_kaon", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spKaon, palette[3], markerStyle, propagator);
    PlotForParticle("dQdx_kaon", PlotParticledQdx, "\\mathrm{d}Q/\\mathrm{d}x \\text{ (ADC/cm)}", false, spKaon, palette[3], markerStyle, propagator);

    PlotForParticles("beta", PlotParticleBetas, "\\beta", false, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);
    PlotForParticles("energy", PlotParticleEnergies, "T \\text{ (MeV)}", false, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);
    PlotForParticles("dEdx", PlotParticledEdx, "-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}", true, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);
    PlotForParticles("dQdx", PlotParticledQdx, "\\mathrm{d}Q/\\mathrm{d}x \\text{ (ADC/cm)}", false, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);

    TCanvas * pKappaCanvas = PlotForParticles("kappa", PlotParticleKappas, "\\kappa", true, spMuon, spProton, spPion, spKaon, palette, markerStyle, propagator);

    TText *pTextLandau = new TText(11., 4.e-4, "Landau");
    pTextLandau->SetTextFont(42);
    pTextLandau->SetTextSize(0.04);
    pTextLandau->Draw();

    TLine *pLineLandau = new TLine(0., 0.01, pKappaCanvas->GetUxmax(), 0.01);
    pLineLandau->SetLineColor(kBlack);
    pLineLandau->SetLineWidth(2);
    pLineLandau->SetLineStyle(2);
    pLineLandau->Draw();

    TText *pTextConvolution = new TText(11., 4.e-2, "Convolution");
    pTextConvolution->SetTextFont(42);
    pTextConvolution->SetTextSize(0.04);
    pTextConvolution->Draw();

    TLine *pLineConvolution = new TLine(0., 0.3, pKappaCanvas->GetUxmax(), 0.3);
    pLineConvolution->SetLineColor(kBlack);
    pLineConvolution->SetLineWidth(2);
    pLineConvolution->SetLineStyle(2);
    pLineConvolution->Draw();

    TText *pTextLognormal = new TText(11., 1.2, "Lognormal");
    pTextLognormal->SetTextFont(42);
    pTextLognormal->SetTextSize(0.04);
    pTextLognormal->Draw();

    TLine *pLineLogNormal = new TLine(0., 10., pKappaCanvas->GetUxmax(), 10.);
    pLineLogNormal->SetLineColor(kBlack);
    pLineLogNormal->SetLineWidth(2);
    pLineLogNormal->SetLineStyle(2);
    pLineLogNormal->Draw();

    TText *pTextNormal = new TText(11., 50., "Normal");
    pTextNormal->SetTextFont(42);
    pTextNormal->SetTextSize(0.04);
    pTextNormal->Draw();
    
    Pause();
}

void TestFiltering()
{
    // Uniform pior
    auto distribution = bf::ParticleFilter::ParticleDistribution{};

    const float minMass = 50.f;
    const float maxMass = 200.f;
    const float minEnergy = 100.f;
    const float maxEnergy = 1000.f;
    const unsigned int numSteps = 1000;
/*
    for (unsigned int i = 0; i < numSteps; ++i)
    {
        const float mass = minMass + static_cast<float>(i) * (maxMass - minMass) / static_cast<float>(numSteps - 1U);

        for (unsigned int j = 0; j < numSteps; ++j)
        {
            const float energy = minEnergy + static_cast<float>(j) * (maxEnergy - minEnergy) / static_cast<float>(numSteps - 1U);
            distribution.emplace_back(1.f, bf::ParticleHelper::GetParticle(mass, energy));
        }
    }*/

    for (unsigned int i = 0; i < numSteps; ++i)
    {
        const float energy = minEnergy + static_cast<float>(i) * (maxEnergy - minEnergy) / static_cast<float>(numSteps - 1U);
        distribution.emplace_back(1.f, bf::ParticleHelper::GetParticle(bf::PhysicalConstants::m_protonMass, energy));
    }

    // Propagate a muon
    const float deltaX = 0.03f;

    const auto detector   = bf::DetectorHelper::GetMicroBooNEDetector();
    const auto propagator = bf::Propagator{detector};
    const auto spProton   = bf::ParticleHelper::GetProton(750.f);
    propagator.PropagateUntilStopped(deltaX, spProton);

    // Try to reconstruct the muon mass
    TCanvas *pCanvas = new TCanvas("c1", "c1", 800, 600);
    auto filter = bf::ParticleFilter{detector, distribution, deltaX};

    for (const bf::ParticleState &state : spProton->GetHistory())
    {
        distribution = filter.Filter(bf::ObservedParticleState{state.GetPosition(), state.GetdEdx()});
        
        TH1F *pHistogram = new TH1F("proton_energy", "proton_energy", 50, 0.f, 1000.f);
        float avgKineticEnergy = 0.f;

        for (const auto &pair : distribution)
        {
            pHistogram->Fill(pair.second->KineticEnergy());
            avgKineticEnergy += pair.first * pair.second->KineticEnergy();
        }

        std::cout << "True KE          = " << state.GetKineticEnergy() << std::endl;
        std::cout << "Avg predicted KE = " << avgKineticEnergy << std::endl;

        pHistogram->Draw();
        pCanvas->Update();

        Pause();

        delete pHistogram;
    }
}

int main()
{  
    const bool testGeneration = false;
    const bool testFiltering = true;

    TApplication *pApplication = Initialize();

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
    
    return 0;
}
