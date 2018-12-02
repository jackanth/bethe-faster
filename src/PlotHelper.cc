/**
 *  @file   bethe-faster/src/PlotHelper.cc
 *
 *  @brief  Implementation of the plot helper class.
 *
 *  $Log: $
 */

#include "PlotHelper.h"

#include "TAxis.h"
#include "TROOT.h"
#include "TSystem.h"

#include <fcntl.h>
#include <unistd.h>

namespace bf
{

std::size_t g_canvasCount = 0UL; ///< The global canvas counter

TApplication *PlotHelper::InitializeApplication()
{
    TApplication *pApplication = nullptr;

    if (gApplication)
        pApplication = gApplication;

    else
    {
        int   argc = 0;
        char *argv[1UL];
        argv[0]      = nullptr;
        pApplication = new TApplication("BetheFaster", &argc, const_cast<char **>(argv));
        pApplication->SetReturnFromRun(kTRUE);
    }

    return pApplication;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void PlotHelper::SetGlobalPlotStyle()
{
    TStyle *pStyle = new TStyle("BetheFasterStyle", "BetheFasterStyle");

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
    pStyle->SetPadRightMargin(0.08f);
    pStyle->SetPadGridX(1);
    pStyle->SetPadGridY(1);
    pStyle->SetPadTickX(1);
    pStyle->SetPadTickY(1);
    pStyle->SetHistFillColor(0);
    pStyle->SetFrameBorderMode(0);

    Int_t palette[8];
    palette[0] = static_cast<std::int16_t>(TColor::GetColor("#1B9E77"));
    palette[1] = static_cast<std::int16_t>(TColor::GetColor("#D95F02"));
    palette[2] = static_cast<std::int16_t>(TColor::GetColor("#7570B3"));
    palette[3] = static_cast<std::int16_t>(TColor::GetColor("#E7298A"));
    palette[4] = static_cast<std::int16_t>(TColor::GetColor("#66A61E"));
    palette[5] = static_cast<std::int16_t>(TColor::GetColor("#E6AB02"));
    palette[6] = static_cast<std::int16_t>(TColor::GetColor("#A6761D"));
    palette[7] = static_cast<std::int16_t>(TColor::GetColor("#666666"));
    pStyle->SetPalette(8, palette, 0.8f);

    gROOT->SetStyle("BetheFasterStyle");
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph PlotHelper::GetParticleEnergyGraph(const std::shared_ptr<Particle> &spParticle, const bool useResidualRange)
{
    std::vector<double> positions, responses;
    const auto &        history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    const double maxResidualRange = history.back()->GetResidualRange();

    for (const auto &spState : history)
    {
        positions.push_back(useResidualRange ? spState->GetResidualRange() : maxResidualRange - spState->GetResidualRange());
        responses.push_back(spState->GetKineticEnergy());
    }

    return TGraph{static_cast<Int_t>(history.size()), positions.data(), responses.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticleEnergyLine(
    const std::shared_ptr<Particle> &spParticle, const unsigned int colour, const std::int16_t lineWidth, const bool useResidualRange)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticleEnergyGraph(spParticle, useResidualRange);

    // Formatting
    graph.GetYaxis()->SetTitle("T \\text{ (MeV)}");
    graph.GetXaxis()->SetTitle(useResidualRange ? "\\text{Residual range (cm)}" : "x\\text{ (cm)}");
    graph.SetLineColor(PlotHelper::GetSchemeColour(colour));
    graph.SetLineWidth(lineWidth);

    graph.DrawClone("AL");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticleEnergyMarkers(
    const std::shared_ptr<Particle> &spParticle, const unsigned int colour, const std::int16_t markerStyle, const bool useResidualRange)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticleEnergyGraph(spParticle, useResidualRange);

    // Formatting
    graph.GetYaxis()->SetTitle("T \\text{ (MeV)}");
    graph.GetXaxis()->SetTitle(useResidualRange ? "\\text{Residual range (cm)}" : "x\\text{ (cm)}");
    graph.SetMarkerColor(PlotHelper::GetSchemeColour(colour));
    graph.SetMarkerStyle(markerStyle);

    graph.DrawClone("AP");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph PlotHelper::GetParticledEdxVersusXGraph(const std::shared_ptr<Particle> &spParticle, const bool useResidualRange)
{
    const auto &history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    std::vector<double> positions, responses;
    const double        maxResidualRange = history.back()->GetResidualRange();

    for (const auto &spState : history)
    {
        positions.push_back(useResidualRange ? spState->GetResidualRange() : maxResidualRange - spState->GetResidualRange());
        responses.push_back(-spState->GetdEdx());
    }

    return TGraph{static_cast<Int_t>(history.size()), positions.data(), responses.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticledEdxVersusXLine(
    const std::shared_ptr<Particle> &spParticle, const unsigned int colour, const std::int16_t lineWidth, const bool useResidualRange)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticledEdxVersusXGraph(spParticle, useResidualRange);

    // Formatting
    graph.GetYaxis()->SetTitle("-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
    graph.GetXaxis()->SetTitle(useResidualRange ? "\\text{Residual range (cm)}" : "x\\text{ (cm)}");
    graph.SetLineColor(PlotHelper::GetSchemeColour(colour));
    graph.SetLineWidth(lineWidth);
    pCanvas->SetLogy();

    graph.DrawClone("AL");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticledEdxVersusXMarkers(
    const std::shared_ptr<Particle> &spParticle, const unsigned int colour, const std::int16_t markerStyle, const bool useResidualRange)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticledEdxVersusXGraph(spParticle, useResidualRange);

    // Formatting
    graph.GetYaxis()->SetTitle("-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
    graph.GetXaxis()->SetTitle(useResidualRange ? "\\text{Residual range (cm)}" : "x\\text{ (cm)}");
    graph.SetMarkerColor(PlotHelper::GetSchemeColour(colour));
    graph.SetMarkerStyle(markerStyle);
    pCanvas->SetLogy();

    graph.DrawClone("AP");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph PlotHelper::GetParticledEdxVersusTGraph(const std::shared_ptr<Particle> &spParticle)
{
    const auto &history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    std::vector<double> energies, responses;

    for (const auto &spState : history)
    {
        energies.push_back(spState->GetKineticEnergy());
        responses.push_back(-spState->GetdEdx());
    }

    return TGraph{static_cast<Int_t>(history.size()), energies.data(), responses.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticledEdxVersusTLine(const std::shared_ptr<Particle> &spParticle, const unsigned int colour, const std::int16_t lineWidth)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticledEdxVersusTGraph(spParticle);

    // Formatting
    graph.GetYaxis()->SetTitle("-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
    graph.GetXaxis()->SetTitle("T \\text{ (MeV)}");
    graph.SetLineColor(PlotHelper::GetSchemeColour(colour));
    graph.SetLineWidth(lineWidth);
    pCanvas->SetLogy();

    graph.DrawClone("AL");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticledEdxVersusTMarkers(const std::shared_ptr<Particle> &spParticle, const unsigned int colour, const std::int16_t markerStyle)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticledEdxVersusTGraph(spParticle);

    // Formatting
    graph.GetYaxis()->SetTitle("-\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
    graph.GetXaxis()->SetTitle("T \\text{ (MeV)}");
    graph.SetMarkerColor(PlotHelper::GetSchemeColour(colour));
    graph.SetMarkerStyle(markerStyle);
    pCanvas->SetLogy();

    graph.DrawClone("AP");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph PlotHelper::GetParticleKappaVersusXGraph(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle, const bool useResidualRange)
{
    const auto &history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    std::vector<double> positions, kappas;
    const double        maxResidualRange = history.back()->GetResidualRange();

    for (const auto &spState : history)
    {
        positions.push_back(useResidualRange ? spState->GetResidualRange() : maxResidualRange - spState->GetResidualRange());
        kappas.push_back(propagator.CalculateKappa(spParticle->Mass(), spState->GetKineticEnergy(), spState->Getdx()));
    }

    return TGraph{static_cast<Int_t>(history.size()), positions.data(), kappas.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticleKappaVersusXLine(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle,
    const unsigned int colour, const std::int16_t lineWidth, const bool useResidualRange)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticleKappaVersusXGraph(propagator, spParticle, useResidualRange);

    // Formatting
    graph.GetYaxis()->SetTitle("\\kappa");
    graph.GetXaxis()->SetTitle(useResidualRange ? "\\text{Residual range (cm)}" : "x\\text{ (cm)}");
    graph.SetLineColor(PlotHelper::GetSchemeColour(colour));
    graph.SetLineWidth(lineWidth);
    pCanvas->SetLogy();

    graph.DrawClone("AL");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticleKappaVersusXMarkers(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle,
    const unsigned int colour, const std::int16_t markerStyle, const bool useResidualRange)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticleKappaVersusXGraph(propagator, spParticle, useResidualRange);

    // Formatting
    graph.GetYaxis()->SetTitle("\\kappa");
    graph.GetXaxis()->SetTitle(useResidualRange ? "\\text{Residual range (cm)}" : "x\\text{ (cm)}");
    graph.SetMarkerColor(PlotHelper::GetSchemeColour(colour));
    graph.SetMarkerStyle(markerStyle);
    pCanvas->SetLogy();

    graph.DrawClone("AP");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph PlotHelper::GetParticleKappaVersusTGraph(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle)
{
    const auto &history = spParticle->GetHistory();

    if (history.empty())
        throw std::runtime_error{"Cannot plot empty particle history"};

    std::vector<double> positions, kappas;

    for (const auto &spState : history)
    {
        positions.push_back(spState->GetKineticEnergy());
        kappas.push_back(propagator.CalculateKappa(spParticle->Mass(), spState->GetKineticEnergy(), spState->Getdx()));
    }

    return TGraph{static_cast<Int_t>(history.size()), positions.data(), kappas.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticleKappaVersusTLine(
    const Propagator &propagator, const std::shared_ptr<Particle> &spParticle, const unsigned int colour, const std::int16_t lineWidth)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticleKappaVersusTGraph(propagator, spParticle);

    // Formatting
    graph.GetYaxis()->SetTitle("\\kappa");
    graph.GetXaxis()->SetTitle("T \\text{ (MeV)}");
    graph.SetLineColor(PlotHelper::GetSchemeColour(colour));
    graph.SetLineWidth(lineWidth);
    pCanvas->SetLogy();

    graph.DrawClone("AL");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TCanvas *PlotHelper::PlotParticleKappaVersusTMarkers(
    const Propagator &propagator, const std::shared_ptr<Particle> &spParticle, const unsigned int colour, const std::int16_t markerStyle)
{
    const auto     canvasName = "bfCanvas" + std::to_string(g_canvasCount++);
    TCanvas *const pCanvas    = new TCanvas{canvasName.c_str(), canvasName.c_str(), 800, 600};
    TGraph         graph      = PlotHelper::GetParticleKappaVersusTGraph(propagator, spParticle);

    // Formatting
    graph.GetYaxis()->SetTitle("\\kappa");
    graph.GetXaxis()->SetTitle("T \\text{ MeV}");
    graph.SetMarkerColor(PlotHelper::GetSchemeColour(colour));
    graph.SetMarkerStyle(markerStyle);
    pCanvas->SetLogy();

    graph.DrawClone("AP");
    return pCanvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void PlotHelper::Pause()
{
    std::cout << "Press return to continue..." << std::endl;
    int flag = fcntl(1, F_GETFL, 0);

    int key = 0;

    while (true)
    {
        gSystem->ProcessEvents();
        fcntl(1, F_SETFL, flag | O_NONBLOCK);
        key = getchar();

        if ((key == '\n') || (key == '\r'))
            break;

        usleep(1000);
    }

    fcntl(1, F_SETFL, flag);
}

} // namespace bf
