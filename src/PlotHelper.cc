/**
 *  @file   bethe-faster/src/PlotHelper.cc
 *
 *  @brief  Implementation of the plot helper class.
 *
 *  $Log: $
 */

#include "PlotHelper.h"

#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include <unistd.h>
#include <fcntl.h>

namespace bf
{

TApplication *PlotHelper::InitializeApplication()
{
    TApplication *pApplication = nullptr;

    if (gApplication)
        pApplication = gApplication;

    else
    {
        int   argc   = 0;
        char *argv[1UL];
        argv[0] = nullptr;
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
    pStyle->SetPadRightMargin(0.1f);
    pStyle->SetPadGridX(0);
    pStyle->SetPadGridY(0);
    pStyle->SetPadTickX(1);
    pStyle->SetPadTickY(1);
    pStyle->SetHistFillColor(0);
    pStyle->SetFrameBorderMode(0);

    gROOT->SetStyle("BetheFasterStyle");
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph PlotHelper::PlotParticledEdx(const std::shared_ptr<Particle> &spParticle)
{
    std::vector<double> positions, responses;
    const auto &        history = spParticle->GetHistory();

    for (const auto &spState : history)
    {
        positions.push_back(spState->GetResidualRange());
        responses.push_back(spState->GetKineticEnergy());
    }

    return TGraph{static_cast<Int_t>(history.size()), positions.data(), responses.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph PlotHelper::PlotParticleKappa(const Propagator &propagator, const std::shared_ptr<Particle> &spParticle)
{
    std::vector<double> positions, kappas;
    const auto &       history = spParticle->GetHistory();

    for (const auto &spState : history)
    {
        positions.push_back(spState->GetResidualRange());
        kappas.push_back(propagator.CalculateKappa(spParticle->Mass(), spState->GetKineticEnergy(), spState->Getdx()));
    }

    return TGraph{static_cast<Int_t>(history.size()), positions.data(), kappas.data()};
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
