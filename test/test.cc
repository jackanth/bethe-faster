/**
 *  @file   bethe-faster/test/Test.cc
 *
 *  @brief  Implementation of the tests.
 *
 *  $Log: $
 */

#include "Test.h"
#include <iostream>

int main()
{
    try
    {
        bf::PlotHelper::SetGlobalPlotStyle();
        const auto detector = bf::DetectorHelper::GetMicroBooNEDetector();

        std::cout << "Testing generation step" << std::endl;
        // TestGeneration(detector);

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
    auto spMuon        = bf::ParticleHelper::GetStoppedMuon();
    auto spProton      = bf::ParticleHelper::GetStoppedProton();
    auto spChargedPion = bf::ParticleHelper::GetStoppedChargedPion();
    auto spChargedKaon = bf::ParticleHelper::GetStoppedChargedKaon();

    // Propagate particles
    while (spMuon->KineticEnergy() < maxEnergy)
        propagator.PropagateBackwards(spMuon, deltaX, mode);

    while (spProton->KineticEnergy() < maxEnergy)
        propagator.PropagateBackwards(spProton, deltaX, mode);

    while (spChargedPion->KineticEnergy() < maxEnergy)
        propagator.PropagateBackwards(spChargedPion, deltaX, mode);

    while (spChargedKaon->KineticEnergy() < maxEnergy)
        propagator.PropagateBackwards(spChargedKaon, deltaX, mode);

    // Draw plots
    const std::string nameSuffix = "_" + modeName + "_" + std::to_string(deltaX) + "cm_" + std::to_string(maxEnergy) + "MeV.eps";

    if (mode == bf::Propagator::PROPAGATION_MODE::STOCHASTIC)
    {
        bf::PlotHelper::PlotParticleEnergyMarkers(spMuon, 0UL)->SaveAs(("muon_energy" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleEnergyMarkers(spChargedPion, 1UL)->SaveAs(("charged_pion_energy" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleEnergyMarkers(spChargedKaon, 2UL)->SaveAs(("charged_kaon_energy" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleEnergyMarkers(spProton, 3UL)->SaveAs(("proton_energy" + nameSuffix).c_str());

        bf::PlotHelper::PlotParticledEdxVersusXMarkers(spMuon, 0UL)->SaveAs(("muon_dEdx_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusXMarkers(spChargedPion, 1UL)->SaveAs(("charged_pion_dEdx_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusXMarkers(spChargedKaon, 2UL)->SaveAs(("charged_kaon_dEdx_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusXMarkers(spProton, 3UL)->SaveAs(("proton_dEdx_versus_x" + nameSuffix).c_str());

        bf::PlotHelper::PlotParticledEdxVersusTMarkers(spMuon, 0UL)->SaveAs(("muon_dEdx_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusTMarkers(spChargedPion, 1UL)->SaveAs(("charged_pion_dEdx_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusTMarkers(spChargedKaon, 2UL)->SaveAs(("charged_kaon_dEdx_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusTMarkers(spProton, 3UL)->SaveAs(("proton_dEdx_versus_T" + nameSuffix).c_str());

        bf::PlotHelper::PlotParticleKappaVersusXMarkers(propagator, spMuon, 0UL)->SaveAs(("muon_kappa_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusXMarkers(propagator, spChargedPion, 1UL)->SaveAs(("charged_pion_kappa_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusXMarkers(propagator, spChargedKaon, 2UL)->SaveAs(("charged_kaon_kappa_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusXMarkers(propagator, spProton, 3UL)->SaveAs(("proton_kappa_versus_x" + nameSuffix).c_str());

        bf::PlotHelper::PlotParticleKappaVersusTMarkers(propagator, spMuon, 0UL)->SaveAs(("muon_kappa_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusTMarkers(propagator, spChargedPion, 1UL)->SaveAs(("charged_pion_kappa_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusTMarkers(propagator, spChargedKaon, 2UL)->SaveAs(("charged_kaon_kappa_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusTMarkers(propagator, spProton, 3UL)->SaveAs(("proton_kappa_versus_T" + nameSuffix).c_str());
    }

    else
    {
        bf::PlotHelper::PlotParticleEnergyLine(spMuon, 0UL)->SaveAs(("muon_energy" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleEnergyLine(spChargedPion, 1UL)->SaveAs(("charged_pion_energy" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleEnergyLine(spChargedKaon, 2UL)->SaveAs(("charged_kaon_energy" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleEnergyLine(spProton, 3UL)->SaveAs(("proton_energy" + nameSuffix).c_str());

        bf::PlotHelper::PlotParticledEdxVersusXLine(spMuon, 0UL)->SaveAs(("muon_dEdx_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusXLine(spChargedPion, 1UL)->SaveAs(("charged_pion_dEdx_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusXLine(spChargedKaon, 2UL)->SaveAs(("charged_kaon_dEdx_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusXLine(spProton, 3UL)->SaveAs(("proton_dEdx_versus_x" + nameSuffix).c_str());

        bf::PlotHelper::PlotParticledEdxVersusTLine(spMuon, 0UL)->SaveAs(("muon_dEdx_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusTLine(spChargedPion, 1UL)->SaveAs(("charged_pion_dEdx_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusTLine(spChargedKaon, 2UL)->SaveAs(("charged_kaon_dEdx_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticledEdxVersusTLine(spProton, 3UL)->SaveAs(("proton_dEdx_versus_T" + nameSuffix).c_str());

        bf::PlotHelper::PlotParticleKappaVersusXLine(propagator, spMuon, 0UL)->SaveAs(("muon_kappa_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusXLine(propagator, spChargedPion, 1UL)->SaveAs(("charged_pion_kappa_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusXLine(propagator, spChargedKaon, 2UL)->SaveAs(("charged_kaon_kappa_versus_x" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusXLine(propagator, spProton, 3UL)->SaveAs(("proton_kappa_versus_x" + nameSuffix).c_str());

        bf::PlotHelper::PlotParticleKappaVersusTLine(propagator, spMuon, 0UL)->SaveAs(("muon_kappa_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusTLine(propagator, spChargedPion, 1UL)->SaveAs(("charged_pion_kappa_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusTLine(propagator, spChargedKaon, 2UL)->SaveAs(("charged_kaon_kappa_versus_T" + nameSuffix).c_str());
        bf::PlotHelper::PlotParticleKappaVersusTLine(propagator, spProton, 3UL)->SaveAs(("proton_kappa_versus_T" + nameSuffix).c_str());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TestFiltering(const bf::Detector &detector)
{
    // Define the particles
    auto spMuon        = bf::ParticleHelper::GetStoppedMuon();
    auto spChargedPion = bf::ParticleHelper::GetStoppedChargedPion();
    auto spChargedKaon = bf::ParticleHelper::GetStoppedChargedKaon();
    auto spProton      = bf::ParticleHelper::GetStoppedProton();

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
    spParticle->Reset();

    while (spParticle->KineticEnergy() < maxEnergy)
        propagator.PropagateBackwards(spParticle, deltaX);

    const auto particleFilter = bf::ParticleFilter{detector};

    auto muonHypothesis = bf::MassHypothesis{bf::PhysicalConstants::m_muonMass, bf::ParticleHelper::GetMinimumMuonEnergy(), 1000UL};
    auto chargedPionHypothesis =
        bf::MassHypothesis{bf::PhysicalConstants::m_chargedPionMass, bf::ParticleHelper::GetMinimumChargedPionEnergy(), 1000UL};
    auto chargedKaonHypothesis =
        bf::MassHypothesis{bf::PhysicalConstants::m_chargedKaonMass, bf::ParticleHelper::GetMinimumChargedKaonEnergy(), 1000UL};
    auto protonHypothesis = bf::MassHypothesis{bf::PhysicalConstants::m_protonMass, bf::ParticleHelper::GetMinimumProtonEnergy(), 1000UL};

    const auto   muonDistributionHistory = particleFilter.Filter(spParticle->GetHistory(), muonHypothesis);
    const double muonLogLikelihood       = bf::ParticleFilter::CalculateMarginalLikelihood(muonDistributionHistory);

    const auto   chargedPionDistributionHistory = particleFilter.Filter(spParticle->GetHistory(), chargedPionHypothesis);
    const double chargedPionLogLikelihood       = bf::ParticleFilter::CalculateMarginalLikelihood(chargedPionDistributionHistory);

    const auto   chargedKaonDistributionHistory = particleFilter.Filter(spParticle->GetHistory(), chargedKaonHypothesis);
    const double chargedKaonLogLikelihood       = bf::ParticleFilter::CalculateMarginalLikelihood(chargedKaonDistributionHistory);

    const auto   protonDistributionHistory = particleFilter.Filter(spParticle->GetHistory(), protonHypothesis);
    const double protonLogLikelihood       = bf::ParticleFilter::CalculateMarginalLikelihood(protonDistributionHistory);

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
}
