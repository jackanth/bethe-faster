/**
 *  @file   bethe-faster/src/Propagator.cc
 *
 *  @brief  Implementation of the propagator class.
 *
 *  $Log: $
 */

#include "Propagator.h"
#include "PhysicalConstants.h"

#include "gsl/gsl_randist.h"
#include "Math/VavilovAccurate.h"
#include "Math/VavilovFast.h"
#include "Math/DistFunc.h"

#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>

namespace bf
{

Propagator::Propagator(Detector detector) :
    m_detector(std::move_if_noexcept(detector)),
    m_cBar(0.),
    m_xiPartial(0.),
    m_pGenerator(nullptr),
    m_pRandom(nullptr)
{
    if (m_detector.m_plasmaEnergy < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize detector because plasma energy value was too small");

    m_cBar = 2. * std::log(m_detector.m_avgIonizationEnergy / m_detector.m_plasmaEnergy) + 1.;

    if (2. * m_detector.m_atomicMass < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize detector because the atomic mass was too small");

    m_xiPartial = PhysicalConstants::m_kCoefficient * static_cast<double>(m_detector.m_atomicNumber) *
                  m_detector.m_density / (2. * m_detector.m_atomicMass);

    // Set up the random number generators
    const gsl_rng_type * pType;
    gsl_rng_env_setup();
    pType = gsl_rng_default;
    m_pGenerator = gsl_rng_alloc(pType);
    gsl_rng_set(m_pGenerator, static_cast<std::uint32_t>(std::time(nullptr)));

    m_pRandom = new TRandom(static_cast<std::uint32_t>(std::time(nullptr)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

Propagator::~Propagator()
{
    gsl_rng_free(m_pGenerator);

    if (m_pRandom)
        delete m_pRandom;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::Propagate(const std::shared_ptr<Particle> &spParticle, const double deltaX) const
{   
    const double beta  = ParticleHelper::GetParticleBeta(spParticle);
    const double beta2 = beta * beta;

    if (beta < 0.001)
    {
        spParticle->Increment(deltaX, -spParticle->KineticEnergy());
        return 0.f;
    }

    const double xi = this->Xi(beta2, deltaX);
    const double maxEnergyTransfer = this->GetMaxEnergyTransfer(beta, spParticle->Mass());
    const double expectedLoss      = this->GetExpectedEnergyLoss(beta, beta2, xi, maxEnergyTransfer);
    const double kappa             = xi / maxEnergyTransfer;
    
    const double actualLoss = this->SampleDeltaE(expectedLoss, kappa, beta2, xi);
    const double deltaE = std::min(-actualLoss, 0.); // particle may not gain energy

    spParticle->Increment(deltaX, deltaE);
    return deltaE / deltaX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::CalculateObservationProbability(const std::shared_ptr<Particle> &spParticle, const double observeddEdx, const double observeddx) const
{
    return this->CalculateTransitionProbability(observeddEdx * observeddx, spParticle->Mass(), spParticle->KineticEnergy(), observeddx);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::CalculateTransitionProbability(const double deltaE, const double mass, const double previousEnergy, const double deltaX) const
{
    const double energyLoss = -deltaE;
    const double previousBeta  = ParticleHelper::GetParticleBeta(mass, previousEnergy);
    const double previousBeta2 = previousBeta * previousBeta;

    if (previousBeta2 < std::numeric_limits<double>::epsilon())
        return 0.;

    const double previousXi        = this->Xi(previousBeta2, deltaX);
    const double maxEnergyTransfer = this->GetMaxEnergyTransfer(previousBeta, mass);
    const double expectedLoss      = this->GetExpectedEnergyLoss(previousBeta, previousBeta2, previousXi, maxEnergyTransfer);
    const double previousKappa     = previousXi / maxEnergyTransfer;

    if (previousKappa > 10.) // Gaussian approximation
    {
        const double sigma = previousXi * std::sqrt((1. - 0.5 * previousBeta2) / previousKappa);
        return gsl_ran_gaussian_pdf(energyLoss - expectedLoss, sigma);
    }

    const double lambda = (energyLoss - expectedLoss) / previousXi - 1. - previousBeta2 -  std::log(previousKappa) + PhysicalConstants::m_eulerConstant;

    if (previousKappa > 0.01) // Vavilov
        return this->GetVavilovProbabilityDensity(previousKappa, previousBeta2, lambda);

    // kappa < 0.01; Landau approximation
    return gsl_ran_landau_pdf(lambda);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::DensityCorrection(const double beta) const
{  
    const double x =  std::log10(beta) - 0.5 * std::log10(1. - beta * beta);;

    if (x > m_detector.m_sternheimerX1)
        return 2. * std::log(10.) * x - m_cBar;

    if (x > m_detector.m_sternheimerX0)
        return 2. * std::log(10.) * x - m_cBar + m_detector.m_sternheimerA * std::pow(m_detector.m_sternheimerX1 - x, m_detector.m_sternheimerK);

    return m_detector.m_sternheimerDelta0 == 0. ? 0. : m_detector.m_sternheimerDelta0 * std::pow(10., 2. * (x - m_detector.m_sternheimerX0));
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::SampleDeltaE(const double meanEnergyLoss, const double kappa, const double beta2, const double xi) const
{    
    if (kappa > 10.) // Gaussian approximation
    {
        const double sigma = xi * std::sqrt((1. - 0.5 * beta2) / kappa);
        return meanEnergyLoss + gsl_ran_gaussian(m_pGenerator, sigma);
    }

    if (kappa > 0.01) // Vavilov
        return meanEnergyLoss + xi * (this->SampleFromVavilovDistribution(kappa, beta2) + 1. + beta2 + std::log(kappa) - PhysicalConstants::m_eulerConstant);

    // kappa < 0.01; Landau approximation
    const double lambdaBar = -(1. - PhysicalConstants::m_eulerConstant) - beta2 - std::log(kappa);
    const double lambdaMax = 0.60715 + 1.1934 * lambdaBar + (0.67794 + 0.052382 * lambdaBar) * std::exp(0.94753 + 0.74442 * lambdaBar);

    double lambda = 0.f;

    do lambda = gsl_ran_landau(m_pGenerator); while (lambda > lambdaMax);

    // Note that the Landau mode is at approximately -0.22278
    return meanEnergyLoss + xi * (lambda + 1. + beta2 + std::log(kappa) - PhysicalConstants::m_eulerConstant);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::GetMaxEnergyTransfer(const double beta, const double mass) const
{
    const double beta2 = beta * beta;
    const double gamma = 1. / std::sqrt(1 - beta2);

    const double numerator   = 2. * PhysicalConstants::m_electronMass * beta2 * gamma * gamma;
    const double denominator = 1. + 2. * gamma * PhysicalConstants::m_electronMass / mass +
                              (PhysicalConstants::m_electronMass / mass) * (PhysicalConstants::m_electronMass / mass);

    return numerator / denominator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::GetExpectedEnergyLoss(const double beta, const double beta2, const double xi, const double maxEnergyLoss) const 
{
    const double lnArg =
        2. * PhysicalConstants::m_electronMass * beta2 * maxEnergyLoss /
        ((1. - beta2) * m_detector.m_avgIonizationEnergy * m_detector.m_avgIonizationEnergy * 1.e-12);

    return xi * (std::log(lnArg) - 2. * beta2 - this->DensityCorrection(beta));
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::SampleFromVavilovDistribution(const double kappa, const double beta2, const double maxError, const std::size_t maxIterations) const
{
    const auto   vavilovDist   = ROOT::Math::VavilovFast{kappa, beta2};
    double lambda     = vavilovDist.Mode();
    double error      = 1.;
    std::size_t count = 0UL;
    double uniformSample = m_pRandom->Uniform(0., 1.);

    while ((std::abs(error) > maxError) && (count++ < maxIterations))
    {
        const double cdf = vavilovDist.Cdf(lambda);
        const double pdf = vavilovDist.Pdf(lambda);

        if (pdf < std::numeric_limits<double>::epsilon())
        {
            uniformSample = m_pRandom->Uniform(0., 1.);
            lambda        = vavilovDist.Mode();
            error         = 1.;
            continue;
        }

        error = (cdf - uniformSample) / pdf;
        lambda = lambda - error;
    }

    return lambda;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double Propagator::GetVavilovProbabilityDensity(const double kappa, const double beta2, const double lambda) const
{
    const auto vavilovDist = ROOT::Math::VavilovFast{kappa, beta2};
    return vavilovDist.Pdf(lambda);
}

} // namespace bf

