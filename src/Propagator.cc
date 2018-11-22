/**
 *  @file   bethe-faster/src/Propagator.cc
 *
 *  @brief  Implementation of the propagator class.
 *
 *  $Log: $
 */

#include "Propagator.h"
#include "Math/VavilovAccurate.h"
#include "PhysicalConstants.h"

#include "gsl/gsl_randist.h"

#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>

namespace bf
{

Propagator::Propagator(Detector detector) :
    m_detector(std::move_if_noexcept(detector)),
    m_cBar(0.f),
    m_chargeMapC1(0.f),
    m_chargeMapC2(0.f),
    m_xiPartial(0.f),
    m_pGenerator(nullptr)
{
    if (m_detector.m_plasmaEnergy < std::numeric_limits<float>::epsilon())
        throw std::runtime_error("Could not initialize detector because plasma energy value was too small");

    m_cBar = 2.f * std::log(m_detector.m_avgIonizationEnergy / m_detector.m_plasmaEnergy) + 1.f;

    if (m_detector.m_thickness * m_detector.m_ionizationEnergy < std::numeric_limits<float>::epsilon())
        throw std::runtime_error("Could not initialize detector because the product of thickness and ionization energy was too small");

    m_chargeMapC1 = m_detector.m_birksA * m_detector.m_chargeCalibration / m_detector.m_ionizationEnergy;

    if (m_detector.m_thickness * m_detector.m_electricFieldStrength < std::numeric_limits<float>::epsilon())
        throw std::runtime_error(
            "Could not initialize detector because the product of thickness and electric field strength was too small");

    m_chargeMapC2 = m_detector.m_birksK / (m_detector.m_density * m_detector.m_electricFieldStrength);

    if (m_detector.m_thickness * m_detector.m_electricFieldStrength < std::numeric_limits<float>::epsilon())
        throw std::runtime_error(
            "Could not initialize detector because the product of thickness and electric field strength was too small");

    if (2. * m_detector.m_atomicMass < std::numeric_limits<float>::epsilon())
        throw std::runtime_error("Could not initialize detector because the atomic mass was too small");

    m_xiPartial = PhysicalConstants::m_kCoefficient * static_cast<float>(m_detector.m_atomicNumber) * m_detector.m_thickness *
                  m_detector.m_density / (2. * m_detector.m_atomicMass);

    // Set up the random number generator
    const gsl_rng_type * pType;
    gsl_rng_env_setup();
    pType = gsl_rng_default;
    m_pGenerator = gsl_rng_alloc(pType);
    gsl_rng_set(m_pGenerator, std::time(nullptr));
}

//------------------------------------------------------------------------------------------------------------------------------------------

Propagator::~Propagator()
{
    gsl_rng_free(m_pGenerator);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Propagator::Propagate(const std::shared_ptr<Particle> &spParticle, const float deltaX) const
{   
    const float beta              = ParticleHelper::GetParticleBeta(spParticle);
    const float beta2             = beta * beta;

    if (beta2 < std::numeric_limits<float>::epsilon())
        spParticle->Increment(deltaX, -spParticle->KineticEnergy());

    const float xi = this->Xi(beta);

    const float expectedLoss      = this->GetExpectedEnergyLoss(beta, beta2, xi);
    const float maxEnergyTransfer = this->GetMaxEnergyTransfer(beta, spParticle->Mass());
    const float kappa             = xi / maxEnergyTransfer;
    
    const float actualLoss = this->SampleDeltaE(expectedLoss, kappa, beta2, xi, maxEnergyTransfer);
    const float deltaE = std::min(-actualLoss * deltaX / m_detector.m_thickness, 0.f);

    spParticle->Increment(deltaX, deltaE);
    return deltaE / deltaX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Propagator::CalculateObservationProbability(const std::shared_ptr<Particle> &spParticle, const float observeddEdx, const float deltaX) const
{
    return this->CalculateTransitionProbability(observeddEdx * deltaX, spParticle->Mass(), spParticle->KineticEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Propagator::CalculateTransitionProbability(const float deltaE, const float mass, const float previousEnergy) const
{
    const float previousBeta  = ParticleHelper::GetParticleBeta(mass, previousEnergy);
    const float previousBeta2 = previousBeta * previousBeta;

    if (previousBeta2 < std::numeric_limits<float>::epsilon())
        return 0.f;

    const float previousXi        = this->Xi(previousBeta);
    const float expectedLoss      = this->GetExpectedEnergyLoss(previousBeta, previousBeta2, previousXi);
    const float maxEnergyTransfer = this->GetMaxEnergyTransfer(previousBeta, mass);
    const float previousKappa     = previousXi / maxEnergyTransfer;

    if (previousKappa > 10.f) // Gaussian approximation
    {
        const float sigma = std::sqrt((1.f - 0.5f * previousBeta2) / previousKappa);
        const float scaledLoss = (deltaE - expectedLoss) / previousXi;
        return gsl_ran_gaussian_pdf(scaledLoss, sigma);
    }

    /*if (kappa > 0.01f) // Lognormal or convolution approximation
    {
        const float convolutionKappa = std::min(0.3f, kappa);

        const float vavilovMean = PhysicalConstants::m_eulerConstant - 1.f - std::log(vavilovMean) - beta2;
        const float vavilovVariance = (1.f - beta2 / 2.f) / convolutionKappa;
        const float vavilovSkew = (0.5f - beta2 / 3.f) / (convolutionKappa * convolutionKappa * std::pow(vavilovVariance, 1.5f));

        const float scaledMoment2 = vavilovVariance / (vavilovMean * vavilovMean);
        const float scaledMoment3 = vavilovSkew / (vavilovMean * vavilovMean * vavilovMean);
        
        const float paramC = scaledMoment3 / std::pow(scaledMoment2, 1.5f);

        const float s2 = scaledMoment3 * scaledMoment3;
        const float v3 = scaledMoment2 * scaledMoment2 * scaledMoment2;
        const float v6 = v3 * v3;
        const float v9 = v6 * v3;

        const float paramUArg = std::pow(0.5 * (std::sqrt(s2 * v6 + 4.f * v9) + scaledMoment3 * v3), 1.f / 3.f) / std::sqrt(v3);
        const float paramU = paramUArg - 1.f / paramUArg;

        const float lognormalMean = 0.5f * std::log(scaledMoment2 / (paramU * paramU * (1.f + paramU * paramU)));
        const float lognormalVariance = std::log(1.f + paramU * paramU);
        const float lognormalShift = 1.f - std::sqrt(scaledMoment2) / paramU;

        const float deviationFactor = gsl_ran_lognormal(m_pGenerator, lognormalMean, std::sqrt(lognormalVariance)) + lognormalShift;

        if (kappa > 0.3f) // Lognormal approximation
            return meanEnergyLoss * deviationFactor;

        // Convolution approximation
        const float poissonMean = 0.3f - kappa * (1.f + beta2 * std::log(0.3f / kappa));
        const unsigned int numCollisions = gsl_ran_poisson(m_pGenerator, poissonMean);
        float poissonLosses = 0.f;

        for (unsigned int i = 0; i < numCollisions; ++i)
            poissonLosses += this->SampleFromRutherfordDistribution(xi / 0.3f, maxEnergyTransfer, beta2);

        return meanEnergyLoss * deviationFactor + poissonLosses;
    }*/

    // Landau approximation
    const float scaledLoss = (-deltaE - expectedLoss) / previousXi - 1.f - previousBeta2 + PhysicalConstants::m_eulerConstant;
    
    // std::cout << "For particle with mass                 " << spParticle->Mass() << std::endl;
    // std::cout << "For particle with beta                 " << spParticle->Beta() << std::endl;
    // std::cout << "Observation corresponds to dE of       " << deltaE << std::endl;
    // std::cout << "This corresponds to a lambda value of  " << scaledLoss << std::endl;
    // std::cout << "This corresponds to p(lambda) of       " << gsl_ran_landau_pdf(scaledLoss) << std::endl;

    return gsl_ran_landau_pdf(scaledLoss);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Propagator::DensityCorrection(const float beta) const
{
    const float x = std::log10(beta) - 0.5 * std::log10(1.f - beta * beta);

    if (x > m_detector.m_sternheimerX1)
        return 2. * std::log(10.) * x - m_cBar;

    if (x > m_detector.m_sternheimerX0)
        return 2. * std::log(10.) * x - m_cBar + m_detector.m_sternheimerA * std::pow(m_detector.m_sternheimerX1 - x, m_detector.m_sternheimerK);

    return m_detector.m_sternheimerDelta0 == 0.f ? 0.f : m_detector.m_sternheimerDelta0 * std::pow(10.f, 2.f * (x - m_detector.m_sternheimerX0));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Propagator::SampleDeltaE(const float meanEnergyLoss, const float kappa, const float beta2, const float xi, const float maxEnergyTransfer) const
{
    if (kappa > 10.f) // Gaussian approximation
    {
        const float sigma = std::sqrt((1.f - 0.5f * beta2) / kappa);
        return meanEnergyLoss + xi * gsl_ran_gaussian(m_pGenerator, sigma);
    }

   /* if (kappa > 0.01f) // Lognormal or convolution approximation
    {
        const float convolutionKappa = std::min(0.3f, kappa);

        const float vavilovMean = PhysicalConstants::m_eulerConstant - 1.f - std::log(vavilovMean) - beta2;
        const float vavilovVariance = (1.f - beta2 / 2.f) / convolutionKappa;
        const float vavilovSkew = (0.5f - beta2 / 3.f) / (convolutionKappa * convolutionKappa * std::pow(vavilovVariance, 1.5f));

        const float scaledMoment2 = vavilovVariance / (vavilovMean * vavilovMean);
        const float scaledMoment3 = vavilovSkew / (vavilovMean * vavilovMean * vavilovMean);
        
        const float paramC = scaledMoment3 / std::pow(scaledMoment2, 1.5f);

        const float s2 = scaledMoment3 * scaledMoment3;
        const float v3 = scaledMoment2 * scaledMoment2 * scaledMoment2;
        const float v6 = v3 * v3;
        const float v9 = v6 * v3;

        const float paramUArg = std::pow(0.5 * (std::sqrt(s2 * v6 + 4.f * v9) + scaledMoment3 * v3), 1.f / 3.f) / std::sqrt(v3);
        const float paramU = paramUArg - 1.f / paramUArg;

        const float lognormalMean = 0.5f * std::log(scaledMoment2 / (paramU * paramU * (1.f + paramU * paramU)));
        const float lognormalVariance = std::log(1.f + paramU * paramU);
        const float lognormalShift = 1.f - std::sqrt(scaledMoment2) / paramU;

        const float deviationFactor = gsl_ran_lognormal(m_pGenerator, lognormalMean, std::sqrt(lognormalVariance)) + lognormalShift;

        if (kappa > 0.3f) // Lognormal approximation
            return meanEnergyLoss * deviationFactor;

        // Convolution approximation
        const float poissonMean = 0.3f - kappa * (1.f + beta2 * std::log(0.3f / kappa));
        const unsigned int numCollisions = gsl_ran_poisson(m_pGenerator, poissonMean);
        float poissonLosses = 0.f;

        for (unsigned int i = 0; i < numCollisions; ++i)
            poissonLosses += this->SampleFromRutherfordDistribution(xi / 0.3f, maxEnergyTransfer, beta2);

        return meanEnergyLoss * deviationFactor + poissonLosses;
    }*/

    // Landau approximation
    const float lambda = gsl_ran_landau(m_pGenerator);
    const float deltaE = meanEnergyLoss + xi * (lambda + 1.f + beta2 - PhysicalConstants::m_eulerConstant);
   // std::cout << "LAMBDA =  " << lambda << " DELTAE = " << deltaE << std::endl;
    return deltaE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Propagator::SampleFromRutherfordDistribution(const float eMin, const float eMax, const float beta2) const
{
    // Rejection sampling
    const float normConstant = eMin * eMax / (eMax - eMin + beta2 * eMin * std::log(eMin/eMax));
    const float rejectionM = normConstant * (eMax - eMin) * (eMax - beta2 * eMin) / (eMax * eMin * eMin);

    while (true)
    {
        const float sample = gsl_ran_flat(m_pGenerator, eMin, eMax);
        const float rutherfordValue = normConstant * (1.f / (sample * sample) - beta2 / (eMax * sample));

        const float acceptanceProbability = (eMax - eMin) * rutherfordValue / rejectionM;

        if (acceptanceProbability > 1.f || acceptanceProbability < 0.f)
            throw std::runtime_error("Acceptance probability was not within range");

        if (gsl_ran_flat(m_pGenerator, 0.f, 1.f) < acceptanceProbability)
            return sample;
    }

    assert(false && "Unreachable");
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Propagator::GetMaxEnergyTransfer(const float beta, const float mass) const
{
    const float beta2 = beta * beta;
    const float gamma = 1. / std::sqrt(1 - beta2);

    const float numerator   = 2.f * PhysicalConstants::m_electronMass * beta2 * gamma * gamma;
    const float denominator = 1.f + 2.f * gamma * PhysicalConstants::m_electronMass / mass +
                              (PhysicalConstants::m_electronMass / mass) * (PhysicalConstants::m_electronMass / mass);

    return numerator / denominator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Propagator::GetExpectedEnergyLoss(const float beta, const float beta2, const float xi) const 
{
    const float lnArg =
        2.f * PhysicalConstants::m_electronMass * beta2 * xi /
        (m_detector.m_density * m_detector.m_thickness * (1.f - beta2) * m_detector.m_avgIonizationEnergy * m_detector.m_avgIonizationEnergy);

    return -xi * (std::log(lnArg) - 2.f * beta2 - this->DensityCorrection(beta));
}

} // namespace bf

