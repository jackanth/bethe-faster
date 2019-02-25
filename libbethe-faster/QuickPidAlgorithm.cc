/**
 *  @file   bethe-faster/src/QuickPidAlgorithm.cc
 *
 *  @brief  Implementation of the QuickPidAlgorithm class.
 *
 *  $Log: $
 */

#include "QuickPidAlgorithm.h"
#include "PhysicalConstants.h"
#include "PlotHelper.h"

#include <iostream>

namespace bf
{

QuickPidAlgorithm::QuickPidAlgorithm(Detector detector) :
    m_detector{std::move_if_noexcept(detector)},
    m_xiPrime{0.},
    m_chiPartial{0.}
{
    if (2. * m_detector.m_atomicMass < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize quick PID algorithm because the atomic mass was too small");

    m_xiPrime = PhysicalConstants::m_kCoefficient * static_cast<double>(m_detector.m_atomicNumber) * m_detector.m_density /
                       (2. * m_detector.m_atomicMass);

    const double denominator = m_detector.m_avgIonizationEnergy * m_detector.m_avgIonizationEnergy * 1.e-12;

    if (denominator < std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Could not initialize quick PID algorithm because the average ionization energy was too small");

    m_chiPartial = std::log(2. * PhysicalConstants::m_electronMass * m_xiPrime / denominator) + PhysicalConstants::m_vavilovJ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool QuickPidAlgorithm::CalculateSecondOrderBraggGradient(std::vector<HitCharge> hitChargeVector, double &gradient, double &intercept) const
{
    if (hitChargeVector.empty())
        throw std::runtime_error{"No data points were provided"};
    
    // Calculate the maximum coordinate
    double maxCoordinate = 0.;

    for (const auto &hitCharge : hitChargeVector)
    {
        if (hitCharge.Coordinate() > maxCoordinate)
            maxCoordinate = hitCharge.Coordinate();
    }

    // Create the plot vector.
    std::vector<std::pair<double, double>> plotVector;

    for (const auto &hitCharge : hitChargeVector)
    {
        try 
        {
            const double xValue = maxCoordinate - hitCharge.Coordinate();

            const double chi = this->Chi(hitCharge.Extent());
            const double dEdx = -hitCharge.EnergyLossRate();
        
            const double tPrime = this->CalculateTPrime(dEdx, hitCharge.Extent());
            const double logTerm = std::log(1. + 2. * tPrime * tPrime / chi) / 3.;
            const double arctanFactor = std::sqrt(chi / 2.);
            const double arctanTerm = arctanFactor * std::atan(tPrime / arctanFactor);
            const double yValue = logTerm + arctanTerm - tPrime;
            
            plotVector.emplace_back(xValue, yValue);
        }
        catch (...)
        {
            continue;
        }
    }

    const std::size_t   numHitCharges = plotVector.size();
    std::vector<double> slopeMedianVector, interceptMedianVector;

    // Perform repeated median regression
    for (std::size_t i = 0UL; i < numHitCharges; ++i)
    {
        const auto [xi, yi] = plotVector.at(i);
        std::vector<double> slopeVector, interceptVector;

        for (std::size_t j = 0UL; j < numHitCharges; ++j)
        {
            if (i == j)
                continue;

            const auto [xj, yj] = plotVector.at(j);
          

            if (std::fabs(xj - xi) < std::numeric_limits<double>::epsilon())
                continue;

            const double pointIntercept = (xj * yi - xi * yj) / (xj - xi);
            const double pointSlope     = (yj - yi) / (xj - xi);

            slopeVector.push_back(pointSlope);
            interceptVector.push_back(pointIntercept);
        }

        if (slopeVector.empty() || interceptVector.empty())
            continue;

        slopeMedianVector.push_back(QuickPidHelper::CalculateMedian(slopeVector));
        interceptMedianVector.push_back(QuickPidHelper::CalculateMedian(interceptVector));
    }

    if (slopeMedianVector.empty() || interceptMedianVector.empty())
        return false;

    gradient  = QuickPidHelper::CalculateMedian(slopeMedianVector);
    intercept = QuickPidHelper::CalculateMedian(interceptMedianVector);
    return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

TGraph QuickPidAlgorithm::GetSecondOrderBraggGradientGraph(const std::vector<HitCharge> &hitChargeVector) const
{
    if (hitChargeVector.empty())
        throw std::runtime_error{"Cannot plot empty hit charge distribution"};

    double maxCoordinate = 0.;

    for (const auto &hitCharge : hitChargeVector)
    {
        if (hitCharge.Coordinate() > maxCoordinate)
            maxCoordinate = hitCharge.Coordinate();
    }

    std::vector<double> xVector, yVector;

    for (const auto &hitCharge : hitChargeVector)
    {
        try 
        {
            const double chi = this->Chi(hitCharge.Extent());
            const double dEdx = -hitCharge.EnergyLossRate();    
            const double tPrime = this->CalculateTPrime(dEdx, hitCharge.Extent());
            const double logTerm = std::log(1. + 2. * tPrime * tPrime / chi) / 3.;
            const double arctanFactor = std::sqrt(chi / 2.);
            const double arctanTerm = arctanFactor * std::atan(tPrime / arctanFactor);
            const double yValue = logTerm + arctanTerm - tPrime;

            xVector.push_back(maxCoordinate - hitCharge.Coordinate());
            yVector.push_back(yValue);
        }
        catch (...)
        {
            continue;
        }
    }

    return TGraph{static_cast<Int_t>(xVector.size()), xVector.data(), yVector.data()};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

double QuickPidAlgorithm::CalculateScaledKineticEnergy(const double residualRange, const double mass, const double deltaX) const
{
    const double chi       = this->Chi(deltaX);
    const double prefactor = mass / (2. * m_xiPrime);

    // Newton's method
    std::vector<double> initialEstimates = {0.2, 0.05, 0.01, 0.5, 0.75, 1.0};

    for (double scaledKineticEnergy : initialEstimates)
    {
        bool failed = false;
        std::uint32_t loopCounter = 1U;
        double func = this->NewtonMethodFunc(scaledKineticEnergy, residualRange, prefactor, chi);

        while (std::abs(func) > 0.000001)
        {
            double funcDeriv = this->NewtonMethodFuncDeriv(scaledKineticEnergy, prefactor, chi);
            scaledKineticEnergy -= func / funcDeriv;
            func = this->NewtonMethodFunc(scaledKineticEnergy, residualRange, prefactor, chi);

            if (loopCounter++ >= 10'000) 
            {
                failed = true;
                break;
            }
        }

        if (failed)
            continue;

        if (scaledKineticEnergy > -std::numeric_limits<double>::epsilon())
            return scaledKineticEnergy;
    }

    throw std::runtime_error{"Newton's method failed to find the positive root"};
}

//------------------------------------------------------------------------------------------------------------------------------------------

double QuickPidAlgorithm::NewtonMethodFunc(const double scaledKineticEnergy, const double residualRange, const double prefactor, const double chi) const
{
    const double logTerm = std::log(1. + 2. * scaledKineticEnergy * scaledKineticEnergy / chi) / 3.;

    const double arctanPrefactor = std::sqrt(0.5 * chi);
    const double arctanTerm      = arctanPrefactor * std::atan(scaledKineticEnergy / arctanPrefactor);

    return 3. * prefactor * (logTerm + arctanTerm - scaledKineticEnergy) - residualRange;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

double QuickPidAlgorithm::CalculateTPrime(const double dEdx, const double deltaX) const
{
    const double chi = this->Chi(deltaX);

    const auto newtonFunc = [&](const double tPrime) {
        const double bracketedTerm = 2. * std::log(1. + tPrime) + chi;

        if (tPrime * (tPrime + 2.) < std::numeric_limits<double>::epsilon())
            throw std::runtime_error{"Could not calculate T' as denominator was too small"};

        const double prefactor = (tPrime + 1.) * (tPrime + 1.) / (tPrime * (tPrime + 2.));
        return dEdx + m_xiPrime * (prefactor * bracketedTerm - 1.);
    };

    const auto newtonFuncDeriv = [&](const double tPrime) {
        const double denominator = tPrime * tPrime * (tPrime + 2.) * (tPrime + 2.);

        if (denominator < std::numeric_limits<double>::epsilon())
            throw std::runtime_error{"Could not calculate T' as denominator was too small"};

        const double bracketedTerm = chi - tPrime * (tPrime + 2.) + 2. * std::log(1. + tPrime);
        return - 2. * m_xiPrime * (tPrime + 1.) * bracketedTerm / denominator;
    };

    // Newton's method
    std::vector<double> initialEstimates = {0.2, 0.05, 0.01, 0.5, 1.0};

    for (double scaledKineticEnergy : initialEstimates) {
        double func = newtonFunc(scaledKineticEnergy);
        bool failed = false;
        std::uint32_t loopCounter = 1U;

        while (std::abs(func) > 0.001)
        {
            try
            {
                double funcDeriv = newtonFuncDeriv(scaledKineticEnergy);
                scaledKineticEnergy -= func / funcDeriv;
                func = newtonFunc(scaledKineticEnergy);
            }
            catch (...)
            {
                failed = true;
                break;
            }

            if (loopCounter++ >= 10'000) 
            {
                failed = true;
                break;
            }
        }

        if (failed)
            continue;

        if (scaledKineticEnergy > -std::numeric_limits<double>::epsilon())
            return scaledKineticEnergy;
    }

    throw std::runtime_error{"Newton's method failed to find a positive root"};
}

} // namespace bf
