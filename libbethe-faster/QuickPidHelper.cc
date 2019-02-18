/**
 *  @file   bethe-faster/src/QuickPidHelper.cc
 *
 *  @brief  Implementation of the quick PID helper class.
 *
 *  $Log: $
 */

#include "QuickPidHelper.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace bf
{

std::tuple<double, double> QuickPidHelper::CalculateBraggGradient(const std::vector<HitCharge> &hitChargeVector) 
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

    const std::size_t   numHitCharges = hitChargeVector.size();
    std::vector<double> slopeMedianVector, interceptMedianVector;

    // Perform repeated median regression
    for (std::size_t i = 0UL; i < numHitCharges; ++i)
    {
        const auto hitChargeI = hitChargeVector.at(i);

        const double xiDenominator = std::sqrt(maxCoordinate - hitChargeI.Coordinate());

        if (xiDenominator < std::numeric_limits<double>::epsilon())
            continue;

        const auto xi = 1. / xiDenominator;
        const auto yi = hitChargeI.EnergyLossRate();

        std::vector<double> slopeVector, interceptVector;

        for (std::size_t j = 0UL; j < numHitCharges; ++j)
        {
            if (i == j)
                continue;

            const auto hitChargeJ = hitChargeVector.at(j);
            const double xjDenominator = std::sqrt(maxCoordinate - hitChargeJ.Coordinate());

            if (xjDenominator < std::numeric_limits<double>::epsilon())
                continue;

            const auto xj = 1. / xjDenominator;
            const auto yj = hitChargeJ.EnergyLossRate();

            if (std::fabs(xj - xi) < std::numeric_limits<double>::epsilon())
                continue;

            const double intercept = (xj * yi - xi * yj) / (xj - xi);
            const double slope     = (yj - yi) / (xj - xi);

            slopeVector.push_back(slope);
            interceptVector.push_back(intercept);
        }

        if (slopeVector.empty() || interceptVector.empty())
            continue;

        slopeMedianVector.push_back(QuickPidHelper::CalculateMedian(slopeVector));
        interceptMedianVector.push_back(QuickPidHelper::CalculateMedian(interceptVector));
    }

    if (slopeMedianVector.empty() || interceptMedianVector.empty())
        throw std::runtime_error{"Not enough good hits to calculate Bragg gradient"};

    const double slope     = QuickPidHelper::CalculateMedian(slopeMedianVector);
    const double intercept = QuickPidHelper::CalculateMedian(interceptMedianVector);

    return {slope, intercept};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

double QuickPidHelper::CalculateMedian(std::vector<double> vec)
{
    const std::size_t numEntries = vec.size();

    if (numEntries == 0UL)
        throw std::runtime_error{"Cannot calculate median of zero items"};

    std::sort(vec.begin(), vec.end());

    if (numEntries % 2 == 0)
        return (vec.at(numEntries / 2 - 1) + vec.at(numEntries / 2)) / 2.;

    return vec.at(numEntries / 2);
}

} // namespace bf
