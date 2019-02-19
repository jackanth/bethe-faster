/**
 *  @file   bethe-faster/include/QuickPidHelper.h
 *
 *  @brief  Header file for the quick PID helper class.
 *
 *  $Log: $
 */

#ifndef BF_QUICK_PID_HELPER_H
#define BF_QUICK_PID_HELPER_H 1

#include "HitCharge.h"
#include <vector>
#include <tuple>

namespace bf
{

/**
 *  @brief  QuickPidHelper class
 */
class QuickPidHelper
{
public:
    /**
     *  @brief  Deleted constructor
     */
    QuickPidHelper() = delete;

    /**
     *  @brief  Calculate the Bragg gradient and intercept parameters
     *
     *  @param hitChargeVector the vector of hit charge objects
     *
     *  @return the gradient and intercept
     */
    static bool CalculateBraggGradient(const std::vector<HitCharge> &hitChargeVector, double &gradient, double &intercept);

    /**
     *  @brief  Calculate the median of some numbers
     *
     *  @param  vec the vector of numbers
     *
     *  @return the median
     */
    static double CalculateMedian(std::vector<double> vec);
};

} // namespace bf

#endif // #ifndef BF_QUICK_PID_HELPER_H
