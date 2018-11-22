/**
 *  @file   bethe-faster/include/DetectorHelper.h
 *
 *  @brief  Header file for the detector helper class.
 *
 *  $Log: $
 */
#ifndef BF_DETECTOR_HELPER_H
#define BF_DETECTOR_HELPER_H 1

#include "Detector.h"

namespace bf
{

/**
 *  @brief  DetectorHelper class
 */
class DetectorHelper
{
public:
    /**
     *  @brief  Deleted constructor
     */
    DetectorHelper() = delete;

    /**
     *  @brief  Get a MicroBooNE-like detector
     * 
     *  @return the MicroBooNE detector
     */
    static Detector GetMicroBooNEDetector() noexcept;
};

} // namespace bf

#endif // #ifndef BF_DETECTOR_HELPER_H
