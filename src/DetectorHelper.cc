/**
 *  @file   bethe-faster/src/DetectorHelper.cc
 *
 *  @brief  Implementation of the detector helper class.
 *
 *  $Log: $
 */

#include "DetectorHelper.h"

namespace bf
{

Detector DetectorHelper::GetMicroBooNEDetector() noexcept
{
    Detector detector;

    detector.m_density               = 1.40;     // g/cm^3
    detector.m_avgIonizationEnergy   = 188.0;    // eV
    detector.m_atomicNumber          = 18U;      // (no units)
    detector.m_atomicMass            = 39.95;    // g/mol
    detector.m_plasmaEnergy          = 22.89;    // eV
    detector.m_sternheimerX0         = 0.2000;   // (no units)
    detector.m_sternheimerX1         = 3.0000;   // (no units)
    detector.m_sternheimerK          = 3.0000;   // (no units)
    detector.m_sternheimerA          = 0.19559;  // (no units)
    detector.m_sternheimerDelta0     = 0.00;     // (no units)

    return detector;
}

} // namespace bf
