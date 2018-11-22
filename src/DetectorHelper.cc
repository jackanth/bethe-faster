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

    detector.m_density               = 1.40f;     // g/cm^3
    detector.m_thickness             = 0.03f;     // cm
    detector.m_avgIonizationEnergy   = 188.0f;    // eV
    detector.m_atomicNumber          = 18U;       // (no units)
    detector.m_atomicMass            = 39.95f;    // g/mol
    detector.m_plasmaEnergy          = 22.89f;    // eV
    detector.m_sternheimerX0         = 0.2000f;   // (no units)
    detector.m_sternheimerX1         = 3.0000f;   // (no units)
    detector.m_sternheimerK          = 3.0000f;   // (no units)
    detector.m_sternheimerA          = 0.19559f;  // (no units)
    detector.m_sternheimerDelta0     = 0.00f;     // (no units)
    detector.m_birksA                = 0.800f;    // (no units)
    detector.m_birksK                = 0.0468f;   // g kV / (MeV cm^3)
    detector.m_electricFieldStrength = 0.273f;    // keV/cm
    detector.m_chargeCalibration     = 5.076e-3f; // ADC/e
    detector.m_ionizationEnergy      = 23.6e-6f;  // MeV/electron

    return detector;
}

} // namespace bf
