/**
 *  @file   bethe-faster/include/Detector.h
 *
 *  @brief  Header file for the detector struct.
 *
 *  $Log: $
 */
#ifndef BF_DETECTOR_H
#define BF_DETECTOR_H 1

#include <stdexcept>
#include <limits>

namespace bf
{

/**
 *  @brief  Detector struct
 */
struct Detector
{
    /**
     *  @brief  Constructor
     */
    Detector() noexcept;

    float        m_density;               ///< The density (g/cm^3)
    float        m_thickness;             ///< The thickness (cm)
    float        m_avgIonizationEnergy;   ///< The average ionization energy (eV)
    unsigned int m_atomicNumber;          ///< The atomic number
    float        m_atomicMass;            ///< The atomic mass (g/mol)
    float        m_plasmaEnergy;          ///< The plasma energy (eV)
    float        m_sternheimerX0;         ///< The Sternheimer x0 parameter
    float        m_sternheimerX1;         ///< The Sternheimer x1 parameter
    float        m_sternheimerK;          ///< The Sternheimer k parameter
    float        m_sternheimerA;          ///< The Sternheimer a parameter
    float        m_sternheimerDelta0;     ///< The Sternheimer delta0 parameter
    float        m_birksA;                ///< The Birks A parameter
    float        m_birksK;                ///< The Birks k parameter (g kV / (MeV cm^3))
    float        m_electricFieldStrength; ///< The electric field strength (keV/cm)
    float        m_chargeCalibration;     ///< The charge calibration constant (ADC/e)
    float        m_ionizationEnergy;      ///< The ionization energy (MeV/electron)
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline Detector::Detector() noexcept :
    m_density(0.f),
    m_thickness(0.f),
    m_avgIonizationEnergy(0.f),
    m_atomicNumber(0U),
    m_atomicMass(0.f),
    m_plasmaEnergy(0.f),
    m_sternheimerX0(0.f),
    m_sternheimerX1(0.f),
    m_sternheimerK(0.f),
    m_sternheimerA(0.f),
    m_sternheimerDelta0(0.f),
    m_birksA(0.f),
    m_birksK(0.f),
    m_electricFieldStrength(0.f),
    m_chargeCalibration(0.f),
    m_ionizationEnergy(0.f)
{
}

} // namespace bf

#endif // #ifndef BF_DETECTOR_H
