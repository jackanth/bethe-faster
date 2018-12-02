/**
 *  @file   bethe-faster/include/Detector.h
 *
 *  @brief  Header file for the detector struct.
 *
 *  $Log: $
 */

#ifndef BF_DETECTOR_H
#define BF_DETECTOR_H 1

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

    double       m_density;               ///< The density (g/cm^3)
    double       m_avgIonizationEnergy;   ///< The average ionization energy (eV)
    unsigned int m_atomicNumber;          ///< The atomic number
    double       m_atomicMass;            ///< The atomic mass (g/mol)
    double       m_plasmaEnergy;          ///< The plasma energy (eV)
    double       m_sternheimerX0;         ///< The Sternheimer x0 parameter
    double       m_sternheimerX1;         ///< The Sternheimer x1 parameter
    double       m_sternheimerK;          ///< The Sternheimer k parameter
    double       m_sternheimerA;          ///< The Sternheimer a parameter
    double       m_sternheimerDelta0;     ///< The Sternheimer delta0 parameter
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline Detector::Detector() noexcept :
    m_density(0.),
    m_avgIonizationEnergy(0.),
    m_atomicNumber(0U),
    m_atomicMass(0.),
    m_plasmaEnergy(0.),
    m_sternheimerX0(0.),
    m_sternheimerX1(0.),
    m_sternheimerK(0.),
    m_sternheimerA(0.),
    m_sternheimerDelta0(0.)
{
}

} // namespace bf

#endif // #ifndef BF_DETECTOR_H
