/**
 *  @file   bethe-faster/include/PhysicalConstants.h
 *
 *  @brief  Header file for the physical constnats struct.
 *
 *  $Log: $
 */

#ifndef BF_PHYSICAL_CONSTANTS_H
#define BF_PHYSICAL_CONSTANTS_H 1

namespace bf
{

/**
 *  @brief  PhysicalConstants struct
 */
struct PhysicalConstants
{
    /**
     *  @brief  Deleted constructor
     */
    PhysicalConstants() = delete;

    static constexpr double m_kCoefficient    = 0.307075;      ///< dE/dx coefficient (value from PDG) (MeV mol-1 cm^2)
    static constexpr double m_eulerConstant   = 0.5772156649;  ///< no units
    static constexpr double m_electronMass    = 0.5109989461;  ///< MeV
    static constexpr double m_muonMass        = 105.6583745;   ///< MeV
    static constexpr double m_protonMass      = 938.2720813;   ///< MeV
    static constexpr double m_chargedPionMass = 139.57018;     ///< MeV
    static constexpr double m_chargedKaonMass = 493.677;       ///< MeV
    static constexpr double m_landauMode      = -0.22278;      ///< no units
    static constexpr double m_ln10            = 2.30258509299; ///< no units
};

} // namespace bf

#endif // #ifndef BF_PHYSICAL_CONSTANTS_H
