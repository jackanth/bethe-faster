/**
 *  @file   bethe-faster/include/HitCharge.h
 *
 *  @brief  Header file for the hit charge class.
 *
 *  $Log: $
 */

#ifndef BF_HIT_CHARGE_H
#define BF_HIT_CHARGE_H 1

#include <stdexcept>

namespace bf
{
/**
 *  @brief  HitCharge class
 */
class HitCharge
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  coordinate 3D coordinate (cm)
     *  @param  energyLossRate the dE/dx value (MeV/cm)
     *  @param  extent the 3D extent (cm)
     */
    HitCharge(const double coordinate, const double energyLossRate, const double extent) noexcept;

    /**
     *  @brief  Get the 3D coordinate (cm)
     *
     *  @return the coordinate
     */
    double Coordinate() const noexcept;

    /**
     *  @brief  Get the dE/dx value (MeV/cm)
     *
     *  @return the dE/dx value
     */
    double EnergyLossRate() const noexcept;

    /**
     *  @brief  Get the 3D extent (cm)
     *
     *  @return the extent
     */
    double Extent() const noexcept;

    /**
     *  @brief  Set the 3D coordinate (cm)
     *
     *  @param coordinate the coordinate
     */
    void Coordinate(const double coordinate) noexcept;

    /**
     *  @brief  Set the dE/dx value (MeV/cm)
     *
     *  @param energyLossRate the dE/dx value
     */
    void EnergyLossRate(const double energyLossRate) noexcept;

    /**
     *  @brief  Set the 3D extent (cm)
     *
     *  @param extent the extent
     */
    void Extent(const double extent) noexcept;

private:
    double m_coordinate;     ///< The 3D coordinate (cm)
    double m_energyLossRate; ///< The energy loss rate (MeV/cm)
    double m_extent;         ///< The 3D extent (cm)
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline HitCharge::HitCharge(const double coordinate, const double energyLossRate, const double extent) noexcept :
    m_coordinate{coordinate},
    m_energyLossRate{energyLossRate},
    m_extent{extent}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double HitCharge::Coordinate() const noexcept
{
    return m_coordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double HitCharge::EnergyLossRate() const noexcept
{
    return m_energyLossRate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double HitCharge::Extent() const noexcept
{
    return m_extent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HitCharge::Coordinate(const double coordinate) noexcept
{
    m_coordinate = coordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HitCharge::EnergyLossRate(const double energyLossRate) noexcept
{
    m_energyLossRate = energyLossRate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HitCharge::Extent(const double extent) noexcept
{
    m_extent = extent;
}

} // namespace bf

#endif // #ifndef BF_HIT_CHARGE_H
