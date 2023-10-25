/**
 * @file   larcorealg/Geometry/ReadoutIDmapper.h
 * @brief  Mapping between geometry/readout ID and flat index.
 * @ingroup Geometry
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_GEOMETRY_READOUTIDMAPPER_H
#define LARCOREALG_GEOMETRY_READOUTIDMAPPER_H

// LArSoft libraries
#include "larcorealg/Geometry/GeometryIDmapper.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

// C/C++ standard libraries
#include <algorithm> // std::fill(), std::for_each()
#include <array>
#include <cassert>
#include <initializer_list>
#include <stdexcept> // std::out_of_range
#include <string>
#include <utility> // std::forward()
#include <vector>

namespace readout {

  template <typename Index = std::size_t>
  class TPCsetIDmapper;

  template <typename Index = std::size_t>
  class ROPIDmapper;

} // namespace readout

// --- BEGIN Readout ID mappers ------------------------------------------------
/// @name Readout ID mappers
/// @ingroup Geometry
/// @{

/** ****************************************************************************
 * @brief Mapping for TPC set identifiers.
 * @tparam Index (default: `std::size_t`) type of flat index
 * @see `geo::GeoIDmapper`
 *
 * A customized version of `geo::GeoIDmapper` offering TPC set ID-specific
 * interface.
 */
template <typename Index /* = std::size_t */>
class readout::TPCsetIDmapper : public geo::GeoIDmapper<readout::TPCsetID, Index> {

  /// Base class.
  using BaseMapper_t = geo::GeoIDmapper<readout::TPCsetID, Index>;

public:
  // import types
  using ID_t = typename BaseMapper_t::ID_t;
  using index_type = typename BaseMapper_t::index_type;

  // import all constructors from `geo::GeoIDmapper`
  using BaseMapper_t::BaseMapper_t;

  /**
   * @brief Prepares the mapping with the specified sizes.
   * @param nCryo number of cryostats
   * @param nTPCsets number of TPCsets per cryostat
   *
   * The mapping is sized to map `nCryo` cryostats, each with `nTPCsets` TPC
   * sets.
   */
  TPCsetIDmapper(unsigned int nCryo, unsigned int nTPCsets) : BaseMapper_t({nCryo, nTPCsets}) {}

  // --- BEGIN Mapping modification --------------------------------------------
  /// @name Mapping modification
  /// @{

  using BaseMapper_t::resize;

  /**
   * @brief Prepares the mapping for the specified sizes.
   * @param nCryo number of cryostats
   * @param nTPCsets number of TPC sets
   * @see `resizeAs()`
   *
   * The mapping is sized to map `nCryo` cryostats, each with `nTPCsets` TPC
   * sets.
   */
  void resize(unsigned int nCryo, unsigned int nTPCsets)
  {
    BaseMapper_t::resize({nCryo, nTPCsets});
  }

  /// @}
  // --- END Mapping modification ----------------------------------------------

  // --- BEGIN Mapping status query --------------------------------------------
  /// @name Mapping status query
  /// @{

  /// Returns whether this mapping covers the specified cryostat.
  bool hasCryostat(geo::CryostatID const& cryoid) const { return BaseMapper_t::hasElement(cryoid); }

  /// Returns whether this mapping covers the specified TPC set.
  bool hasTPCset(readout::TPCsetID const& tpcsetid) const
  {
    return BaseMapper_t::hasElement(tpcsetid);
  }

  /// @}
  // --- END Mapping status query ----------------------------------------------

}; // readout::TPCsetIDmapper<>

/** ****************************************************************************
 * @brief Mapping for readout plane identifiers.
 * @tparam Index (default: `std::size_t`) type of flat index
 * @see `geo::GeoIDmapper`
 *
 * A customized version of `geo::GeoIDmapper` offering
 * readout plane ID-specific interface.
 */
template <typename Index /* = std::size_t */>
class readout::ROPIDmapper : public geo::GeoIDmapper<readout::ROPID, Index> {

  /// Base class.
  using BaseMapper_t = geo::GeoIDmapper<readout::ROPID, Index>;

public:
  // import types
  using ID_t = typename BaseMapper_t::ID_t;
  using index_type = typename BaseMapper_t::index_type;

  // import all constructors from `geo::GeoIDmapper`
  using BaseMapper_t::BaseMapper_t;

  /**
   * @brief Prepares the mapping with the specified sizes.
   * @param nCryo number of cryostats
   * @param nTPCsets number of TPC sets per cryostat
   * @param nROPs number of readout planes per TPC set
   *
   * The mapping is sized to map `nCryo` cryostats, each with `nTPCsets` TPC
   * sets, each one with `nROPs` readout planes.
   */
  ROPIDmapper(unsigned int nCryo, unsigned int nTPCsets, unsigned int nROPs)
    : BaseMapper_t({nCryo, nTPCsets, nROPs})
  {}

  // --- BEGIN Mapping modification --------------------------------------------
  /// @name Mapping modification
  /// @{

  using BaseMapper_t::resize;

  /**
   * @brief Prepares the mapping for the specified sizes.
   * @param nCryo number of cryostats
   * @param nTPCsets number of TPC sets per cryostat
   * @param nROPs number of readout planes per TPC set
   * @see `resizeAs()`
   *
   * The mapping is sized to map `nCryo` cryostats, each with `nTPCsets` TPC
   * sets, each one with `nROPs` readout planes.
   */
  void resize(unsigned int nCryo, unsigned int nTPCsets, unsigned int nROPs)
  {
    BaseMapper_t::resize({nCryo, nTPCsets, nROPs});
  }

  /// @}
  // --- END Mapping modification ----------------------------------------------

  // --- BEGIN Mapping status query --------------------------------------------
  /// @name Mapping status query
  /// @{

  /// Returns whether this mapping covers the specified cryostat.
  bool hasCryostat(geo::CryostatID const& cryoid) const { return BaseMapper_t::hasElement(cryoid); }

  /// Returns whether this mapping covers the specified TPC set.
  bool hasTPC(readout::TPCsetID const& tpcsetid) const
  {
    return BaseMapper_t::hasElement(tpcsetid);
  }

  /// Returns whether this mapping covers the specified readout plane.
  bool hasPlane(readout::ROPID const& ropid) const { return BaseMapper_t::hasElement(ropid); }

  /// @}
  // --- END Mapping status query ----------------------------------------------

}; // readout::ROPIDmapper<>

/// @}
// --- END Readout ID mappers --------------------------------------------------
//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_READOUTIDMAPPER_H
