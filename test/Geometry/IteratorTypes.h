#ifndef TEST_GEOMETRY_ITERATORTYPES_H
#define TEST_GEOMETRY_ITERATORTYPES_H

// LArSoft libraries
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/details/element_iterators.h"
#include "larcorealg/Geometry/details/id_iterators.h"
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

namespace geo {
  /**
   * @brief Forward iterator browsing all element IDs in the detector
   *
   * Prefer asking GeometryCore object for iterators rather than constructing them anew.
   * Stand-alone example (not recommended):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * geo::cryostat_id_iterator const cbegin{geom, geom->begin_wire()};
   * geo::cryostat_id_iterator const cend{geom, geom->end_wire()};
   * for (auto iCryostat = cbegin; iCryostat != cend; ++iCryostat) {
   *   geo::CryostatID const& cid = *iCryostat;
   *   std::cout << "We are at: " << cid << std::endl;
   *   // ...
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  using cryostat_id_iterator = details::id_iterator<CryostatID, details::GeometryIterationPolicy>;
  using TPC_id_iterator = details::id_iterator<TPCID, details::GeometryIterationPolicy>;
  using plane_id_iterator = details::id_iterator<PlaneID, details::ReadoutIterationPolicy>;
  using wire_id_iterator = details::id_iterator<WireID, details::ReadoutIterationPolicy>;

  using TPCset_id_iterator =
    details::id_iterator<readout::TPCsetID, details::ReadoutIterationPolicy>;
  using ROP_id_iterator = details::id_iterator<readout::ROPID, details::ReadoutIterationPolicy>;

  using cryostat_id_sentinel = details::id_sentinel<CryostatID>;
  using TPC_id_sentinel = details::id_sentinel<TPCID>;
  using plane_id_sentinel = details::id_sentinel<PlaneID>;
  using wire_id_sentinel = details::id_sentinel<WireID>;

  using TPCset_id_sentinel = details::id_sentinel<readout::TPCsetID>;
  using ROP_id_sentinel = details::id_sentinel<readout::ROPID>;

  /**
   * @brief Forward iterator browsing all elements in the detector
   *
   * The comments from the ID iterators above are valid here as well.
   * This object has a different dereferencing operator that obtains
   * the elements directly, or throws on failure.
   */
  using cryostat_iterator =
    details::element_iterator_for<GeometryCore, CryostatGeo, details::GeometryIterationPolicy>;
  using TPC_iterator =
    details::element_iterator_for<GeometryCore, TPCGeo, details::GeometryIterationPolicy>;
  using plane_iterator =
    details::element_iterator_for<WireReadoutGeom, PlaneGeo, details::ReadoutIterationPolicy>;
  using wire_iterator =
    details::element_iterator_for<WireReadoutGeom, WireGeo, details::ReadoutIterationPolicy>;

  using cryostat_sentinel = details::element_sentinel_for<CryostatGeo>;
  using TPC_sentinel = details::element_sentinel_for<TPCGeo>;
  using plane_sentinel = details::element_sentinel_for<PlaneGeo>;
  using wire_sentinel = details::element_sentinel_for<WireGeo>;
}

#endif // LARCOREALG_GEOMETRY_ITERATORTYPES_H
