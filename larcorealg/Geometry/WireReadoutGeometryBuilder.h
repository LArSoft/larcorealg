#ifndef LARCOREALG_GEOMETRY_WIREREADOUTGEOMETRYBUILDER_H
#define LARCOREALG_GEOMETRY_WIREREADOUTGEOMETRYBUILDER_H

// LArSoft libraries
#include "larcorealg/Geometry/GeoNodePath.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/fwd.h"

// C++ standard library
#include <cstddef>
#include <map>
#include <vector>

namespace geo {

  /**
   * @brief Manages the extraction of LArSoft geometry information from ROOT.
   *
   * The general interface only provides abstraction for the high level objects (cryostats
   * and auxiliary detectors).  The implementations can use a finer internal structure to
   * address single subcomponents (e.g. wire planes).
   *
   * Builder objects can be configured via FHiCL parameters.
   *
   * This is an abstract interface.
   *
   *
   * Customization of geometry objects
   * ----------------------------------
   *
   * The builders return collections of LArSoft geometry objects dynamically allocated. In
   * this way, in future it will be possible to easily customize those objects for
   * detector-specific needs.  Note that as of LArSoft `v08_06_00`, no polymorphism is
   * actually implemented.
   *
   */
  class WireReadoutGeometryBuilder {
  public:
    // --- BEGIN Data types ----------------------------------------------------
    /// Identification of a single node in ROOT geometry.
    using Path_t = GeoNodePath;

    // --- END Data types ------------------------------------------------------

    // --- BEGIN Constructors and destructor -----------------------------------
    /// Virtual destructor.
    virtual ~WireReadoutGeometryBuilder() = default;

    // --- END Constructors and destructor -------------------------------------

    // --- BEGIN Plane information ------------------------------------------
    /// @name Cryostat information
    /// @{

    /// Collection of wire-plane information objects.
    using Planes_t = std::map<std::size_t, std::vector<PlaneGeo>>;

    /**
     * @brief Looks for all cryostats under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed planes, keyed by TPC TGeoVolume pointer
     *
     * The planes contain all their inner elements.  The current node itself of the
     * path is also considered as plane candidate, then it is descended into.
     */
    Planes_t extractPlanes(Path_t path) const { return doExtractPlanes(path); }

    /// @}
    // --- END Plane information --------------------------------------------

  private:
    virtual Planes_t doExtractPlanes(Path_t& path) const = 0;

  }; // class WireReadoutGeometryBuilder

} // namespace geo

#endif // LARCOREALG_GEOMETRY_WIREREADOUTGEOMETRYBUILDER_H
