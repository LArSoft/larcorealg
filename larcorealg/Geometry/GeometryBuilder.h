/**
 * @file   larcorealg/Geometry/GeometryBuilder.h
 * @brief  Interface for geometry extractor classes.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYBUILDER_H
#define LARCOREALG_GEOMETRY_GEOMETRYBUILDER_H

// LArSoft libraries
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeoNodePath.h"

// ROOT libraries
#include "TGeoNode.h"

// C++ standard library
#include <algorithm> // std::transform()
#include <iterator>  // std::back_inserter()
#include <string>
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
  class GeometryBuilder {
  public:
    // --- BEGIN Data types ----------------------------------------------------
    /// Identification of a single node in ROOT geometry.
    using Path_t = geo::GeoNodePath;

    /// Type of direct collection of geometry objects.
    template <typename GeoObj>
    using GeoColl_t = std::vector<GeoObj>;

    // --- END Data types ------------------------------------------------------

    // --- BEGIN Constructors and destructor -----------------------------------
    /// Virtual destructor.
    virtual ~GeometryBuilder() = default;

    // --- END Constructors and destructor -------------------------------------

    // --- BEGIN Cryostat information ------------------------------------------
    /// @name Cryostat information
    /// @{

    /// Collection of cryostat information objects.
    using Cryostats_t = GeoColl_t<geo::CryostatGeo>;

    /**
     * @brief Looks for all cryostats under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed cryostats
     *
     * The cryostats contain all their inner elements.  The current node itself of the
     * path is also considered as cryostat candidate, then it is descended into.
     */
    Cryostats_t extractCryostats(Path_t path) const { return doExtractCryostats(path); }

    /// @}
    // --- END Cryostat information --------------------------------------------

    // --- BEGIN Auxiliary detector information --------------------------------
    /// @name Auxiliary detector information
    /// @{

    /// Collection of auxiliary detector information objects.
    using AuxDets_t = GeoColl_t<geo::AuxDetGeo>;

    /**
     * @brief Looks for all auxiliary detectors under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed auxiliary detectors
     *
     * The auxiliary detectors contain all their inner elements.  The current node itself
     * of the path is also considered as auxiliary detector candidate, then it is
     * descended into.
     */
    AuxDets_t extractAuxiliaryDetectors(Path_t path) const
    {
      return doExtractAuxiliaryDetectors(path);
    }

    /// @}
    // --- END Auxiliary detector information ----------------------------------

    // --- END Static utility methods ------------------------------------------

  private:
    virtual Cryostats_t doExtractCryostats(Path_t& path) const = 0;
    virtual AuxDets_t doExtractAuxiliaryDetectors(Path_t& path) const = 0;

  }; // class GeometryBuilder

} // namespace geo

#endif // LARCOREALG_GEOMETRY_GEOMETRYBUILDER_H
