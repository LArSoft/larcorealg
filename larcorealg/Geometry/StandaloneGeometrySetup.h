/**
 * @file   larcorealg/Geometry/StandaloneGeometrySetup.h
 * @brief  Utilities for one-line geometry initialization.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 22, 2017
 * @ingroup Geometry
 *
 * The main entry point for initializing the geometry is `SetupGeometry()`.
 *
 */

#ifndef LARCOREALG_GEOMETRY_STANDALONEGEOMETRYSETUP_H
#define LARCOREALG_GEOMETRY_STANDALONEGEOMETRYSETUP_H

// LArSoft libraries
#include "larcorealg/Geometry/AuxDetGeoObjectSorterStandard.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcorealg/Geometry/WireReadoutSorterStandard.h"
#include "larcorealg/Geometry/WireReadoutStandardGeom.h"

// art-provided libraries
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory> // std::make_unique()
#include <set>
#include <string>

namespace lar::standalone {

  // --- BEGIN Geometry group ------------------------------------------------
  /// @ingroup Geometry
  /// @{

  //--------------------------------------------------------------------------
  /**
   * @brief  Initializes a LArSoft geometry object.
   * @param pset parameters for geometry configuration
   * @param wireReadoutGeom channel mapping object to be used, already constructed
   * @return the geometry object, fully initialized
   * @see SetupGeometry()
   *
   * This function creates, sets up and returns a geometry object using the
   * specified channel mapping.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * // create a channel mapping algorithm
   * std::make_unique<geo::StandardWireReadoutGeom> wireReadoutGeom
   *   (pset.get<fhicl::ParameterSet>("SortingParameters"));
   *
   * std::unique_ptr<geo::GeometryCore> geom
   *   = SetupGeometryWithChannelMapping(pset, wireReadoutGeom);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * If no set up is required for channel mapping after construction, the use
   * of `SetupGeometry()` is preferred over this function.
   *
   *
   * Configuration parameters
   * =========================
   *
   * It is expected that a standard `geo::Geometry` service configuration will
   * correctly set up the geometry.
   *
   * In addition to the parameters documented in `geo::GeometryCore`, the
   * following parameters are supported:
   *
   * - *RelativePath* (string, default: no path): this path is prepended to
   *   the geometry file names before searching from them; the path string
   *   does not affect the file name
   * - *GDML* (string, mandatory): path of the GDML file to be served to
   *   GEANT4 *   for detector simulation. The full file is composed out of
   *   the optional relative path specified by `RelativePath` path and the
   *   base name specified in `GDML` parameter; this path is searched for in
   *   the directories configured in the `FW_SEARCH_PATH` environment
   *   variable;
   * - *ROOT* (string, mandatory): currently overridden by `GDML` parameter,
   *   whose value is used instead;
   *   this path is assembled in the same way as the one for `GDML` parameter,
   *   except that no alternative (wireless) geometry is used even if
   *   `DisableWiresInG4` is specified (see below); this file is used to load
   *   the geometry used in the internal simulation and reconstruction,
   *   basically everywhere except for the GEANT4 simulation
   * - *DisableWiresInG4* (boolean, default: false): if true, GEANT4 is loaded
   *   with an alternative geometry from a file with the standard name as
   *   configured with the /GDML/ parameter, but with an additional `_nowires`
   *   appended before the `.gdml` suffix
   * - *SortingParameters* (a parameter set; default: empty): this
   *   configuration is directly passed to the channel mapping algorithm (see
   *   `geo::WireReadoutGeom`); its content is dependent on the chosen
   *   implementation of `geo::WireReadoutGeom`
   */
  std::unique_ptr<geo::GeometryCore> GeometryFor(fhicl::ParameterSet const& pset,
                                                 std::unique_ptr<geo::GeoObjectSorter> sorter);

  std::unique_ptr<geo::AuxDetGeometryCore> AuxDetGeometryFor(
    fhicl::ParameterSet const& pset,
    std::unique_ptr<geo::AuxDetGeoObjectSorter> sorter,
    std::unique_ptr<geo::AuxDetInitializer> initializer);

  namespace detail {
    template <typename T>
    auto make_unique_maybe_default(fhicl::ParameterSet const& pset)
    {
      if (pset.is_empty()) {
        if constexpr (std::is_constructible_v<T>) { return std::make_unique<T>(); }
      }
      return std::make_unique<T>(pset);
    }
  }

  template <typename ObjectSorter = geo::GeoObjectSorterStandard>
  std::unique_ptr<geo::GeometryCore> SetupGeometry(fhicl::ParameterSet const& pset)
  {
    auto sorting_parameters = pset.get<fhicl::ParameterSet>("SortingParameters", {});
    return GeometryFor(pset, detail::make_unique_maybe_default<ObjectSorter>(sorting_parameters));
  }

  template <typename ObjectSorter = geo::AuxDetGeoObjectSorterStandard>
  std::unique_ptr<geo::AuxDetGeometryCore> SetupAuxDetGeometry(
    fhicl::ParameterSet const& pset,
    std::unique_ptr<geo::AuxDetInitializer> initializer = nullptr)
  {
    auto sorting_parameters = pset.get<fhicl::ParameterSet>("SortingParameters", {});
    return AuxDetGeometryFor(pset,
                             detail::make_unique_maybe_default<ObjectSorter>(sorting_parameters),
                             std::move(initializer));
  }

  template <typename ObjectSorter = geo::WireReadoutSorterStandard,
            typename WireGeom = geo::WireReadoutStandardGeom>
  std::unique_ptr<geo::WireReadoutGeom> SetupReadout(fhicl::ParameterSet const& pset,
                                                     geo::GeometryCore const* geom)
  {
    auto sorting_parameters = pset.get<fhicl::ParameterSet>("SortingParameters", {});
    return std::make_unique<WireGeom>(
      pset, geom, detail::make_unique_maybe_default<ObjectSorter>(sorting_parameters));
  }

  // --- END Geometry group --------------------------------------------------
  /// @}

} // namespace lar::standalone

#endif // LARCOREALG_GEOMETRY_STANDALONEGEOMETRYSETUP_H
