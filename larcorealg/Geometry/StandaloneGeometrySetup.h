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
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"

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
   * @param channelMap channel mapping object to be used, already constructed
   * @return the geometry object, fully initialized
   * @see SetupGeometry()
   *
   * This function creates, sets up and returns a geometry object using the
   * specified channel mapping.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * // create a channel mapping algorithm
   * std::make_unique<geo::StandardChannelMapAlg> channelMap
   *   (pset.get<fhicl::ParameterSet>("SortingParameters"));
   *
   * std::unique_ptr<geo::GeometryCore> geom
   *   = SetupGeometryWithChannelMapping(pset, channelMap);
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
   *   `geo::ChannelMapAlg`); its content is dependent on the chosen
   *   implementation of `geo::ChannelMapAlg`
   */
  std::unique_ptr<geo::GeometryCore> GeometryFor(
    fhicl::ParameterSet const& pset,
    std::unique_ptr<geo::GeoObjectSorter const> sorter);

  template <typename ObjectSorter>
  std::unique_ptr<geo::GeometryCore> SetupGeometryWithoutChannelMapping(
    fhicl::ParameterSet const& pset)
  {
    return GeometryFor(
      pset, std::make_unique<ObjectSorter>(pset.get<fhicl::ParameterSet>("SortingParameters", {})));
  }

  //--------------------------------------------------------------------------
  /**
   * @brief  Initializes a LArSoft geometry object.
   * @tparam ChannelMapClass type of `geo::ChannelMapAlg` to be used
   * @tparam Args (optional) arguments for the construction of channel mapping
   * @param pset complete set of parameters for geometry configuration
   * @param args arguments to the channel mapping object constructor
   * @return the geometry object, fully initialized
   *
   * This function creates, sets up and returns a geometry object using the
   * specified channel mapping.
   * This is a simplified version of `SetupGeometryWithChannelMapping()`,
   * that can be used if no special treatment is needed for the channel
   * mapping after construction and before it is made to interact with the
   * geometry.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * // read FHiCL configuration from a configuration file:
   * fhicl::ParameterSet pset;
   * cet::filepath_lookup_after1 policy("FHICL_FILE_PATH");
   * pset = fhicl::_ParameterSet::make(configPath, policy);
   *
   * // set up message facility
   * mf::StartMessageFacility
   *   (pset.get<fhicl::ParameterSet>("services.message"));
   *
   * // geometry setup
   * std::unique_ptr<geo::GeometryCore> geom
   *   = SetupGeometry<geo::StandardChannelMapAlg>(pset);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Note that this function constructs the channel mapping object using a
   * constructor with arguments a parameter set and in addition, optionally,
   * any other argument specified in `args`.
   *
   */

  struct GeometryAndChannelMapAlg {
    std::unique_ptr<geo::GeometryCore> geometry;
    std::unique_ptr<geo::ChannelMapAlg> channelMapAlg;
  };

  template <typename ChannelMapClass, typename ObjectSorter, typename... Args>
  GeometryAndChannelMapAlg SetupGeometry(fhicl::ParameterSet const& pset, Args&&... args)
  {
    auto geom = SetupGeometryWithoutChannelMapping<ObjectSorter>(pset);
    auto channelMap = std::make_unique<ChannelMapClass>(geom.get(), std::forward<Args>(args)...);
    return {move(geom), move(channelMap)};
  }

  // --- END Geometry group --------------------------------------------------
  /// @}

} // namespace lar::standalone

#endif // LARCOREALG_GEOMETRY_STANDALONEGEOMETRYSETUP_H
