#ifndef LARCOREALG_GEOMETRY_GEOMETRYBUILDERSTANDARD_H
#define LARCOREALG_GEOMETRY_GEOMETRYBUILDERSTANDARD_H

// LArSoft libraries
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/GeometryExtractor.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TableFragment.h"

// C++ standard library
#include <string>

namespace geo {

  /**
   * @brief Extracts of LArSoft geometry information from ROOT.
   *
   * The builder manages several components, each devoted to the extraction of a specific
   * type of geometry object (e.g. cryostat, or wire plane within a TPC).
   *
   * Further customization notes
   * ===========================
   *
   * This builder does not extend the interface of `geo::GeometryBuilder`, but it defines
   * an interface that other builder classes can override to customize single elements of
   * the build. As long as the interface is complied to, the different components are
   * interchangeable.
   *
   * If instead a different interface is needed for one component, the parent component
   * needs to be customised too. For example, if the signature of `doExtractPlanes()` is
   * changed, also `doMakePlane()` needs to be customized to correctly call the
   * previous. In that case, take care of deleting the inherited interface to avoid
   * confusion and errors.
   *
   *
   * Technical notes on customization
   * --------------------------------
   *
   * The internal structure of the builder follows the pattern already employed in the
   * base class.  The base class defines both the public interface and the implementation,
   * but it separates the two leaving the former as non-virtual functions, and the latter
   * as virtual functions accessible only by derived classes.
   */

  class GeometryBuilderStandard : public GeometryBuilder {
  public:
    /// Configuration parameters.
    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::TableFragment<GeometryExtractor::Config> extractor;
      fhicl::Atom<std::string> opDetGeoName{
        Name("opDetGeoName"),
        Comment("the start of the name of optical detector GDML nodes"),
        "volOpDetSensitive" // default
      };

    }; // struct Config

    explicit GeometryBuilderStandard(fhicl::Table<Config> const& config);

    using AuxDetSensitive_t = GeoColl_t<AuxDetSensitiveGeo>;
    using OpDets_t = GeoColl_t<OpDetGeo>;
    using TPCs_t = GeoColl_t<TPCGeo>;

    GeometryExtractor fExtractObjects;
    std::string fOpDetGeoName; /// Name of the optical detector nodes.

    // --- BEGIN Auxiliary detector information --------------------------------
    /// @name Auxiliary detector information
    /// @{

    // extractAuxiliaryDetectors() and AuxDets_t are inherited public interface

    /// Constructs a `geo::AuxDetGeo` from the current node of the `path`.
    AuxDetGeo makeAuxDet(Path_t& path) const;

    /// Core implementation of `extractCryostats()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    AuxDets_t doExtractAuxiliaryDetectors(Path_t& path) const override;

    /// @}
    // --- END Auxiliary detector information ----------------------------------

    // --- BEGIN Auxiliary detector sensitive volume information ---------------
    /// @name Auxiliary detector sensitive volume information
    /// @{

    /**
     * @brief Looks for all auxiliary detectors under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed auxiliary detectors
     *
     * The auxiliary detectors contain all their inner elements.
     * The current node itself of the path is also considered as auxiliary
     * detector candidate, then it is descended into.
     *
     * @note Multithreading note: `path` is allowed to change during processing.
     */
    AuxDetSensitive_t extractAuxDetSensitive(Path_t& path) const;

    /// @}
    // --- END Auxiliary detector sensitive volume information -----------------

    // --- BEGIN Cryostat information ------------------------------------------
    /// @name Cryostat information
    /// @{

    // extractCryostats() and Cryostats_t are inherited public interface

    /// Constructs a `geo::CryostatGeo` from the current node of the `path`.
    CryostatGeo makeCryostat(Path_t& path) const;

    /// Core implementation of `extractCryostats()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    Cryostats_t doExtractCryostats(Path_t& path) const override;

    /// @}
    // --- END Cryostat information --------------------------------------------

    // --- BEGIN Optical detector information ----------------------------------
    /// @name Optical detector information
    /// @{

    /**
     * @brief Looks for all optical detectors under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed optical detector objects
     */
    OpDets_t extractOpDets(Path_t& path) const;

    /// @}
    // --- END Optical detector information ------------------------------------

    // --- BEGIN TPC information -----------------------------------------------
    /// @name TPC information
    /// @{

    /**
     * @brief Looks for all TPCs under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed TPC objects
     *
     * Each TPC has its own wire planes already in.
     */
    TPCs_t extractTPCs(Path_t& path) const;

    /// Constructs a `geo::TPCGeo` from the current node of the `path`.
    TPCGeo makeTPC(Path_t& path) const;

    /// @}
    // --- END TPC information -------------------------------------------------

  }; // class GeometryBuilderStandard

} // namespace geo

#endif // LARCOREALG_GEOMETRY_GEOMETRYBUILDERSTANDARD_H
