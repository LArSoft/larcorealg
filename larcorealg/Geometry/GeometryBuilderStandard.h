/**
 * @file   larcorealg/Geometry/GeometryBuilderStandard.h
 * @brief  Standard implementation of geometry extractor.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/Geometry/GeometryBuilder.h`,
 *         `larcorealg/Geometry/GeometryBuilderStandard.cxx`
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYBUILDERSTANDARD_H
#define LARCOREALG_GEOMETRY_GEOMETRYBUILDERSTANDARD_H

// LArSoft libraries
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// support libraries
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
class TGeoNode;

// C++ standard library
#include <functional>
#include <limits> // std::numeric_limits<>
#include <string>

namespace geo {

  /**
   * @brief Extracts of LArSoft geometry information from ROOT.
   *
   * The builder manages several components, each devoted to the extraction of a specific
   * type of geometry object (e.g. cryostat, or wire plane within a TPC).
   *
   *
   * Further customization notes
   * ============================
   *
   * This builder does not extend the interface of `geo::GeometryBuilder`, but it defines
   * a protected interface that other builder classes could override to customize single
   * elements of the build. As long as the interface is complied to, the different
   * components are interchangeable.
   *
   * If instead a different interface is needed for one component, the parent component
   * needs to be customised too. For example, if the signature of `doExtractPlanes()` is
   * changed, also `doMakePlane()` needs to be customized to correctly call the
   * previous. In that case, take care of deleting the inherited interface to avoid
   * confusion and errors.
   *
   *
   * Technical notes on customization
   * ---------------------------------
   *
   * The internal structure of the builder follows the pattern already employed in the
   * base class.  The base class defines both the public interface and the implementation,
   * but it separates the two leaving the former as non-virtual functions, and the latter
   * as virtual functions accessible only by derived classes.
   *
   * The `geo::GeometryBuilderStandard` class replicates this pattern in a more hidden
   * level.  The general flow of the algorithm is a top-down crawl of the geometry tree
   * structure, where the top objects (cryostats and auxiliary detectors) are discovered
   * and built, and each of these objects takes care of discovering its own relevant
   * components. Therefore e.g. the cryostat algorithm will, once found a candidate
   * cryostat, descend into it to discover TPCs and optical detectors. This nested
   * discovery is delegated to other algorithms, and e.g. the TPC algorithm will take care
   * of creating a TPC and populating it with wire planes whose discovery is again
   * delegated to another algorithm.
   *
   * The interface of these algorithms is fixed and is part of the protected class
   * interface, in a way mirroring `geo::GeometryBuilder` in that it does not rely on
   * virtuality, but entirely protected. The implementation is also in the protected
   * space.
   *
   * Each component type has five elements:
   *
   * * a type describing a collection of the object of this component; this is integral
   *   part of the protected interface
   * * an interface to create an object for a single component, called `makeXxx()`,
   *   expected to rely on its implementation method `doMakeXxx()`
   * * an interface to discover all the components of this type, and creating them; it is
   *   called `extractXxx()` and it logically relies on `makeXxx()`, expected to rely on
   *   its implementation method `doExtractXxx()`
   * * a virtual implementation method of the component creation routine, called
   *   `doMakeXxx()` and expected to invoke the `extractYyy()` interface of all the
   *   subcomponents nested inside it
   * * a virtual implementation method of the component discovery routine, called
   *   `doExtractXxx()` and expected to invoke the `makeYyy()` interface of all the
   *   subcomponents nested inside it
   *
   * The discovery interface and the collection type of two of these components are
   * directly part of the public interface inherited from `geo::GeometryBuilder`.
   *
   */
  class GeometryBuilderStandard : public GeometryBuilder {
  public:
    /// Configuration parameters.
    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<Path_t::Depth_t> maxDepth{
        Name("maxDepth"),
        Comment("maximum number of level of the geometry structure to descend"),
        std::numeric_limits<Path_t::Depth_t>::max() // default
      };

      fhicl::Atom<std::string> opDetGeoName{
        Name("opDetGeoName"),
        Comment("the start of the name of optical detector GDML nodes"),
        "volOpDetSensitive" // default
      };

    }; // struct Config

    GeometryBuilderStandard(Config const& config);

  private:
    /// Maximum level to descend into in the path.
    Path_t::Depth_t fMaxDepth = std::numeric_limits<Path_t::Depth_t>::max();

    /// Name of the optical detector nodes.
    std::string fOpDetGeoName = "volOpDetSensitive";

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

    using AuxDetSensitive_t = GeoColl_t<AuxDetSensitiveGeo>;

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

    /// Constructs a `geo::AuxDetSensitiveGeo` from the current node of the `path`.
    AuxDetSensitiveGeo makeAuxDetSensitive(Path_t& path) const;

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

    using OpDets_t = GeoColl_t<OpDetGeo>;

    /**
     * @brief Looks for all optical detectors under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed optical detector objects
     */
    OpDets_t extractOpDets(Path_t& path) const;

    /// Constructs a `geo::OpDetGeo` from the current node of the `path`.
    OpDetGeo makeOpDet(Path_t& path) const;

    /// @}
    // --- END Optical detector information ------------------------------------

    // --- BEGIN TPC information -----------------------------------------------
    /// @name TPC information
    /// @{

    using TPCs_t = GeoColl_t<TPCGeo>;

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

    // --- BEGIN Plane information ---------------------------------------------
    /// @name Wire plane information
    /// @{

    using Planes_t = GeoColl_t<PlaneGeo>;

    /**
     * @brief Looks for all wire planes under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed wire plane objects
     *
     * Each plane has its own wires already in.
     */
    Planes_t extractPlanes(Path_t& path) const;

    /// Constructs a `geo::PlaneGeo` from the current node of the `path`.
    PlaneGeo makePlane(Path_t& path) const;

    /// @}
    // --- END Plane information -----------------------------------------------

    // --- BEGIN Wire information ----------------------------------------------
    /// @name Wire information
    /// @{

    using Wires_t = GeoColl_t<WireGeo>;

    /**
     * @brief Looks for all wires under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed wires
     *
     */
    Wires_t extractWires(Path_t& path) const;

    /// Constructs a `geo::WireGeo` from the current node of the `path`.
    WireGeo makeWire(Path_t& path) const;

    /// @}
    // --- END Wire information ------------------------------------------------

  private:
    /**
     * @brief Boilerplate implementation of `doExtractXxxx()` methods.
     * @tparam ObjGeo the geometry object being extracted (e.g. `geo::WireGeo`)
     * @param path the path to the node describing the object
     * @param IsObj function to identify if a node is of the right type
     * @param MakeObj class method creating the target object from a path
     * @return a fully constructed object of type `ObjGeo`
     *
     * This implementation first evaluates if the current node in the specified path is
     * suitable to create a `ObjGeo`; if not, then it descends into the node daughters and
     * recursively to their descendents.  For each candidate node, a `ObjGeo` is
     * created. All descendents of the candidates are ignored.
     *
     * @note Multithreading note: `path` is allowed to change during processing.
     */
    template <typename ObjGeo>
    GeoColl_t<ObjGeo> doExtractGeometryObjects(Path_t& path,
                                               std::function<bool(TGeoNode const&)> IsObj,
                                               ObjGeo (GeometryBuilderStandard::*MakeObj)(Path_t&)
                                                 const) const;

  }; // class GeometryBuilderStandard

} // namespace geo

#endif // LARCOREALG_GEOMETRY_GEOMETRYBUILDERSTANDARD_H
