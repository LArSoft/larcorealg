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
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// support libraries
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TGeoNode.h"

// C++ standard library
#include <string>
#include <string_view>
#include <limits> // std::numeric_limits<>


namespace geo {

  /**
   * @brief Extracts of LArSoft geometry information from ROOT.
   *
   * The builder manages several components, each devoted to the extraction of
   * a specific type of geometry object (e.g. cryostat, or wire plane within a
   * TPC).
   *
   *
   * Further customization notes
   * ============================
   *
   * This builder does not extend the interface of `geo::GeometryBuilder`, but
   * it defines a protected interface that other builder classes could override
   * to customize single elements of the build. As long as the interface is
   * complied to, the different components are interchangeable.
   *
   * If instead a different interface is needed for one component, the parent
   * component needs to be customised too. For example, if the signature of
   * `doExtractPlanes()` is changed, also `doMakePlane()` needs to be
   * customized to correctly call the previous. In that case, take care of
   * deleting the inherited interface to avoid confusion and errors.
   *
   *
   * Technical notes on customization
   * ---------------------------------
   *
   * The internal structure of the builder follows the pattern already employed
   * in the base class.
   * The base class defines both the public interface and the implementation,
   * but it separates the two leaving the former as non-virtual functions,
   * and the latter as virtual functions accessible only by derived classes.
   *
   * The `geo::GeometryBuilderStandard` class replicates this pattern in a more
   * hidden level.
   * The general flow of the algorithm is a top-down crawl of the geometry tree
   * structure, where the top objects (cryostats and auxiliary detectors) are
   * discovered and built, and each of these objects takes care of discovering
   * its own relevant components. Therefore e.g. the cryostat algorithm will,
   * once found a candidate cryostat, descend into it to discover TPCs and
   * optical detectors. This nested discovery is delegated to other algorithms,
   * and e.g. the TPC algorithm will take care of creating a TPC and populating
   * it with wire planes whose discovery is again delegated to another
   * algorithm.
   *
   * The interface of these algorithms is fixed and is part of the protected
   * class interface, in a way mirroring `geo::GeometryBuilder` in that it does
   * not rely on virtuality, but entirely protected. The implementation is also
   * in the protected space.
   *
   * Each component type has five elements:
   *
   * * a type describing a collection of the object of this component;
   *   this is integral part of the protected interface
   * * an interface to create an object for a single component, called
   *   `makeXxx()`, expected to rely on its implementation method `doMakeXxx()`
   * * an interface to discover all the components of this type, and creating
   *   them; it is called `extractXxx()` and it logically relies on `makeXxx()`,
   *   expected to rely on its implementation method `doExtractXxx()`
   * * a virtual implementation method of the component creation routine, called
   *   `doMakeXxx()` and expected to invoke the `extractYyy()` interface of
   *   all the subcomponents nested inside it
   * * a virtual implementation method of the component discovery routine,
   *   called `doExtractXxx()` and expected to invoke the `makeYyy()` interface
   *   of all the subcomponents nested inside it
   *
   * The discovery interface and the collection type of two of these components
   * are directly part of the public interface inherited from
   * `geo::GeometryBuilder`.
   *
   */
  class GeometryBuilderStandard: public GeometryBuilder {

      public:

    /// Configuration parameters.
    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<Path_t::Depth_t> maxDepth {
        Name("maxDepth"),
        Comment("maximum number of level of the geometry structure to descend"),
        std::numeric_limits<Path_t::Depth_t>::max() // default
        };

      fhicl::Atom<std::string> opDetGeoName {
        Name("opDetGeoName"),
        Comment("the start of the name of optical detector GDML nodes"),
        "volOpDetSensitive" // default
        };

    }; // struct Config

    GeometryBuilderStandard(Config const& config);

    //
    // we don't expand the public interface here
    //

      protected:

    /// Maximum level to descend into in the path.
    Path_t::Depth_t fMaxDepth = std::numeric_limits<Path_t::Depth_t>::max();

    /// Name of the optical detector nodes.
    std::string fOpDetGeoName = "volOpDetSensitive";


    // --- BEGIN Auxiliary detector information --------------------------------
    /// @name Auxiliary detector information
    /// @{

    // extractAuxiliaryDetectors() and AuxDets_t are inherited public interface

    /// Constructs a `geo::AuxDetGeo` from the current node of the `path`.
    geo::AuxDetGeo makeAuxDet(Path_t& path)
      { return doMakeAuxDet(path); }


    /// Core implementation of `extractCryostats()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    virtual AuxDets_t doExtractAuxiliaryDetectors(Path_t& path) override;

    /// Core implementation of `extractAuxiliaryDetectors()`.
    virtual geo::AuxDetGeo doMakeAuxDet(Path_t& path);

    /// @}
    // --- END Auxiliary detector information ----------------------------------


    // --- BEGIN Auxiliary detector sensitive volume information ---------------
    /// @name Auxiliary detector sensitive volume information
    /// @{

    using AuxDetSensitive_t = GeoColl_t<geo::AuxDetSensitiveGeo>;

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
    AuxDetSensitive_t extractAuxDetSensitive(Path_t& path)
      { auto localPath = path; return doExtractAuxDetSensitive(localPath); }

    /// Constructs a `geo::AuxDetSensitiveGeo` from the current node of the
    /// `path`.
    geo::AuxDetSensitiveGeo makeAuxDetSensitive(Path_t& path)
      { return doMakeAuxDetSensitive(path); }


    /// Core implementation of `extractAuxDetSensitive()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    virtual AuxDetSensitive_t doExtractAuxDetSensitive(Path_t& path);

    /// Core implementation of `makeAuxDetSensitive()`.
    virtual geo::AuxDetSensitiveGeo doMakeAuxDetSensitive(Path_t& path);

    /// @}
    // --- END Auxiliary detector sensitive volume information -----------------


    // --- BEGIN Cryostat information ------------------------------------------
    /// @name Cryostat information
    /// @{

    // extractCryostats() and Cryostats_t are inherited public interface

    /// Constructs a `geo::CryostatGeo` from the current node of the `path`.
    geo::CryostatGeo makeCryostat(Path_t& path)
      { return doMakeCryostat(path); }


    /// Core implementation of `extractCryostats()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    virtual Cryostats_t doExtractCryostats(Path_t& path) override;

    /// Core implementation of `extractAuxDetSensitive()`.
    virtual geo::CryostatGeo doMakeCryostat(Path_t& path);

    /// @}
    // --- END Cryostat information --------------------------------------------


    // --- BEGIN TPC information -----------------------------------------------
    /// @name TPC information
    /// @{

    using TPCs_t = GeoColl_t<geo::TPCGeo>;

    /**
     * @brief Looks for all TPCs under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed TPC objects
     *
     * Each TPC has its own wire planes already in.
     */
    TPCs_t extractTPCs(Path_t& path)
      { auto localPath = path; return doExtractTPCs(localPath); }

    /// Constructs a `geo::TPCGeo` from the current node of the `path`.
    geo::TPCGeo makeTPC(Path_t& path)
      { return doMakeTPC(path); }


    /// Core implementation of `extractTPCs()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    virtual TPCs_t doExtractTPCs(Path_t& path);

    /// Core implementation of `makeTPC()`.
    virtual geo::TPCGeo doMakeTPC(Path_t& path);

    /// @}
    // --- END TPC information -------------------------------------------------


    // --- BEGIN Plane information ---------------------------------------------
    /// @name Wire plane information
    /// @{

    using Planes_t = GeoColl_t<geo::PlaneGeo>;

    /**
     * @brief Looks for all wire planes under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed wire plane objects
     *
     * Each plane has its own wires already in.
     */
    Planes_t extractPlanes(Path_t& path)
      { auto localPath = path; return doExtractPlanes(localPath); }

    /// Constructs a `geo::PlaneGeo` from the current node of the `path`.
    geo::PlaneGeo makePlane(Path_t& path)
      { return doMakePlane(path); }


    /// Core implementation of `extractPlanes()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    virtual Planes_t doExtractPlanes(Path_t& path);

    /// Core implementation of `makePlanes()`.
    virtual geo::PlaneGeo doMakePlane(Path_t& path);

    /// @}
    // --- END Plane information -----------------------------------------------


    // --- BEGIN Wire information ----------------------------------------------
    /// @name Wire information
    /// @{

    using Wires_t = GeoColl_t<geo::WireGeo>;

    /**
     * @brief Looks for all wires under the specified path.
     * @param path path pointing to the starting node
     * @return a list of fully constructed wires
     *
     */
    Wires_t extractWires(Path_t& path)
      { auto localPath = path; return doExtractWires(localPath); }

    /// Constructs a `geo::WireGeo` from the current node of the `path`.
    geo::WireGeo makeWire(Path_t& path)
      { return doMakeWire(path); }


    /// Core implementation of `extractWires()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    virtual Wires_t doExtractWires(Path_t& path);

    /// Core implementation of `makeWire()`.
    virtual geo::WireGeo doMakeWire(Path_t& path);

    /// @}
    // --- END Wire information ------------------------------------------------


    // --- BEGIN Optical detector information ----------------------------------
    /// @name Optical detector information
    /// @{

    using OpDets_t = GeoColl_t<geo::OpDetGeo>;

    /**
     * @brief Looks for all optical detectors under the specified path.
     * @param path path pointing to the starting node
     * @param refFrame is non-null, use this as reference for directions
     * @return a list of fully constructed optical detector objects
     */
    OpDets_t extractOpDets(Path_t& path)
      { auto localPath = path; return doExtractOpDets(localPath); }

    /// Constructs a `geo::OpDetGeo` from the current node of the `path`.
    geo::OpDetGeo makeOpDet(Path_t& path)
      { return doMakeOpDet(path); }


    /// Core implementation of `extractOpDets()`.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    virtual OpDets_t doExtractOpDets(Path_t& path);

    /// Core implementation of `makeOpDet()`.
    virtual geo::OpDetGeo doMakeOpDet(Path_t& path);

    /// @}
    // --- END Optical detector information ------------------------------------



    // --- BEGIN Note type identification --------------------------------------
    /**
     * @name Node type identification
     *
     * These are implementation details of `doExtractGeometryObjects()` and its
     * users. They can be made virtual if the needs arises.
     */
    /// @{
    /// Returns whether the specified node is recognised as auxiliary detector.
    bool isAuxDetNode(TGeoNode const& node) const;

    /// Returns whether the specified node is recognised as sensitive volume of
    /// auxiliary detector.
    bool isAuxDetSensitiveNode(TGeoNode const& node) const;

    /// Returns whether the specified node is recognised as a cryostat.
    bool isCryostatNode(TGeoNode const& node) const;

    /// Returns whether the specified node is recognised as a optical detector.
    bool isOpDetNode(TGeoNode const& node) const;

    /// Returns whether the specified node is recognised as a TPC.
    bool isTPCNode(TGeoNode const& node) const;

    /// Returns whether the specified node is recognised as a wire plane.
    bool isPlaneNode(TGeoNode const& node) const;

    /// Returns whether the specified node is recognised as a wire.
    bool isWireNode(TGeoNode const& node) const;

    /// @}
    // --- END Note type identification ----------------------------------------


    /// Returns whether the start of `s` matches the full `key`.
    /// @note Remove this when C++20 is adopted (`s.starts_with(key)`).
    static bool starts_with
      (std::string_view const& s, std::string_view const& key)
      { return s.compare(0, key.size(), key) == 0; }


      private:

    /// Trait for type check on `MakeObj` template parameter.
    template <typename ObjGeo, typename T>
    struct IsMakeObjMemberFunctionType;
    
    /**
     * @brief Boilerplate implementation of `doExtractXxxx()` methods.
     * @tparam ObjGeo the geometry object being extracted (e.g. `geo::WireGeo`)
     * @tparam IsObj function to identify if a node is of the right type
     * @tparam MakeObj class method creating the target object from a path
     * @tparam MakeObjArgs type of additional arguments for `MakeObj`
     * @param path the path to the node describing the object
     * @return a fully constructed object of type `ObjGeo`
     *
     * This implementation first evaluates if the current node in the specified
     * path is suitable to create a `ObjGeo`; if not, then it descends into the
     * node daughters and recursively to their descendents.
     * For each candidate node, a `ObjGeo` is created. All descendents of the
     * candidates are ignored.
     *
     * @note Multithreading note: `path` is allowed to change during processing.
     */
    template <
      typename ObjGeo,
      bool (geo::GeometryBuilderStandard::*IsObj)(TGeoNode const&) const,
      auto MakeObj,
      typename... MakeObjArgs
      >
    GeoColl_t<ObjGeo> doExtractGeometryObjects
      (Path_t& path, MakeObjArgs&&... args);

    
  }; // class GeometryBuilderStandard

} // namespace geo


#endif // LARCOREALG_GEOMETRY_GEOMETRYBUILDERSTANDARD_H
