#ifndef LARCOREALG_GEOMETRY_WIREREADOUTGEOMBUILDERSTANDARD_H
#define LARCOREALG_GEOMETRY_WIREREADOUTGEOMBUILDERSTANDARD_H

// LArSoft libraries
#include "larcorealg/Geometry/GeometryExtractor.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/WireReadoutGeometryBuilder.h"

// support libraries
#include "fhiclcpp/fwd.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TableFragment.h"

// ROOT libraries
class TGeoNode;

// C++ standard library
#include <functional>
#include <limits> // std::numeric_limits<>

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
  class WireReadoutGeomBuilderStandard : public WireReadoutGeometryBuilder {
  public:
    struct Config {
      fhicl::TableFragment<GeometryExtractor::Config> extractor;
    };
    explicit WireReadoutGeomBuilderStandard(fhicl::Table<Config> const& config);
    explicit WireReadoutGeomBuilderStandard(fhicl::ParameterSet const& pset);

  private:
    GeometryExtractor fExtractObjects;

    Planes_t doExtractPlanes(Path_t& path) const override;
    PlaneGeo makePlane(Path_t& path) const;

    std::vector<WireGeo> extractWires(Path_t& path) const;
    WireGeo makeWire(Path_t const& path) const;
  }; // class GeometryBuilderStandard

} // namespace geo

#endif // LARCOREALG_GEOMETRY_WIREREADOUTGEOMBUILDERSTANDARD_H
