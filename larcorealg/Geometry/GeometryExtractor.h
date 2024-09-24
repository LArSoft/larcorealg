/**
 * @file   larcorealg/Geometry/GeometryExtractor.h
 * @brief  Implementation of geometry extractor.
 * @see    `larcorealg/Geometry/GeometryBuilder.h`,
 *         `larcorealg/Geometry/GeometryBuilderStandard.cxx`
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYEXTRACTOR_H
#define LARCOREALG_GEOMETRY_GEOMETRYEXTRACTOR_H

// LArSoft libraries
#include "larcorealg/Geometry/GeoNodePath.h"

// art libraries
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TGeoNode.h"

// C++ standard library
#include <functional>
#include <limits> // std::numeric_limits<>

namespace geo {

  /**
   * @brief Object for extracting geometry objects from the GDML file.
   *
   * The general flow of the algorithm is a top-down crawl of the geometry tree structure,
   * where the top objects (cryostats and auxiliary detectors) are discovered and built,
   * and each of these objects takes care of discovering its own relevant
   * components. Therefore e.g. the cryostat algorithm will, once found a candidate
   * cryostat, descend into it to discover TPCs and optical detectors. This nested
   * discovery is delegated to other algorithms, and e.g. the TPC algorithm will take care
   * of creating a TPC and populating it with wire planes whose discovery is again
   * delegated to another algorithm.
   */

  class GeometryExtractor {
    using Path_t = GeoNodePath;

  public:
    struct Config {
      fhicl::Atom<Path_t::Depth_t> maxDepth{
        fhicl::Name("maxDepth"),
        fhicl::Comment("maximum number of level of the geometry structure to descend"),
        std::numeric_limits<Path_t::Depth_t>::max() // default
      };
    };

    explicit GeometryExtractor(Config const& config) : fMaxDepth{config.maxDepth()} {}

    /**
     * @brief Boilerplate implementation of geometry-extraction methods.
     * @tparam ObjGeo the geometry object being extracted (e.g. `geo::WireGeo`)
     * @param path the path to the node describing the object
     * @param IsObj function to identify if a node is of the right type
     * @param captureObject callable object invoked to create the target object from a path
     *
     * This implementation first evaluates if the current node in the specified path is
     * suitable to create a `ObjGeo`; if not, then it descends into the node daughters and
     * recursively to their descendents.  For each candidate node, a `ObjGeo` is
     * created. All descendents of the candidates are ignored.
     *
     * @note The recursive descent is an inefficiency whenever it is know that the
     * specified element will not be part of a given branch of the geometry tree.
     *
     * @note Multithreading note: `path` is allowed to change during processing.
     */
    template <typename FT>
    void operator()(Path_t& path,
                    std::function<bool(TGeoNode const&)> IsObj,
                    FT captureObject) const;

  private:
    /// Maximum level to descend into in the path.
    Path_t::Depth_t fMaxDepth;
  };

  // =================================================================================

  template <typename FT>
  void GeometryExtractor::operator()(Path_t& path,
                                     std::function<bool(TGeoNode const&)> const IsObj,
                                     FT captureObject) const
  {
    TGeoNode const* current = path.current();
    if (IsObj(*current)) {
      captureObject(path);
      return;
    }

    // descend into the next layer down, concatenate the results and return them
    if (path.depth() >= fMaxDepth) return; // yep, this is empty

    int const n = current->GetNdaughters();
    for (int i = 0; i < n; ++i) {
      path.append(current->GetDaughter(i));
      operator()(path, IsObj, captureObject);
      path.pop();
    }
  }

} // namespace geo

#endif // LARCOREALG_GEOMETRY_GEOMETRYEXTRACTOR_H
