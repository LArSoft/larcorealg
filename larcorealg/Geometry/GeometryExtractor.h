/**
 * @file   larcorealg/Geometry/GeometryBuilderStandard.h
 * @brief  Standard implementation of geometry extractor.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
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
     * @brief Boilerplate implementation of `doExtractXxxx()` methods.
     * @tparam ObjGeo the geometry object being extracted (e.g. `geo::WireGeo`)
     * @param path the path to the node describing the object
     * @param IsObj function to identify if a node is of the right type
     * @param MakeObj class method creating the target object from a path
     *
     * This implementation first evaluates if the current node in the specified path is
     * suitable to create a `ObjGeo`; if not, then it descends into the node daughters and
     * recursively to their descendents.  For each candidate node, a `ObjGeo` is
     * created. All descendents of the candidates are ignored.
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
