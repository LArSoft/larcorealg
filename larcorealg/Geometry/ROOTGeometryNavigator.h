/**
 * @file   larcorealg/Geometry/ROOTGeometryNavigator.h
 * @brief  Class representing a path in ROOT geometry.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/Geometry/GeometryBuilder.h`,
 *         `larcorealg/Geometry/GeoNodePath.cxx`
 */

#ifndef LARCOREALG_GEOMETRY_ROOTGEOMETRYNAVIGATOR_H
#define LARCOREALG_GEOMETRY_ROOTGEOMETRYNAVIGATOR_H

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h" // geo::...::makeFromCoords()
#include "larcorealg/Geometry/GeoNodePath.h"

// ROOT libraries
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

// C++ standard library
#include <cassert>
#include <utility> // std::forward()

namespace geo {

  class ROOTGeometryNavigator;

} // namespace geo

//------------------------------------------------------------------------------
/**
 * @brief Executes an operation on all the nodes of the ROOT geometry.
 *
 * For example, to collect the path (see `geo::GeoNodePath`) of all the volumes
 * named `"volTPC"` in the geometry loaded in LArSoft:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::string const volumeName = "volTPC";
 *
 * auto const& geom = *(lar::providerFrom<geo::Geometry>());
 *
 * geo::ROOTGeometryNavigator navigator { *(geom.ROOTGeoManager()) };
 *
 * // the operation executed on all nodes accumulates the paths in `volumePaths`
 * std::vector<geo::GeoNodePath> volumePaths;
 * auto findVolume = [&volumePaths, volumeName](auto& path)
 *   {
 *     if (path.current().GetVolume()->GetName() == volumeName)
 *       volumePaths.push_back(path);
 *     return true;
 *   };
 *
 * geo::ROOTGeometryNavigator navigator { *(geom.ROOTGeoManager()) };
 *
 * navigator.apply(findVolume);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 */
class geo::ROOTGeometryNavigator {

  TGeoNode const* fTopNode = nullptr;

public:
  /// Constructor: picks the manager.
  ROOTGeometryNavigator(TGeoManager const& manager) : fTopNode(manager.GetTopNode()) {}

  /**
   * @brief Applies the specified operation to all nodes under the `path`.
   * @tparam Op type of operation (see description)
   * @param path the path to the first node to operate on
   * @param op operation to be applied
   * @return whether all nodes in the path were processed
   *
   * The operation `Op` must be a callable accepting a `geo::GeoNodePath` immutable
   * argument and returning a value convertible to boolean.  If a call to `op` results
   * into a `false` value, the recursion is terminated and `false` is returned. `path`
   * will be pointing to the last node already processed.
   *
   * The node at the head of the path is processed first, then for each daughter node,
   * first the daughter itself then its own daughters, recursively.
   */
  template <typename Op>
  bool apply(GeoNodePath& path, Op&& op) const;

  /**
   * @brief Applies the specified operation to all nodes under `node`.
   * @tparam Op type of operation (see description)
   * @param node the node to start from
   * @param op operation to be applied
   * @return whether all nodes in the path were processed
   * @see `apply(geo::GeoNodePath&, Op&&) const`
   *
   * The operation `Op` must be a callable accepting a `geo::GeoNodePath` immutable
   * argument.
   */
  template <typename Op>
  bool apply(TGeoNode const* node, Op&& op) const;

  /**
   * @brief Applies the specified operation to all nodes.
   * @tparam Op type of operation (see description)
   * @param op operation to be applied
   * @return whether all nodes in the path were processed
   *
   * The operation `Op` must be a callable accepting a `geo::GeoNodePath` immutable
   * argument.
   */
  template <typename Op>
  bool apply(Op&& op) const;

}; // geo::ROOTGeometryNavigator

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Op>
bool geo::ROOTGeometryNavigator::apply(GeoNodePath& path, Op&& op) const
{
  if (!op(path)) return false;

  TGeoNode const* node = path.current();
  TGeoVolume const* pVolume = node->GetVolume();
  if (pVolume) { // is it even possible not to?
    int const nDaughters = pVolume->GetNdaughters();
    for (int iDaughter : util::counter<int>(nDaughters)) {
      TGeoNode const* pDaughter = pVolume->GetNode(iDaughter);
      if (!pDaughter) continue; // fishy...

      path.append(pDaughter);
      if (!apply(path, std::forward<Op>(op))) return false;
      path.pop();
    } // for
  }   // if we have a volume

  return true;
}

//------------------------------------------------------------------------------
template <typename Op>
bool geo::ROOTGeometryNavigator::apply(TGeoNode const* node, Op&& op) const
{
  GeoNodePath path{node};
  return apply(path, std::forward<Op>(op));
}

//------------------------------------------------------------------------------
template <typename Op>
bool geo::ROOTGeometryNavigator::apply(Op&& op) const
{
  assert(fTopNode);
  return apply(fTopNode, std::forward<Op>(op));
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_ROOTGEOMETRYNAVIGATOR_H
