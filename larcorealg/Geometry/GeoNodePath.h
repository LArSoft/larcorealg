/**
 * @file   larcorealg/Geometry/GeoNodePath.h
 * @brief  Class representing a path in ROOT geometry.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/Geometry/GeometryBuilder.h`,
 *         `larcorealg/Geometry/GeoNodePath.cxx`
 */

#ifndef LARCOREALG_GEOMETRY_GEONODEPATH_H
#define LARCOREALG_GEOMETRY_GEONODEPATH_H

// LArSoft libraries
#include "larcorealg/Geometry/LocalTransformation.h"

// ROOT libraries
#include "TGeoNode.h"

// C++ standard library
#include <cstddef> // std::size_t
#include <string>
#include <vector>

namespace geo {

  /**
   * @brief Representation of a node and its ancestry.
   *
   * A `GeoNodePath` contains a sequence of nodes, from the `root()` node down to a
   * `current()` one.
   *
   * It behaves like a `stack` in that it inserts and removes elements at the "top", which
   * is also what defines the current node.
   *
   */
  class GeoNodePath {
  public:
    // FIXME: a better factorization is desirable
    using Entry = GeoNodePathEntry;

    // --- BEGIN Data types ----------------------------------------------------
    /// Type used to represent the depth of the path.
    using Depth_t = std::size_t;

    // --- END Data types ------------------------------------------------------

    // --- BEGIN Constructors and destructor -----------------------------------
    /// Sets all the the specified nodes into the current path.
    explicit GeoNodePath(TGeoNode const* topNode);
    // --- END Constructors and destructor -------------------------------------

    // --- BEGIN Query and access ----------------------------------------------
    /// Returns whether there is a current node.
    bool empty() const;

    /// Returns the depth of the path (elements including up to the current).
    Depth_t depth() const;

    /// Returns the current node. Undefined if the path is empty.
    TGeoNode const* current() const;

    /// Returns the current node. Undefined if the path is empty.
    Entry current_entry() const;

    /// Returns the parent entry of the current entry, or null if there is no parent.
    Entry parent_entry() const;

    // --- END Query and access ------------------------------------------------

    // --- BEGIN Content management --------------------------------------------
    /// Adds a node to the current path.
    void append(TGeoNode const* node);

    /// Removes the current node from the path, moving the current one up.
    void pop();
    // --- END Content management ----------------------------------------------

    /// Returns the total transformation to the current node, as a `Matrix`.
    template <typename Matrix = TGeoHMatrix>
    Matrix currentTransformation() const;

    /// Prints the full path (as node names) into a string.
    operator std::string() const;

  private:
    std::vector<Entry> fNodes; ///< Local path of pointers to ROOT geometry nodes.

  }; // class GeoNodePath

} // namespace geo

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Matrix /* = TGeoHMatrix */>
Matrix geo::GeoNodePath::currentTransformation() const
{
  return transformationFromPath<Matrix>(fNodes.begin(), fNodes.end());
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_GEONODEPATH_H
