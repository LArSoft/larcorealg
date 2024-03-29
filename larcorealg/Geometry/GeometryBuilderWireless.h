/**
 * @file   larcorealg/Geometry/GeometryBuilderWireless.h
 * @brief  Implementation of wireless geometry extractor.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 26, 2019
 * @see    `larcorealg/Geometry/GeometryBuilder.h`,
 *         `larcorealg/Geometry/GeometryBuilderStandard.cxx`
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYBUILDERWIRELESS_H
#define LARCOREALG_GEOMETRY_GEOMETRYBUILDERWIRELESS_H

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderStandard.h"

namespace geo {

  /**
   * @brief Geometry builder which ignores wires on wire planes.
   *
   * This builder works like `geo::GeometryBuilderStandard`, with the exception
   * that it does not consider the wires on the wire plane objects: wires may
   * or may not exist.
   *
   */
  class GeometryBuilderWireless : public geo::GeometryBuilderStandard {

  public:
    // import all constructors from base class
    using geo::GeometryBuilderStandard::GeometryBuilderStandard;

    //
    // we don't expand the public interface here
    //

  protected:
    // --- BEGIN Wire information ----------------------------------------------
    /// @name Wire information
    /// @{

    /// Core implementation of `extractWires()`: no wires returned whatsoever.
    ///
    /// The actual algorithm is specialization of `doExtractGeometryObjects()`.
    virtual Wires_t doExtractWires(Path_t&) { return {}; }

    /// @}
    // --- END Wire information ------------------------------------------------

  }; // class GeometryBuilderWireless

} // namespace geo

#endif // LARCOREALG_GEOMETRY_GEOMETRYBUILDERWIRELESS_H
