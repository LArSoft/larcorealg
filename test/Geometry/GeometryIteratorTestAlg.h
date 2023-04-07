/**
 * @file   GeometryIteratorTestAlg.h
 * @brief  Tests the correct iteration of the geo::Geometry iterators
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   May 7th, 2015
 */

#ifndef GEO_GEOMETRYITERATORTESTALG_H
#define GEO_GEOMETRYITERATORTESTALG_H

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class GeometryIteratorTestAlg {
  public:
    explicit GeometryIteratorTestAlg(GeometryCore const* new_geo,
                                     WireReadoutGeom const* wireReadoutGeom)
      : geom{new_geo}, wireReadoutGeom{wireReadoutGeom}
    {}

    /// Executes the test
    unsigned int Run() const;

    /// @{
    /// @name ID iterator tests
    void CryostatIDIteratorsTest() const;
    void TPCIDIteratorsTest() const;
    void PlaneIDIteratorsTest() const;
    void WireIDIteratorsTest() const;
    void TPCsetIDIteratorsTest() const;
    void ROPIDIteratorsTest() const;
    /// @}

    /// @{
    /// @name Element iterator tests
    void CryostatIteratorsTest() const;
    void TPCIteratorsTest() const;
    void PlaneIteratorsTest() const;
    void WireIteratorsTest() const;
    /// @}

  private:
    GeometryCore const* geom;
    WireReadoutGeom const* wireReadoutGeom;
  };

} // namespace geo

#endif // GEO_GEOMETRYITERATORTESTALG_H
