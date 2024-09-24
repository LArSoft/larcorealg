/**
 * @file   GeometryGeoIDTestAlg.h
 * @brief  Tests the correct assignment of IDs to detector geometry objects
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   October 31, 2016
 */

#ifndef TEST_GEOMETRY_GEOMETRYGEOIDTESTALG_H
#define TEST_GEOMETRY_GEOMETRYGEOIDTESTALG_H

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class GeometryGeoIDTestAlg {
  public:
    explicit GeometryGeoIDTestAlg(GeometryCore const* new_geom, WireReadoutGeom const* channel_map)
      : geom{new_geom}, wireReadoutGeom{channel_map}
    {}

    /// Executes the test
    unsigned int Run() const;

    /// @name All the ID iterator tests
    /// @{
    void CryostatGeoIDTest() const;
    void TPCGeoIDTest() const;
    void PlaneGeoIDTest() const;
    void WireGeoIDTest() const;
    /// @}

  private:
    GeometryCore const* geom = nullptr;
    WireReadoutGeom const* wireReadoutGeom = nullptr;
  }; // class GeometryGeoIDTestAlg

} // namespace geo

#endif // TEST_GEOMETRY_GEOMETRYGEOIDTESTALG_H
