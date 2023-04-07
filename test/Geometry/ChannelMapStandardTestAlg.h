/**
 * @file   ChannelMapStandardTestAlg.h
 * @brief  Tests the standard channel mapping algorithm.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 26th, 2015
 */

#ifndef GEO_CHANNELMAPSTANDARDTESTALG_H
#define GEO_CHANNELMAPSTANDARDTESTALG_H

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class ChannelMapStandardTestAlg {
  public:
    explicit ChannelMapStandardTestAlg(GeometryCore const* new_geo,
                                       WireReadoutGeom const* wireReadoutGeom)
      : geom{new_geo}, wireReadoutGeom{wireReadoutGeom}
    {}

    /// Executes the test
    unsigned int Run();

    void TPCsetMappingTest() const;
    void ROPMappingTest() const;
    void ChannelMappingTest() const;

  private:
    GeometryCore const* geom;
    WireReadoutGeom const* wireReadoutGeom;

  }; // class ChannelMapStandardTestAlg

} // namespace geo

#endif // GEO_CHANNELMAPSTANDARDTESTALG_H
