////////////////////////////////////////////////////////////////////////
/// @file  GeoObjectSorterStandard.h
/// @brief Standard algorithm class for sorting of geo::XXXGeo objects
/// @ingroup Geometry
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_GEOOBJECTSORTERSTANDARD_H
#define LARCOREALG_GEOMETRY_GEOOBJECTSORTERSTANDARD_H

#include "fhiclcpp/fwd.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/fwd.h"

namespace geo {

  /// @ingroup Geometry
  class GeoObjectSorterStandard : public GeoObjectSorter {
  public:
    GeoObjectSorterStandard();
    explicit GeoObjectSorterStandard(fhicl::ParameterSet const&);

  private:
    bool compareCryostats(CryostatGeo const& c1, CryostatGeo const& c2) const override;
    bool compareTPCs(TPCGeo const& t1, TPCGeo const& t2) const override;
  };

}

#endif // LARCOREALG_GEOMETRY_GEOOBJECTSORTERSTANDARD_H
