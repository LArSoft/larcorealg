////////////////////////////////////////////////////////////////////////
/// @file  WireReadoutSorterStandard.h
/// @brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
/// @ingroup Geometry
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_WIREREADOUTSORTERSTANDARD_H
#define LARCOREALG_GEOMETRY_WIREREADOUTSORTERSTANDARD_H

#include "fhiclcpp/fwd.h"
#include "larcorealg/Geometry/WireReadoutSorter.h"
#include "larcorealg/Geometry/fwd.h"

namespace geo {
  class WireReadoutSorterStandard : public WireReadoutSorter {
  public:
    WireReadoutSorterStandard();
    explicit WireReadoutSorterStandard(fhicl::ParameterSet const&);

  private:
    bool compareWires(WireGeo const& w1, WireGeo const& w2) const override;
  };
}

#endif // LARCOREALG_GEOMETRY_WIREREADOUTSORTERSTANDARD_H
