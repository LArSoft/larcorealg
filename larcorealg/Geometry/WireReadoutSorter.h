////////////////////////////////////////////////////////////////////////
/// @file  GeoObjectSorter.h
/// @brief Interface to algorithm class for sorting geo::XXXGeo objects
/// @ingroup Geometry
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_WIREREADOUTSORTER_H
#define LARCOREALG_GEOMETRY_WIREREADOUTSORTER_H

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class WireReadoutSorter {
  public:
    virtual ~WireReadoutSorter() = default;

    virtual bool compareWires(WireGeo const& w1, WireGeo const& w2) const = 0;
  };

}

#endif // LARCOREALG_GEOMETRY_WIREREADOUTSORTER_H
