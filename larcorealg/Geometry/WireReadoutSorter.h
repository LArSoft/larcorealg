////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.h
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_WIREREADOUTSORTER_H
#define GEO_WIREREADOUTSORTER_H

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class WireReadoutSorter {
  public:
    virtual ~WireReadoutSorter() = default;

    virtual bool compareWires(WireGeo const& w1, WireGeo const& w2) const = 0;
  };

}

#endif // GEO_WIREREADOUTSORTER_H
