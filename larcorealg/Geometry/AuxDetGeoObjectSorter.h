////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.h
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETGEOOBJECTSORTER_H
#define GEO_AUXDETGEOOBJECTSORTER_H

#include <vector>

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  /// \ingroup Geometry
  class AuxDetGeoObjectSorter {
  public:
    virtual ~AuxDetGeoObjectSorter() = default;

    virtual void SortAuxDets(std::vector<AuxDetGeo>& adgeo) const = 0;
    virtual void SortAuxDetSensitive(std::vector<AuxDetSensitiveGeo>& adsgeo) const = 0;
  };

}

#endif // GEO_GEOOBJECTSORTER_H
