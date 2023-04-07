////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.h
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTER_H
#define GEO_GEOOBJECTSORTER_H

#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include <functional>
#include <vector>

namespace geo {

  class GeoObjectSorter {
  public:
    virtual ~GeoObjectSorter() = default;

    virtual bool compareAuxDets(AuxDetGeo const& ad1, AuxDetGeo const& ad2) const = 0;
    virtual bool compareAuxDetSensitives(AuxDetSensitiveGeo const& ads1,
                                         AuxDetSensitiveGeo const& ads2) const = 0;
    virtual bool compareCryostats(CryostatGeo const& c1, CryostatGeo const& c2) const = 0;
    virtual bool compareTPCs(TPCGeo const& t1, TPCGeo const& t2) const = 0;

    virtual bool compareOpDets(OpDetGeo const& od1, OpDetGeo const& od2) const;
  };

}

#endif // GEO_GEOOBJECTSORTER_H
