////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.h
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETGEOOBJECTSORTER_H
#define GEO_AUXDETGEOOBJECTSORTER_H

#include "larcorealg/Geometry/fwd.h"

#include <vector>

namespace geo {

  /// \ingroup Geometry
  class AuxDetGeoObjectSorter {
  public:
    virtual ~AuxDetGeoObjectSorter() = default;
    void sort(std::vector<AuxDetGeo>& ads) const;
    void sort(std::vector<AuxDetSensitiveGeo>& adss) const;

  private:
    virtual bool compareAuxDets(AuxDetGeo const& ad1, AuxDetGeo const& ad2) const = 0;
    virtual bool compareAuxDetSensitives(AuxDetSensitiveGeo const& ads1,
                                         AuxDetSensitiveGeo const& ads2) const = 0;
  };

}

#endif // GEO_AUXDETGEOOBJECTSORTER_H
