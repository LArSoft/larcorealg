////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterStandard.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETGEOOBJECTSORTERSTANDARD_H
#define GEO_AUXDETGEOOBJECTSORTERSTANDARD_H

#include "fhiclcpp/fwd.h"
#include "larcorealg/Geometry/AuxDetGeoObjectSorter.h"
#include "larcorealg/Geometry/fwd.h"

namespace geo {

  /// @ingroup Geometry
  class AuxDetGeoObjectSorterStandard : public AuxDetGeoObjectSorter {
  public:
    AuxDetGeoObjectSorterStandard();
    explicit AuxDetGeoObjectSorterStandard(fhicl::ParameterSet const&);

  private:
    bool compareAuxDets(AuxDetGeo const& ad1, AuxDetGeo const& ad2) const override;
    bool compareAuxDetSensitives(AuxDetSensitiveGeo const& ad1,
                                 AuxDetSensitiveGeo const& ad2) const override;
  };

}

#endif // GEO_AUXDETGEOOBJECTSORTERSTANDARD_H
