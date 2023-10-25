////////////////////////////////////////////////////////////////////////
/// @file  GeoObjectSorterStandard.h
/// @brief Standard algorithm class for sorting of geo::AuxDet objects
/// @ingroup Geometry
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_AUXDETGEOOBJECTSORTERSTANDARD_H
#define LARCOREALG_GEOMETRY_AUXDETGEOOBJECTSORTERSTANDARD_H

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

#endif // LARCOREALG_GEOMETRY_AUXDETGEOOBJECTSORTERSTANDARD_H
