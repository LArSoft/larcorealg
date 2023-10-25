////////////////////////////////////////////////////////////////////////
/// @file  GeoObjectSorterStandard.cxx
/// @brief Interface to algorithm class for sorting standard geo::XXXGeo objects
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/AuxDetGeoObjectSorterStandard.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

#include <cstdlib>
#include <string>

namespace geo {

  AuxDetGeoObjectSorterStandard::AuxDetGeoObjectSorterStandard() = default;
  AuxDetGeoObjectSorterStandard::AuxDetGeoObjectSorterStandard(fhicl::ParameterSet const&) {}

  //----------------------------------------------------------------------------
  bool AuxDetGeoObjectSorterStandard::compareAuxDets(AuxDetGeo const& ad1,
                                                     AuxDetGeo const& ad2) const
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string const ad1name = ad1.TotalVolume()->GetName();
    std::string const ad2name = ad2.TotalVolume()->GetName();

    // assume volume name is "volAuxDet##"
    int ad1Num = std::atoi(ad1name.substr(9, ad1name.size()).c_str());
    int ad2Num = std::atoi(ad2name.substr(9, ad2name.size()).c_str());

    return ad1Num < ad2Num;
  }

  //----------------------------------------------------------------------------
  bool AuxDetGeoObjectSorterStandard::compareAuxDetSensitives(AuxDetSensitiveGeo const& ad1,
                                                              AuxDetSensitiveGeo const& ad2) const
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string const ad1name = ad1.TotalVolume()->GetName();
    std::string const ad2name = ad2.TotalVolume()->GetName();

    // assume volume name is "volAuxDetSensitive##"
    int ad1Num = std::atoi(ad1name.substr(18, ad1name.size()).c_str());
    int ad2Num = std::atoi(ad2name.substr(18, ad2name.size()).c_str());

    return ad1Num < ad2Num;
  }

}
