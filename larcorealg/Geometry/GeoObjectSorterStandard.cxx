////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterStandard.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include <utility>

namespace {
  constexpr double DistanceTol = 0.001; // cm
}

namespace geo {

  GeoObjectSorterStandard::GeoObjectSorterStandard() = default;
  GeoObjectSorterStandard::GeoObjectSorterStandard(fhicl::ParameterSet const&) {}

  //----------------------------------------------------------------------------
  bool GeoObjectSorterStandard::compareAuxDets(AuxDetGeo const& ad1, AuxDetGeo const& ad2) const
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string const ad1name = ad1.TotalVolume()->GetName();
    std::string const ad2name = ad2.TotalVolume()->GetName();

    // assume volume name is "volAuxDet##"
    int ad1Num = atoi(ad1name.substr(9, ad1name.size()).c_str());
    int ad2Num = atoi(ad2name.substr(9, ad2name.size()).c_str());

    return ad1Num < ad2Num;
  }

  //----------------------------------------------------------------------------
  bool GeoObjectSorterStandard::compareAuxDetSensitives(AuxDetSensitiveGeo const& ad1,
                                                        AuxDetSensitiveGeo const& ad2) const
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string const ad1name = ad1.TotalVolume()->GetName();
    std::string const ad2name = ad2.TotalVolume()->GetName();

    // assume volume name is "volAuxDetSensitive##"
    int ad1Num = atoi(ad1name.substr(18, ad1name.size()).c_str());
    int ad2Num = atoi(ad2name.substr(18, ad2name.size()).c_str());

    return ad1Num < ad2Num;
  }

  //----------------------------------------------------------------------------
  bool GeoObjectSorterStandard::compareCryostats(CryostatGeo const& c1, CryostatGeo const& c2) const
  {
    CryostatGeo::LocalPoint_t const local{0., 0., 0.};
    auto const xyz1 = c1.toWorldCoords(local);
    auto const xyz2 = c2.toWorldCoords(local);
    return xyz1.X() < xyz2.X();
  }

  //----------------------------------------------------------------------------
  bool GeoObjectSorterStandard::compareTPCs(TPCGeo const& t1, TPCGeo const& t2) const
  {
    return t1.GetCenter().X() < t2.GetCenter().X();
  }

}
