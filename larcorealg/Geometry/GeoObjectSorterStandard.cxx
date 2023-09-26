////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterStandard.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

namespace geo {

  GeoObjectSorterStandard::GeoObjectSorterStandard() = default;
  GeoObjectSorterStandard::GeoObjectSorterStandard(fhicl::ParameterSet const&) {}

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
