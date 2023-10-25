////////////////////////////////////////////////////////////////////////
/// @file  GeoObjectSorter.cxx
/// @brief Interface to algorithm class for sorting geo::XXXGeo objects
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

#include <algorithm>
#include <utility>

namespace {
  template <typename T, typename Arg>
  auto bind(bool (T::*ft)(Arg const&, Arg const&) const, T const* t)
  {
    return [t, ft](auto const& a, auto const& b) { return (t->*ft)(a, b); };
  }
}

namespace geo {

  void GeoObjectSorter::sort(std::vector<CryostatGeo>& cryostats) const
  {
    std::sort(cryostats.begin(), cryostats.end(), bind(&GeoObjectSorter::compareCryostats, this));
  }

  void GeoObjectSorter::sort(std::vector<TPCGeo>& tpcs) const
  {
    std::sort(tpcs.begin(), tpcs.end(), bind(&GeoObjectSorter::compareTPCs, this));
  }

  void GeoObjectSorter::sort(std::vector<OpDetGeo>& ods) const
  {
    std::sort(ods.begin(), ods.end(), bind(&GeoObjectSorter::compareOpDets, this));
  }

  bool GeoObjectSorter::compareOpDets(OpDetGeo const& od1, OpDetGeo const& od2) const
  {
    auto const [xyz1, xyz2] = std::make_pair(od1.GetCenter(), od2.GetCenter());

    if (xyz1.Z() != xyz2.Z()) return xyz1.Z() > xyz2.Z();
    if (xyz1.Y() != xyz2.Y()) return xyz1.Y() > xyz2.Y();
    return xyz1.X() > xyz2.X();
  }
}
