////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.cxx
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include <utility>

namespace geo {
  bool GeoObjectSorter::compareOpDets(OpDetGeo const& od1, OpDetGeo const& od2) const
  {
    auto const [xyz1, xyz2] = std::make_pair(od1.GetCenter(), od2.GetCenter());

    if (xyz1.Z() != xyz2.Z()) return xyz1.Z() > xyz2.Z();
    if (xyz1.Y() != xyz2.Y()) return xyz1.Y() > xyz2.Y();
    return xyz1.X() > xyz2.X();
  }
}
