////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.cxx
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include <algorithm>
#include <vector>

namespace geo{

  static bool sortorderOpDets(const OpDetGeo& t1, const OpDetGeo& t2)
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.};
    t1.LocalToWorld(local, xyz1);
    t2.LocalToWorld(local, xyz2);

    if(xyz1[2] != xyz2[2])
      return xyz1[2] < xyz2[2];
    else if(xyz1[1] != xyz2[1])
      return xyz1[1] < xyz2[1];
    else
      return xyz1[0] < xyz2[0];
  }

  void GeoObjectSorter::SortOpDets(std::vector<geo::OpDetGeo> & opdet) const
  {
    std::sort(opdet.begin(), opdet.end(), sortorderOpDets);
  }
}
