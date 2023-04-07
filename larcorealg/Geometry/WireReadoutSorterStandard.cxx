#include "larcorealg/Geometry/WireReadoutSorterStandard.h"
#include "larcorealg/Geometry/WireGeo.h"

#include <utility>

namespace {
  constexpr double DistanceTol = 0.001; // cm
}

namespace geo {

  WireReadoutSorterStandard::WireReadoutSorterStandard() = default;
  WireReadoutSorterStandard::WireReadoutSorterStandard(fhicl::ParameterSet const&) {}

  //----------------------------------------------------------------------------
  bool WireReadoutSorterStandard::compareWires(WireGeo const& w1, WireGeo const& w2) const
  {
    auto const [xyz1, xyz2] = std::make_pair(w1.GetCenter(), w2.GetCenter());

    //sort by z first
    if (std::abs(xyz1.Z() - xyz2.Z()) > DistanceTol) return xyz1.Z() < xyz2.Z();

    //if same z sort by y
    if (std::abs(xyz1.Y() - xyz2.Y()) > DistanceTol) return xyz1.Y() < xyz2.Y();

    //if same y sort by x
    return xyz1.X() < xyz2.X();
  }

}
