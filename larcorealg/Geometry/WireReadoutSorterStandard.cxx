#include "larcorealg/Geometry/WireReadoutSorterStandard.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/Geometry/WireGeo.h"

#include <utility>

namespace {
  constexpr lar::util::RealComparisons<double> cmp{0.001 /*threshold in cm*/};
}

namespace geo {

  WireReadoutSorterStandard::WireReadoutSorterStandard() = default;
  WireReadoutSorterStandard::WireReadoutSorterStandard(fhicl::ParameterSet const&) {}

  //----------------------------------------------------------------------------
  bool WireReadoutSorterStandard::compareWires(WireGeo const& w1, WireGeo const& w2) const
  {
    auto const [c1, c2] = std::pair{w1.GetCenter(), w2.GetCenter()};

    // sort by z first
    if (cmp.nonEqual(c1.Z(), c2.Z())) return c1.Z() < c2.Z();

    // if same z sort by y
    if (cmp.nonEqual(c1.Y(), c2.Y())) return c1.Y() < c2.Y();

    // if same y sort by x
    return c1.X() < c2.X();
  }

}
