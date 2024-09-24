#include "larcorealg/Geometry/AuxDetGeoObjectSorter.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

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

  void AuxDetGeoObjectSorter::sort(std::vector<AuxDetGeo>& ads) const
  {
    std::sort(ads.begin(), ads.end(), bind(&AuxDetGeoObjectSorter::compareAuxDets, this));
  }

  void AuxDetGeoObjectSorter::sort(std::vector<AuxDetSensitiveGeo>& adss) const
  {
    std::sort(
      adss.begin(), adss.end(), bind(&AuxDetGeoObjectSorter::compareAuxDetSensitives, this));
  }

}
