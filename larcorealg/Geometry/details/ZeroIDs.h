#ifndef LARCOREALG_GEOMETRY_DETAILS_ZEREOIDS_H
#define LARCOREALG_GEOMETRY_DETAILS_ZEREOIDS_H

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace geo::details {
  inline constexpr CryostatID cryostat_zero{0};
  inline constexpr TPCID tpc_zero{cryostat_zero, 0};
}

#endif // LARCOREALG_GEOMETRY_DETAILS_ZEREOIDS_H
