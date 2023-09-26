#ifndef GEO_DETAILS_MAYBE_DEFAULT_DETECTOR_NAME_H
#define GEO_DETAILS_MAYBE_DEFAULT_DETECTOR_NAME_H

// ==================================================
// Given a filename, remove the suffix (if present)
// ==================================================

#include "fhiclcpp/fwd.h"

#include <string>

namespace geo::details {
  std::string maybe_default_detector_name(fhicl::ParameterSet const& pset,
                                          std::string const& filename);
}

#endif // GEO_DETAILS_MAYBE_DEFAULT_DETECTOR_NAME_H
