#include "larcorealg/Geometry/details/maybe_default_detector_name.h"

#include "fhiclcpp/ParameterSet.h"

#include "boost/filesystem/path.hpp"

namespace geo::details {
  std::string maybe_default_detector_name(fhicl::ParameterSet const& pset,
                                          std::string const& filename)
  {
    if (auto name = pset.get_if_present<std::string>("Name")) { return *name; }
    return boost::filesystem::path{filename}.stem().native();
  }
}
