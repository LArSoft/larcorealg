// LArSoft libraries
#include "larcorealg/Geometry/WireReadoutGeomBuilderStandard.h"

// ROOT libraries
#include "TGeoNode.h"

// C++ standard library
#include <algorithm> // std::move(...)
#include <string_view>

using namespace std::literals;

namespace {
  template <typename Dest, typename Src>
  void extendCollection(Dest& dest, Src&& src)
  {
    std::move(src.begin(), src.end(), std::back_inserter(dest));
  }

  /// Returns whether the start of `s` matches the full `key`.
  /// @note Remove this when C++20 is adopted (`s.starts_with(key)`).
  bool starts_with(std::string_view const s, std::string_view const key)
  {
    return s.compare(0, key.size(), key) == 0;
  }

  bool isPlaneNode(TGeoNode const& node) { return starts_with(node.GetName(), "volTPCPlane"sv); }
  bool isWireNode(TGeoNode const& node) { return starts_with(node.GetName(), "volTPCWire"sv); }
}

geo::WireReadoutGeomBuilderStandard::WireReadoutGeomBuilderStandard(
  fhicl::Table<Config> const& config)
  : fExtractObjects(config().extractor())
{}

geo::WireReadoutGeomBuilderStandard::WireReadoutGeomBuilderStandard(fhicl::ParameterSet const& pset)
  : WireReadoutGeomBuilderStandard{fhicl::Table<Config>{pset}}
{}
//------------------------------------------------------------------------------
auto geo::WireReadoutGeomBuilderStandard::doExtractPlanes(Path_t& path) const -> Planes_t
{
  Planes_t result;
  fExtractObjects(path, isPlaneNode, [&result, this](Path_t& path) {
    auto const [tpc, tpc_hash] = path.parent_entry();
    assert(tpc);
    result[tpc_hash].push_back(makePlane(path));
  });
  return result;
}

//------------------------------------------------------------------------------
geo::PlaneGeo geo::WireReadoutGeomBuilderStandard::makePlane(Path_t& path) const
{
  return {path.current(), path.currentTransformation<TransformationMatrix>(), extractWires(path)};
}

//------------------------------------------------------------------------------
std::vector<geo::WireGeo> geo::WireReadoutGeomBuilderStandard::extractWires(Path_t& path) const
{
  std::vector<WireGeo> result;
  fExtractObjects(
    path, isWireNode, [&result, this](Path_t const& path) { result.push_back(makeWire(path)); });
  return result;
}

//------------------------------------------------------------------------------
geo::WireGeo geo::WireReadoutGeomBuilderStandard::makeWire(Path_t const& path) const
{
  return {path.current(), path.currentTransformation<TransformationMatrix>()};
}
