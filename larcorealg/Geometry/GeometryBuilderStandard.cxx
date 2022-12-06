/**
 * @file   larcorealg/Geometry/GeometryBuilderStandard.cxx
 * @brief  Standard implementation of geometry extractor (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/Geometry/GeometryBuilderStandard.h`
 */

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderStandard.h"

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

  bool isAuxDetNode(TGeoNode const& node) { return starts_with(node.GetName(), "volAuxDet"sv); }
  bool isAuxDetSensitiveNode(TGeoNode const& node)
  {
    return std::string_view(node.GetName()).find("Sensitive") != std::string_view::npos;
  }
  bool isCryostatNode(TGeoNode const& node) { return starts_with(node.GetName(), "volCryostat"sv); }

  bool isOpDetNode(TGeoNode const& node, std::string_view const opDetGeoName)
  {
    return starts_with(node.GetName(), opDetGeoName);
  }
  bool isTPCNode(TGeoNode const& node) { return starts_with(node.GetName(), "volTPC"sv); }
  bool isPlaneNode(TGeoNode const& node) { return starts_with(node.GetName(), "volTPCPlane"sv); }
  bool isWireNode(TGeoNode const& node) { return starts_with(node.GetName(), "volTPCWire"sv); }
}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::GeometryBuilderStandard(Config const& config)
  : fMaxDepth(config.maxDepth()), fOpDetGeoName(config.opDetGeoName())
{}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::AuxDets_t geo::GeometryBuilderStandard::doExtractAuxiliaryDetectors(
  Path_t& path) const
{
  return doExtractGeometryObjects(path, isAuxDetNode, &GeometryBuilderStandard::makeAuxDet);
}

//------------------------------------------------------------------------------
geo::AuxDetGeo geo::GeometryBuilderStandard::makeAuxDet(Path_t& path) const
{
  return {path.current(),
          path.currentTransformation<TransformationMatrix>(),
          extractAuxDetSensitive(path)};
}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::AuxDetSensitive_t
geo::GeometryBuilderStandard::extractAuxDetSensitive(Path_t& path) const
{
  return doExtractGeometryObjects(
    path, isAuxDetSensitiveNode, &GeometryBuilderStandard::makeAuxDetSensitive);
}

//------------------------------------------------------------------------------
geo::AuxDetSensitiveGeo geo::GeometryBuilderStandard::makeAuxDetSensitive(Path_t& path) const
{
  return {path.current(), path.currentTransformation<TransformationMatrix>()};
}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::Cryostats_t geo::GeometryBuilderStandard::doExtractCryostats(
  Path_t& path) const
{
  return doExtractGeometryObjects(path, isCryostatNode, &GeometryBuilderStandard::makeCryostat);
}

//------------------------------------------------------------------------------
geo::CryostatGeo geo::GeometryBuilderStandard::makeCryostat(Path_t& path) const
{
  return {path.current(),
          path.currentTransformation<TransformationMatrix>(),
          extractTPCs(path),
          extractOpDets(path)};
}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::OpDets_t geo::GeometryBuilderStandard::extractOpDets(
  Path_t& path) const
{
  return doExtractGeometryObjects(
    path,
    [this](auto const& node) { return isOpDetNode(node, fOpDetGeoName); },
    &GeometryBuilderStandard::makeOpDet);
}

//------------------------------------------------------------------------------
geo::OpDetGeo geo::GeometryBuilderStandard::makeOpDet(Path_t& path) const
{
  return {path.current(), path.currentTransformation<TransformationMatrix>()};
}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::TPCs_t geo::GeometryBuilderStandard::extractTPCs(Path_t& path) const
{
  return doExtractGeometryObjects(path, isTPCNode, &GeometryBuilderStandard::makeTPC);
}

//------------------------------------------------------------------------------
geo::TPCGeo geo::GeometryBuilderStandard::makeTPC(Path_t& path) const
{
  return {path.current(), path.currentTransformation<TransformationMatrix>(), extractPlanes(path)};
}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::Planes_t geo::GeometryBuilderStandard::extractPlanes(
  Path_t& path) const
{
  return doExtractGeometryObjects(path, isPlaneNode, &GeometryBuilderStandard::makePlane);
}

//------------------------------------------------------------------------------
geo::PlaneGeo geo::GeometryBuilderStandard::makePlane(Path_t& path) const
{
  return {path.current(), path.currentTransformation<TransformationMatrix>(), extractWires(path)};
}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::Wires_t geo::GeometryBuilderStandard::extractWires(Path_t& path) const
{
  return doExtractGeometryObjects(path, isWireNode, &GeometryBuilderStandard::makeWire);
}

//------------------------------------------------------------------------------
geo::WireGeo geo::GeometryBuilderStandard::makeWire(Path_t& path) const
{
  return {path.current(), path.currentTransformation<TransformationMatrix>()};
}

//------------------------------------------------------------------------------
template <typename ObjGeo>
geo::GeometryBuilder::GeoColl_t<ObjGeo> geo::GeometryBuilderStandard::doExtractGeometryObjects(
  Path_t& path,
  std::function<bool(TGeoNode const&)> const IsObj,
  ObjGeo (GeometryBuilderStandard::*MakeObj)(Path_t&) const) const
{
  GeoColl_t<ObjGeo> objs;

  //
  // if this is a wire, we are set
  //
  if (IsObj(*path.current())) {
    objs.push_back((this->*MakeObj)(path));
    return objs;
  }

  //
  // descend into the next layer down, concatenate the results and return them
  //
  if (path.depth() >= fMaxDepth) return objs; // yep, this is empty

  TGeoVolume const* volume = path.current()->GetVolume();
  int const n = volume->GetNdaughters();
  for (int i = 0; i < n; ++i) {
    path.append(volume->GetNode(i));
    extendCollection(objs, doExtractGeometryObjects(path, IsObj, MakeObj));
    path.pop();
  } // for

  return objs;
}

//------------------------------------------------------------------------------
