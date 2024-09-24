/**
 * @file   larcorealg/Geometry/GeometryBuilderStandard.cxx
 * @brief  Standard implementation of geometry extractor (implementation file).
 * @see    `larcorealg/Geometry/GeometryBuilderStandard.h`
 */

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderStandard.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// ROOT libraries
#include "TGeoNode.h"

// C++ standard library
#include <algorithm>
#include <string_view>
#include <tuple>

using namespace std::literals;

namespace {
  /// Returns whether the start of `s` matches the full `key`.
  /// @note Remove this when C++20 is adopted (`s.starts_with(key)`).
  bool starts_with(std::string_view const s, std::string_view const key)
  {
    return s.compare(0, key.size(), key) == 0;
  }

  bool isAuxDetNode(TGeoNode const& node)
  {
    return starts_with(node.GetName(), "volAuxDet"sv);
  }
  bool isAuxDetSensitiveNode(TGeoNode const& node)
  {
    return std::string_view(node.GetName()).find("Sensitive") != std::string_view::npos;
  }
  bool isCryostatNode(TGeoNode const& node)
  {
    return starts_with(node.GetName(), "volCryostat"sv);
  }

  bool isOpDetNode(TGeoNode const& node, std::string_view const opDetGeoName)
  {
    return starts_with(node.GetName(), opDetGeoName);
  }
  bool isTPCNode(TGeoNode const& node)
  {
    return starts_with(node.GetName(), "volTPC"sv);
  }
  bool isPlaneNode(TGeoNode const& node)
  {
    return starts_with(node.GetName(), "volTPCPlane"sv);
  }

  geo::DriftSign PosOrNeg(bool const condition)
  {
    return condition ? geo::DriftSign::Positive : geo::DriftSign::Negative;
  }

  //......................................................................
  geo::DriftAxis DriftAxisWithSign(geo::Point_t const& TPCcenter,
                                   std::vector<geo::PlaneGeo> const& planes)
  {
    auto const PlaneCenter = planes[0].GetBoxCenter(); // any will do

    auto const driftVector = PlaneCenter - TPCcenter; // approximation!

    if ((std::abs(driftVector.X()) > std::abs(driftVector.Y())) &&
        (std::abs(driftVector.X()) > std::abs(driftVector.Z()))) {
      return {geo::Coordinate::X, PosOrNeg(driftVector.X() > 0)};
    }
    if (std::abs(driftVector.Y()) > std::abs(driftVector.Z())) {
      return {geo::Coordinate::Y, PosOrNeg(driftVector.Y() > 0)};
    }
    return {geo::Coordinate::Z, PosOrNeg(driftVector.Z() > 0)};
  }

  //------------------------------------------------------------------------------
  geo::AuxDetSensitiveGeo makeAuxDetSensitive(geo::GeoNodePath const& path)
  {
    return {path.current(), path.currentTransformation<geo::TransformationMatrix>()};
  }

  geo::OpDetGeo makeOpDet(geo::GeoNodePath const& path)
  {
    return {path.current(), path.currentTransformation<geo::TransformationMatrix>()};
  }

  geo::PlaneGeo makePlane(geo::GeoNodePath const& path)
  {
    // Planes are constructed only to provide information for the TPC constructor.  We
    // therefore can supply an empty wires collection with them.
    return {path.current(),
            path.currentTransformation<geo::TransformationMatrix>(),
            std::vector<geo::WireGeo>{}};
  }

  double driftDistance(geo::Point_t result,
                       geo::DriftAxis const drift_axis_with_sign,
                       TGeoBBox const* active_volume_box,
                       geo::Point_t last_plane_box_center)
  {
    auto const [axis, sign] = drift_axis_with_sign;
    auto const [offX, offY, offZ] = std::tuple{
      active_volume_box->GetDX(), active_volume_box->GetDY(), active_volume_box->GetDZ()};

    switch (axis) {
    case geo::Coordinate::X: {
      double const cathode_x = result.X() - to_int(sign) * offX;
      return std::abs(cathode_x - last_plane_box_center.X());
    }
    case geo::Coordinate::Y: {
      double const cathode_y = result.Y() - to_int(sign) * offY;
      return std::abs(cathode_y - last_plane_box_center.Y());
    }
    case geo::Coordinate::Z: {
      double const cathode_z = result.Z() - to_int(sign) * offZ;
      return std::abs(cathode_z - last_plane_box_center.Z());
    }
    }

    // unreachable
    return {};
  }
}

//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::GeometryBuilderStandard(fhicl::Table<Config> const& config)
  : fExtractObjects(config().extractor()), fOpDetGeoName(config().opDetGeoName())
{}

//------------------------------------------------------------------------------
auto geo::GeometryBuilderStandard::doExtractAuxiliaryDetectors(Path_t& path) const -> AuxDets_t
{
  std::vector<AuxDetGeo> result;
  fExtractObjects(
    path, isAuxDetNode, [&result, this](Path_t& path) { result.push_back(makeAuxDet(path)); });
  return result;
}

//------------------------------------------------------------------------------
geo::AuxDetGeo geo::GeometryBuilderStandard::makeAuxDet(Path_t& path) const
{
  return {path.current(),
          path.currentTransformation<TransformationMatrix>(),
          extractAuxDetSensitive(path)};
}

//------------------------------------------------------------------------------
auto geo::GeometryBuilderStandard::extractAuxDetSensitive(Path_t& path) const -> AuxDetSensitive_t
{
  std::vector<AuxDetSensitiveGeo> result;
  fExtractObjects(path, isAuxDetSensitiveNode, [&result](Path_t const& path) {
    result.push_back(makeAuxDetSensitive(path));
  });
  return result;
}

//------------------------------------------------------------------------------
auto geo::GeometryBuilderStandard::doExtractCryostats(Path_t& path) const -> Cryostats_t
{
  std::vector<CryostatGeo> result;
  fExtractObjects(
    path, isCryostatNode, [&result, this](Path_t& path) { result.push_back(makeCryostat(path)); });
  return result;
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
auto geo::GeometryBuilderStandard::extractOpDets(Path_t& path) const -> OpDets_t
{
  std::vector<OpDetGeo> result;
  fExtractObjects(
    path,
    [this](auto const& node) { return isOpDetNode(node, fOpDetGeoName); },
    [&result](Path_t const& path) { result.push_back(makeOpDet(path)); });
  return result;
}

//------------------------------------------------------------------------------
auto geo::GeometryBuilderStandard::extractTPCs(Path_t& path) const -> TPCs_t
{
  std::vector<TPCGeo> result;
  fExtractObjects(
    path, isTPCNode, [&result, this](Path_t& path) { result.push_back(makeTPC(path)); });
  return result;
}

//------------------------------------------------------------------------------
geo::TPCGeo geo::GeometryBuilderStandard::makeTPC(Path_t& path) const
{
  auto const [tpc, hash_value] = path.current_entry();
  auto tpc_matrix = path.currentTransformation<TransformationMatrix>();

  std::vector<PlaneGeo> planes;
  fExtractObjects(
    path, isPlaneNode, [&planes](Path_t& path) { planes.push_back(makePlane(path)); });

  // we don't keep the active volume information... just store its center:
  auto const* active_volume_node = TPCGeo::NodeForActiveVolume(tpc);
  assert(active_volume_node);

  auto const daughter_matrix =
    tpc_matrix * makeTransformationMatrix(*active_volume_node->GetMatrix());
  auto const active_volume_center =
    vect::toPoint(daughter_matrix(ROOT::Math::Transform3D::Point{}));

  auto const driftAxisWithSign = DriftAxisWithSign(active_volume_center, planes);
  auto const driftAxis = vect::normalize(planes[0].GetBoxCenter() - active_volume_center);
  double distance_to_last_plane{std::numeric_limits<double>::min()};
  geo::PlaneGeo const* last_plane = nullptr;
  for (auto const& plane : planes) {
    double const distance = vect::dot(plane.GetBoxCenter() - active_volume_center, driftAxis);
    if (distance > distance_to_last_plane) {
      distance_to_last_plane = distance;
      last_plane = &plane;
    }
  }

  auto const* active_volume_box =
    static_cast<TGeoBBox const*>(active_volume_node->GetVolume()->GetShape());
  auto const drift_distance = driftDistance(
    active_volume_center, driftAxisWithSign, active_volume_box, last_plane->GetBoxCenter());

  return {tpc, hash_value, std::move(tpc_matrix), driftAxisWithSign, drift_distance};
}
