////////////////////////////////////////////////////////////////////////
/// \file  WireReadoutGeom.cxx
/// \brief Interface to geometry for wire readouts
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcorealg/CoreUtils/NumericUtils.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/GeometryBuilderStandard.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/Intersections.h"
#include "larcorealg/Geometry/WireReadoutGeomBuilderStandard.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TGeoManager.h"

#include <algorithm>
#include <cassert>

namespace {
  /// Throws an exception ("GeometryCore" category) unless pid1 and pid2 are on different
  /// planes of the same TPC (ID validity is not checked)
  void CheckIndependentPlanesOnSameTPC(geo::PlaneID const& pid1,
                                       geo::PlaneID const& pid2,
                                       const char* caller)
  {
    if (pid1.asTPCID() != pid2.asTPCID()) {
      throw cet::exception("GeometryCore")
        << caller << " needs two planes on the same TPC (got " << std::string(pid1) << " and "
        << std::string(pid2) << ")\n";
    }
    if (pid1 == pid2) {
      throw cet::exception("GeometryCore")
        << caller << " needs two different planes, got " << std::string(pid1) << " twice\n";
    }
  }

  //......................................................................
  std::vector<double> Pitches(std::vector<geo::PlaneGeo> const& planes)
  {
    std::vector<double> result;
    result.reserve(planes.size());
    geo::Point_t refPlaneCenter = planes[0].GetCenter();
    for (size_t p = 0; p < planes.size(); ++p) {
      geo::Point_t const& center = planes[p].GetCenter();
      result[p] = (p == 0) ? 0.0 : result[p - 1] + std::abs(center.X() - refPlaneCenter.X());
      refPlaneCenter = center;
    }
    return result;
  }

  //......................................................................
  std::set<geo::View_t> Views(std::vector<geo::PlaneGeo> const& planes)
  {
    std::set<geo::View_t> result;
    std::transform(planes.cbegin(),
                   planes.cend(),
                   std::inserter(result, result.begin()),
                   [](auto const& plane) { return plane.View(); });
    return result;
  }
}

namespace geo {

  WireReadoutGeom::WireReadoutGeom(GeometryCore const* geom,
                                   std::unique_ptr<WireReadoutGeometryBuilder> builder,
                                   std::unique_ptr<WireReadoutSorter> sorter)
    : Iteration{details::ReadoutIterationPolicy{geom, this},
                details::ToReadoutGeometryElement{this}}
    , fGeom{geom}
  {
    GeoNodePath path{geom->ROOTGeoManager()->GetTopNode()};
    auto planesPerTPC = builder->extractPlanes(path);

    Compare<WireGeo> compareWires = [&sorter](auto const& a, auto const& b) {
      return sorter->compareWires(a, b);
    };

    // Update views
    std::set<View_t> updatedViews;
    for (auto const& tpc : fGeom->Iterate<TPCGeo>()) {
      auto& planes = fPlanes[tpc.ID()] = planesPerTPC.at(tpc.Hash());
      SortPlanes(tpc.ID(), planes);
      SortSubVolumes(planes, compareWires);
      UpdateAfterSorting(tpc, planes);

      auto const& TPCviews = ::Views(planes);
      updatedViews.insert(TPCviews.cbegin(), TPCviews.cend());
      fPlane0Pitch[tpc.ID()] = Pitches(planes);
    }
    allViews = std::move(updatedViews);
  }
  WireReadoutGeom::~WireReadoutGeom() = default;

  //......................................................................
  void WireReadoutGeom::SortPlanes(TPCID const& tpcid, std::vector<PlaneGeo>& planes) const
  {
    // Sort planes by increasing drift distance.
    //
    // This function should work in bootstrap mode, relying on least things as
    // possible. Therefore we compute here a proxy of the drift axis.

    // Determine the drift axis (or close to): from TPC center to plane center

    // Instead of using the plane center, which might be not available yet, we use the
    // plane box center, which only needs the geometry description to be available.

    // We use the first plane -- it does not make any difference.
    auto const TPCcenter = fGeom->TPC(tpcid).GetCenter();
    auto const driftAxis = vect::normalize(planes[0].GetBoxCenter() - TPCcenter);

    auto by_distance = [&TPCcenter, &driftAxis](auto const& a, auto const& b) {
      return vect::dot(a.GetBoxCenter() - TPCcenter, driftAxis) <
             vect::dot(b.GetBoxCenter() - TPCcenter, driftAxis);
    };
    cet::sort_all(planes, by_distance);
  }

  //......................................................................
  // Returns distance between plane 0 to each of the remaining planes not the distance
  // between two consecutive planes
  double WireReadoutGeom::Plane0Pitch(TPCID const& tpcid, unsigned int p) const
  {
    return fPlane0Pitch.at(tpcid)[p];
  }

  //......................................................................
  double WireReadoutGeom::PlanePitch(TPCID const& tpcid, unsigned int p1, unsigned int p2) const
  {
    auto const& plane0_pitch = fPlane0Pitch.at(tpcid);
    return std::abs(plane0_pitch[p2] - plane0_pitch[p1]);
  }

  //......................................................................
  // Sort the PlaneGeo objects and the WireGeo objects inside
  void WireReadoutGeom::SortSubVolumes(std::vector<PlaneGeo>& planes,
                                       Compare<WireGeo> compareWires) const
  {
    for (PlaneGeo& plane : planes)
      plane.SortWires(compareWires);
  }

  //......................................................................
  void WireReadoutGeom::UpdateAfterSorting(TPCGeo const& tpc, std::vector<PlaneGeo>& planes)
  {
    // ask the planes to update; also check
    for (unsigned int p = 0; p < planes.size(); ++p) {
      planes[p].UpdateAfterSorting(PlaneID(tpc.ID(), p), tpc);

      // check that the plane normal is opposite to the TPC drift direction
      assert(lar::util::makeVector3DComparison(1e-5).equal(-(planes[p].GetNormalDirection()),
                                                           tpc.DriftDir()));
    }
  }

  //......................................................................
  // Number of different views, or wire orientations
  unsigned int WireReadoutGeom::Nviews() const { return MaxPlanes(); }

  //......................................................................
  unsigned int WireReadoutGeom::MaxPlanes() const
  {
    unsigned int maxPlanes = 0;
    for (auto const& [tpc_id, planes] : fPlanes) {
      unsigned int maxPlanesInTPC = planes.size();
      if (maxPlanesInTPC > maxPlanes) maxPlanes = maxPlanesInTPC;
    }
    return maxPlanes;
  }

  //......................................................................
  unsigned int WireReadoutGeom::MaxWires() const
  {
    unsigned int maxWires = 0;
    for (PlaneGeo const& plane : Iterate<PlaneGeo>()) {
      unsigned int maxWiresInPlane = plane.Nwires();
      if (maxWiresInPlane > maxWires) maxWires = maxWiresInPlane;
    }
    return maxWires;
  }

  //......................................................................
  bool WireReadoutGeom::HasPlane(PlaneID const& planeid) const
  {
    return PlanePtr(planeid) != nullptr;
  }

  //......................................................................
  unsigned int WireReadoutGeom::Nplanes(TPCID const& tpcid) const
  {
    auto it = fPlanes.find(tpcid);
    if (it == fPlanes.cend()) { return 0; }
    return it->second.size();
  }

  //......................................................................
  PlaneGeo const& WireReadoutGeom::FirstPlane(TPCID const& tpcid) const
  {
    if (auto const* plane_ptr = PlanePtr({tpcid, 0})) { return *plane_ptr; }
    throw cet::exception("WireReadoutGeom")
      << "TPC with ID " << tpcid << " does not have a first plane.";
  }

  //......................................................................
  PlaneGeo const& WireReadoutGeom::Plane(TPCID const& tpcid, View_t const view) const
  {
    auto it = fPlanes.find(tpcid);
    if (it == fPlanes.cend()) {
      throw cet::exception("WireReadoutGeom") << "No TPC with ID " << tpcid << " supported.";
    }
    for (PlaneGeo const& plane : it->second) {
      if (plane.View() == view) { return plane; }
    }
    throw cet::exception("WireReadoutGeom")
      << "TPCGeo " << tpcid << " has no plane for view #" << (size_t)view << "\n";
  }

  //......................................................................
  PlaneGeo const& WireReadoutGeom::Plane(PlaneID const& planeid) const
  {
    if (auto const* plane_ptr = PlanePtr(planeid)) { return *plane_ptr; }
    throw cet::exception("WireReadoutGeom") << "Plane with ID " << planeid << " is not supported.";
  }

  //......................................................................
  PlaneGeo const* WireReadoutGeom::PlanePtr(PlaneID const& planeid) const
  {
    auto const [tpc_id, plane] = std::make_pair(planeid.parentID(), planeid.Plane);
    auto it = fPlanes.find(tpc_id);
    if (it == fPlanes.cend()) { return nullptr; }
    if (std::size_t const n = it->second.size(); n <= plane) { return nullptr; }
    return &it->second[plane];
  }

  //......................................................................
  std::set<View_t> WireReadoutGeom::Views(TPCID const& tpcid) const
  {
    auto it = fPlanes.find(tpcid);
    if (it == fPlanes.cend()) { return {}; }
    return ::Views(it->second);
  }

  //......................................................................
  unsigned int WireReadoutGeom::Nwires(PlaneID const& planeid) const
  {
    PlaneGeo const* pPlane = PlanePtr(planeid);
    return pPlane ? pPlane->NElements() : 0;
  }

  //......................................................................
  bool WireReadoutGeom::HasWire(WireID const& wireid) const
  {
    PlaneGeo const* pPlane = PlanePtr(wireid);
    return pPlane ? pPlane->HasWire(wireid) : false;
  }

  //......................................................................
  WireGeo const* WireReadoutGeom::WirePtr(WireID const& wireid) const
  {
    PlaneGeo const* pPlane = PlanePtr(wireid);
    return pPlane ? pPlane->WirePtr(wireid) : nullptr;
  }

  //......................................................................
  WireGeo const& WireReadoutGeom::Wire(WireID const& wireid) const
  {
    return Plane(wireid).Wire(wireid);
  }

  //......................................................................
  void WireReadoutGeom::WireEndPoints(WireID const& wireid, double* xyzStart, double* xyzEnd) const
  {
    auto const [start, end] = WireEndPoints(wireid);

    xyzStart[0] = start.X();
    xyzStart[1] = start.Y();
    xyzStart[2] = start.Z();
    xyzEnd[0] = end.X();
    xyzEnd[1] = end.Y();
    xyzEnd[2] = end.Z();

    if (xyzEnd[2] < xyzStart[2]) {
      // Ensure that "End" has higher z-value than "Start"
      std::swap(xyzStart[0], xyzEnd[0]);
      std::swap(xyzStart[1], xyzEnd[1]);
      std::swap(xyzStart[2], xyzEnd[2]);
    }
    if (xyzEnd[1] < xyzStart[1] && std::abs(xyzEnd[2] - xyzStart[2]) < 0.01) {
      // If wire is vertical ensure that "End" has higher y-value than "Start"
      std::swap(xyzStart[0], xyzEnd[0]);
      std::swap(xyzStart[1], xyzEnd[1]);
      std::swap(xyzStart[2], xyzEnd[2]);
    }
  }

  //......................................................................
  std::optional<WireIDIntersection> WireReadoutGeom::WireIDsIntersect(WireID const& wid1,
                                                                      WireID const& wid2) const
  {
    if (!WireIDIntersectionCheck(wid1, wid2)) { return std::nullopt; }

    // get the endpoints to see if wires intersect
    Segment const w1 = WireEndPoints(wid1);
    Segment const w2 = WireEndPoints(wid2);

    // TODO extract the coordinates in the right way;
    // is it any worth, since then the result is in (y, z), whatever it means?
    WireIDIntersection result;
    bool const cross = IntersectLines(start(w1).Y(),
                                      start(w1).Z(),
                                      finish(w1).Y(),
                                      finish(w1).Z(),
                                      start(w2).Y(),
                                      start(w2).Z(),
                                      finish(w2).Y(),
                                      finish(w2).Z(),
                                      result.y,
                                      result.z);
    if (!cross) { return std::nullopt; }

    bool const within = lar::util::PointWithinSegments(start(w1).Y(),
                                                       start(w1).Z(),
                                                       finish(w1).Y(),
                                                       finish(w1).Z(),
                                                       start(w2).Y(),
                                                       start(w2).Z(),
                                                       finish(w2).Y(),
                                                       finish(w2).Z(),
                                                       result.y,
                                                       result.z);

    result.TPC = (within ? wid1.TPC : TPCID::InvalidID);
    return result;
  }

  //......................................................................
  bool WireReadoutGeom::WireIDsIntersect(WireID const& wid1,
                                         WireID const& wid2,
                                         Point_t& intersection) const
  {
    // This is not a real 3D intersection: the wires do not cross, since they are required
    // to belong to two different planes.
    //
    // We take the point on the first wire which is closest to the other one.
    static_assert(std::numeric_limits<decltype(intersection.X())>::has_infinity,
                  "the vector coordinate type can't represent infinity!");
    constexpr auto infinity = std::numeric_limits<decltype(intersection.X())>::infinity();

    if (!WireIDIntersectionCheck(wid1, wid2)) {
      intersection = {infinity, infinity, infinity};
      return false;
    }

    WireGeo const& wire1 = Wire(wid1);
    WireGeo const& wire2 = Wire(wid2);

    // Distance of the intersection point from the center of the two wires:
    IntersectionPointAndOffsets<Point_t> intersectionAndOffset =
      WiresIntersectionAndOffsets(wire1, wire2);
    intersection = intersectionAndOffset.point;

    return std::abs(intersectionAndOffset.offset1) <= wire1.HalfL() &&
           std::abs(intersectionAndOffset.offset2) <= wire2.HalfL();
  }

  //......................................................................
  // This method returns the distance between wires in the specified view it assumes all
  // planes of a given view have the same pitch
  double WireReadoutGeom::WireAngleToVertical(View_t view, TPCID const& tpcid) const
  {
    for (PlaneGeo const& plane : Iterate<PlaneGeo>(tpcid)) {
      if (plane.View() == view) return plane.ThetaZ();
    } // for
    throw cet::exception("WireReadoutGeom")
      << "WireAngleToVertical(): no view \"" << PlaneGeo::ViewName(view) << "\" (#" << ((int)view)
      << ") in " << std::string(tpcid);
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutGeom::NOpChannels(unsigned int NOpDets) const
  {
    // By default just return the number of optical detectos
    return NOpDets;
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutGeom::MaxOpChannel(unsigned int NOpDets) const
  {
    // By default just return the number of optical detectos
    return NOpChannels(NOpDets);
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutGeom::NOpHardwareChannels(unsigned int /*opDet*/) const
  {
    // By default, 1 channel per optical detector
    return 1;
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutGeom::OpChannel(unsigned int detNum, unsigned int /* channel */) const
  {
    return detNum;
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutGeom::OpDetFromOpChannel(unsigned int opChannel) const
  {
    return opChannel;
  }

  //----------------------------------------------------------------------------
  OpDetGeo const& WireReadoutGeom::OpDetGeoFromOpChannel(unsigned int opChannel) const
  {
    return fGeom->OpDetGeoFromOpDet(OpDetFromOpChannel(opChannel));
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutGeom::HardwareChannelFromOpChannel(unsigned int /* opChannel */) const
  {
    return 0;
  }

  //----------------------------------------------------------------------------
  bool WireReadoutGeom::IsValidOpChannel(unsigned int opChannel, unsigned int NOpDets) const
  {
    // Check channel number
    if (opChannel >= NOpChannels(NOpDets)) return false;

    // Check opdet number
    unsigned int opdet = OpDetFromOpChannel(opChannel);
    if (opdet >= NOpDets) return false;

    // Check hardware channel number
    unsigned int hChan = HardwareChannelFromOpChannel(opChannel);
    return hChan < NOpHardwareChannels(opdet);
  }

  SigType_t WireReadoutGeom::SignalType(raw::ChannelID_t const channel) const
  {
    return SignalTypeForChannelImpl(channel);
  }

  SigType_t WireReadoutGeom::SignalType(readout::ROPID const& ropid) const
  {
    return SignalTypeForROPIDImpl(ropid);
  }

  SigType_t WireReadoutGeom::SignalTypeForROPIDImpl(readout::ROPID const& ropid) const
  {
    return SignalType(FirstChannelInROP(ropid));
  }

  //......................................................................
  std::vector<raw::ChannelID_t> WireReadoutGeom::ChannelsInTPCs() const
  {
    std::vector<raw::ChannelID_t> channels;
    channels.reserve(Nchannels());

    for (auto const& ts : Iterate<readout::TPCsetID>()) {
      for (auto const t : TPCsetToTPCs(ts)) {
        for (auto const& wire : Iterate<WireID>(t)) {
          channels.push_back(PlaneWireToChannel(wire));
        }
      }
    }
    std::sort(channels.begin(), channels.end());
    auto last = std::unique(channels.begin(), channels.end());
    channels.erase(last, channels.end());
    return channels;
  }

  //......................................................................
  unsigned int WireReadoutGeom::NOpChannels() const { return NOpChannels(fGeom->NOpDets()); }

  //......................................................................
  unsigned int WireReadoutGeom::MaxOpChannel() const { return MaxOpChannel(fGeom->NOpDets()); }

  //......................................................................
  bool WireReadoutGeom::IsValidOpChannel(int opChannel) const
  {
    return IsValidOpChannel(opChannel, fGeom->NOpDets());
  }

  //......................................................................
  SigType_t WireReadoutGeom::SignalType(PlaneID const& pid) const
  {
    // map wire plane -> readout plane -> first channel, then use SignalType(channel)

    auto const ropid = WirePlaneToROP(pid);
    if (!ropid.isValid) {
      throw cet::exception("WireReadoutGeom") << "SignalType(): Mapping of wire plane "
                                              << std::string(pid) << " to readout plane failed!\n";
    }
    return SignalType(ropid);
  }

  //......................................................................
  View_t WireReadoutGeom::View(readout::ROPID const& ropid) const
  {
    auto const pid = FirstWirePlaneInROP(ropid);
    return pid ? Plane(pid).View() : kUnknown;
  }

  View_t WireReadoutGeom::View(raw::ChannelID_t const channel) const
  {
    return (channel == raw::InvalidChannelID) ? kUnknown : View(ChannelToROP(channel));
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t WireReadoutGeom::NearestChannel(Point_t const& worldPos,
                                                   PlaneID const& planeid) const
  {
    // This method is supposed to return a channel number rather than a wire number.
    // Perform the conversion here (although, maybe faster if we deal in wire numbers
    // rather than channel numbers?)

    // NOTE on failure both NearestChannel() and upstream:
    // * according to documentation, should return invalid channel
    // * in the actual code throw an exception because of a BUG
    //
    // The following implementation automatically becomes in fact compliant to the
    // documentation if upstreams becomes compliant to.  When that happens, just delete
    // this comment.
    WireID const wireID = Plane(planeid).NearestWireID(worldPos);
    return wireID ? PlaneWireToChannel(wireID) : raw::InvalidChannelID;
  }

  //......................................................................
  std::optional<WireIDIntersection> WireReadoutGeom::ChannelsIntersect(
    raw::ChannelID_t const c1,
    raw::ChannelID_t const c2) const
  {
    // [GP] these errors should be exceptions, and this function is deprecated because it
    // violates interoperability
    std::vector<WireID> chan1wires = ChannelToWire(c1);
    if (chan1wires.empty()) {
      mf::LogError("ChannelsIntersect")
        << "1st channel " << c1 << " maps to no wire (is it a real one?)";
      return std::nullopt;
    }
    std::vector<WireID> chan2wires = ChannelToWire(c2);
    if (chan2wires.empty()) {
      mf::LogError("ChannelsIntersect")
        << "2nd channel " << c2 << " maps to no wire (is it a real one?)";
      return std::nullopt;
    }

    if (chan1wires.size() > 1) {
      mf::LogWarning("ChannelsIntersect")
        << "1st channel " << c1 << " maps to " << chan2wires.size() << " wires; using the first!";
      return std::nullopt;
    }
    if (chan2wires.size() > 1) {
      mf::LogError("ChannelsIntersect")
        << "2nd channel " << c2 << " maps to " << chan2wires.size() << " wires; using the first!";
      return std::nullopt;
    }

    return WireIDsIntersect(chan1wires[0], chan2wires[0]);
  }

  readout::TPCsetID WireReadoutGeom::FindTPCsetAtPosition(Point_t const& worldLoc) const
  {
    return TPCtoTPCset(fGeom->FindTPCAtPosition(worldLoc));
  }

  //--------------------------------------------------------------------
  bool WireReadoutGeom::WireIDIntersectionCheck(WireID const& wid1, WireID const& wid2) const
  {
    if (wid1.asTPCID() != wid2) {
      mf::LogError("WireIDIntersectionCheck")
        << "Comparing two wires on different TPCs: return failure.";
      return false;
    }
    if (wid1.Plane == wid2.Plane) {
      mf::LogError("WireIDIntersectionCheck")
        << "Comparing two wires in the same plane: return failure";
      return false;
    }
    if (!HasWire(wid1)) {
      mf::LogError("WireIDIntersectionCheck")
        << "1st wire " << wid1 << " does not exist (max wire number: " << Nwires(wid1.planeID())
        << ")";
      return false;
    }
    if (!HasWire(wid2)) {
      mf::LogError("WireIDIntersectionCheck")
        << "2nd wire " << wid2 << " does not exist (max wire number: " << Nwires(wid2.planeID())
        << ")";
      return false;
    }
    return true;
  }

  //----------------------------------------------------------------------------
  PlaneID WireReadoutGeom::ThirdPlane(PlaneID const& pid1, PlaneID const& pid2) const
  {
    // how many planes in the TPC pid1 belongs to:
    unsigned const int nPlanes = Nplanes(pid1);
    if (nPlanes != 3) {
      throw cet::exception("GeometryCore")
        << "ThirdPlane() supports only TPCs with 3 planes, and I see " << nPlanes << " instead\n";
    }

    PlaneID::PlaneID_t target_plane = nPlanes;
    for (PlaneID::PlaneID_t iPlane = 0; iPlane < nPlanes; ++iPlane) {
      if ((iPlane == pid1.Plane) || (iPlane == pid2.Plane)) continue;
      if (target_plane != nPlanes) {
        throw cet::exception("GeometryCore")
          << "ThirdPlane() found too many planes that are not " << std::string(pid1) << " nor "
          << std::string(pid2) << "! (first " << target_plane << ", then " << iPlane << ")\n";
      } // if we had a target already
      target_plane = iPlane;
    } // for
    if (target_plane == nPlanes) {
      throw cet::exception("GeometryCore")
        << "ThirdPlane() can't find a plane that is not " << std::string(pid1) << " nor "
        << std::string(pid2) << "!\n";
    }

    return PlaneID(pid1, target_plane);
  }

  //----------------------------------------------------------------------------
  double WireReadoutGeom::ThirdPlaneSlope(PlaneID const& pid1,
                                          double slope1,
                                          PlaneID const& pid2,
                                          double slope2,
                                          PlaneID const& output_plane) const
  {
    CheckIndependentPlanesOnSameTPC(pid1, pid2, "ThirdPlaneSlope()");

    // We need the "wire coordinate direction" for each plane.  This is perpendicular to
    // the wire orientation.  PlaneGeo::PhiZ() defines the right orientation too.
    return ComputeThirdPlaneSlope(
      Plane(pid1).PhiZ(), slope1, Plane(pid2).PhiZ(), slope2, Plane(output_plane).PhiZ());
  }

  //----------------------------------------------------------------------------
  double WireReadoutGeom::ThirdPlaneSlope(PlaneID const& pid1,
                                          double slope1,
                                          PlaneID const& pid2,
                                          double slope2) const
  {
    return ThirdPlaneSlope(pid1, slope1, pid2, slope2, ThirdPlane(pid1, pid2));
  }

  //----------------------------------------------------------------------------
  double WireReadoutGeom::ThirdPlane_dTdW(PlaneID const& pid1,
                                          double slope1,
                                          PlaneID const& pid2,
                                          double slope2,
                                          PlaneID const& output_plane) const
  {
    CheckIndependentPlanesOnSameTPC(pid1, pid2, "ThirdPlane_dTdW()");

    double angle[3], pitch[3];
    PlaneGeo const* const planes[3] = {&Plane(pid1), &Plane(pid2), &Plane(output_plane)};

    // We need wire pitch and "wire coordinate direction" for each plane.  The latter is
    // perpendicular to the wire orientation.  PlaneGeo::PhiZ() defines the right
    // orientation too.
    for (size_t i = 0; i < 3; ++i) {
      angle[i] = planes[i]->PhiZ();
      pitch[i] = planes[i]->WirePitch();
    }

    return ComputeThirdPlane_dTdW(
      angle[0], pitch[0], slope1, angle[1], pitch[1], slope2, angle[2], pitch[2]);
  }

  //----------------------------------------------------------------------------
  double WireReadoutGeom::ThirdPlane_dTdW(PlaneID const& pid1,
                                          double slope1,
                                          PlaneID const& pid2,
                                          double slope2) const
  {
    return ThirdPlane_dTdW(pid1, slope1, pid2, slope2, ThirdPlane(pid1, pid2));
  }

  //----------------------------------------------------------------------------
  // Given slopes dTime/dWire in two planes, return with the slope in the 3rd plane.
  // Requires slopes to be in the same metrics, e.g. converted in a distances ratio.

  // Note: Uses equation in H. Greenlee's talk:
  //       https://cdcvs.fnal.gov/redmine/attachments/download/1821/larsoft_apr20_2011.pdf
  //       slide 2
  double WireReadoutGeom::ComputeThirdPlaneSlope(double angle1,
                                                 double slope1,
                                                 double angle2,
                                                 double slope2,
                                                 double angle3)
  {
    // note that, if needed, the trigonometric functions can be pre-calculated.

    // Can't resolve very small slopes
    if ((std::abs(slope1) < 0.001) && (std::abs(slope2)) < 0.001) return 0.001;

    // We need the "wire coordinate direction" for each plane.  This is perpendicular to
    // the wire orientation.
    double slope3 = 0.001;
    if (std::abs(slope1) > 0.001 && std::abs(slope2) > 0.001) {
      slope3 =
        (+(1. / slope1) * std::sin(angle3 - angle2) - (1. / slope2) * std::sin(angle3 - angle1)) /
        std::sin(angle1 - angle2);
    }
    if (slope3 != 0.)
      slope3 = 1. / slope3;
    else
      slope3 = 999.;

    return slope3;
  }

  //----------------------------------------------------------------------------
  double WireReadoutGeom::ComputeThirdPlane_dTdW(double angle1,
                                                 double pitch1,
                                                 double dTdW1,
                                                 double angle2,
                                                 double pitch2,
                                                 double dTdW2,
                                                 double angle_target,
                                                 double pitch_target)
  {
    // we need to convert dt/dw into homogeneous coordinates, and then back;
    // slope = [dT * (TDCperiod / driftVelocity)] / [dW * wirePitch]
    // The coefficient of dT is assumed to be the same for all the planes, and it finally
    // cancels out. Pitches cancel out only if they are all the same.
    return pitch_target *
           ComputeThirdPlaneSlope(angle1, dTdW1 / pitch1, angle2, dTdW2 / pitch2, angle_target);
  }

}
