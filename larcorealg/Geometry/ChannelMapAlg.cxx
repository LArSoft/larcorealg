////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapAlg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace geo {

  ChannelMapAlg::ChannelMapAlg(GeometryCore const* geom)
    : Iteration{details::ReadoutIterationPolicy{geom, this}}, fGeom{geom}
  {}
  ChannelMapAlg::~ChannelMapAlg() = default;

  //----------------------------------------------------------------------------
  unsigned int ChannelMapAlg::NOpChannels(unsigned int NOpDets) const
  {
    // By default just return the number of optical detectos
    return NOpDets;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapAlg::MaxOpChannel(unsigned int NOpDets) const
  {
    // By default just return the number of optical detectos
    return NOpChannels(NOpDets);
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapAlg::NOpHardwareChannels(unsigned int /*opDet*/) const
  {
    // By default, 1 channel per optical detector
    return 1;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapAlg::OpChannel(unsigned int detNum, unsigned int /* channel */) const
  {
    return detNum;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapAlg::OpDetFromOpChannel(unsigned int opChannel) const { return opChannel; }

  //----------------------------------------------------------------------------
  OpDetGeo const& ChannelMapAlg::OpDetGeoFromOpChannel(unsigned int opChannel) const
  {
    return fGeom->OpDetGeoFromOpDet(OpDetFromOpChannel(opChannel));
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapAlg::HardwareChannelFromOpChannel(unsigned int /* opChannel */) const
  {
    return 0;
  }

  //----------------------------------------------------------------------------
  bool ChannelMapAlg::IsValidOpChannel(unsigned int opChannel, unsigned int NOpDets) const
  {
    // Check channel number
    if (opChannel >= NOpChannels(NOpDets)) return false;

    // Check opdet number
    unsigned int opdet = OpDetFromOpChannel(opChannel);
    if (opdet >= NOpDets) return false;

    // Check hardware channel number
    unsigned int hChan = HardwareChannelFromOpChannel(opChannel);
    if (hChan >= NOpHardwareChannels(opdet)) return false;

    return true;
  }

  //----------------------------------------------------------------------------
  size_t ChannelMapAlg::NearestAuxDet(Point_t const& point,
                                      std::vector<geo::AuxDetGeo> const& auxDets,
                                      double tolerance) const
  {
    for (size_t a = 0; a < auxDets.size(); ++a) {
      auto const localPoint = auxDets[a].toLocalCoords(point);

      double const HalfCenterWidth = 0.5 * (auxDets[a].HalfWidth1() + auxDets[a].HalfWidth2());

      if (localPoint.Z() >= -(auxDets[a].Length() / 2 + tolerance) &&
          localPoint.Z() <= (auxDets[a].Length() / 2 + tolerance) &&
          localPoint.Y() >= -auxDets[a].HalfHeight() - tolerance &&
          localPoint.Y() <= auxDets[a].HalfHeight() + tolerance &&
          // if AuxDet a is a box, then HalfSmallWidth = HalfWidth
          localPoint.X() >= -HalfCenterWidth +
                              localPoint.Z() * (HalfCenterWidth - auxDets[a].HalfWidth2()) /
                                (0.5 * auxDets[a].Length()) -
                              tolerance &&
          localPoint.X() <= HalfCenterWidth -
                              localPoint.Z() * (HalfCenterWidth - auxDets[a].HalfWidth2()) /
                                (0.5 * auxDets[a].Length()) +
                              tolerance)
        return a;

    } // for loop over AudDet a

    // throw an exception because we couldn't find the sensitive volume
    throw cet::exception("ChannelMap") << "Can't find AuxDet for position (" << point.X() << ","
                                       << point.Y() << "," << point.Z() << ")\n";
  }

  //----------------------------------------------------------------------------
  size_t ChannelMapAlg::NearestSensitiveAuxDet(Point_t const& point,
                                               std::vector<geo::AuxDetGeo> const& auxDets,
                                               double tolerance) const
  {
    size_t auxDetIdx = NearestAuxDet(point, auxDets, tolerance);
    geo::AuxDetGeo const& adg = auxDets[auxDetIdx];

    for (size_t a = 0; a < adg.NSensitiveVolume(); ++a) {
      geo::AuxDetSensitiveGeo const& adsg = adg.SensitiveVolume(a);
      auto const localPoint = adsg.toLocalCoords(point);

      double const HalfCenterWidth = 0.5 * (adsg.HalfWidth1() + adsg.HalfWidth2());

      if (localPoint.Z() >= -(adsg.Length() / 2 + tolerance) &&
          localPoint.Z() <= (adsg.Length() / 2 + tolerance) &&
          localPoint.Y() >= -adsg.HalfHeight() - tolerance &&
          localPoint.Y() <= adsg.HalfHeight() + tolerance &&
          // if AuxDet a is a box, then HalfSmallWidth = HalfWidth
          localPoint.X() >=
            -HalfCenterWidth +
              localPoint.Z() * (HalfCenterWidth - adsg.HalfWidth2()) / (0.5 * adsg.Length()) -
              tolerance &&
          localPoint.X() <=
            HalfCenterWidth -
              localPoint.Z() * (HalfCenterWidth - adsg.HalfWidth2()) / (0.5 * adsg.Length()) +
              tolerance)
        return a;
    } // for loop over AuxDetSensitive a

    // throw an exception because we couldn't find the sensitive volume
    throw cet::exception("Geometry") << "Can't find AuxDetSensitive for position (" << point.X()
                                     << "," << point.Y() << "," << point.Z() << ")\n";
  }

  //----------------------------------------------------------------------------
  size_t ChannelMapAlg::ChannelToAuxDet(std::vector<geo::AuxDetGeo> const& /* auxDets */,
                                        std::string const& detName,
                                        uint32_t const& /*channel*/) const
  {
    // loop over the map of AuxDet names to Geo object numbers to determine which auxdet
    // we have.  If no name in the map matches the provided string, throw an exception
    for (auto itr : fADNameToGeo)
      if (itr.first.compare(detName) == 0) return itr.second;

    throw cet::exception("Geometry") << "No AuxDetGeo matching name: " << detName;
  }

  //----------------------------------------------------------------------------
  // the first member of the pair is the index in the auxDets vector for the AuxDetGeo,
  // the second member is the index in the vector of AuxDetSensitiveGeos for that AuxDetGeo
  std::pair<size_t, size_t> ChannelMapAlg::ChannelToSensitiveAuxDet(
    std::vector<geo::AuxDetGeo> const& auxDets,
    std::string const& detName,
    uint32_t const& channel) const
  {
    size_t adGeoIdx = ChannelToAuxDet(auxDets, detName, channel);

    // look for the index of the sensitive volume for the given channel
    if (fADChannelToSensitiveGeo.count(adGeoIdx) > 0) {

      auto itr = fADChannelToSensitiveGeo.find(adGeoIdx);

      // get the vector of channels to AuxDetSensitiveGeo index
      if (channel < itr->second.size()) return std::make_pair(adGeoIdx, itr->second[channel]);

      throw cet::exception("Geometry")
        << "Given AuxDetSensitive channel, " << channel
        << ", cannot be found in vector associated to AuxDetGeo index: " << adGeoIdx
        << ". Vector has size " << itr->second.size();
    }

    throw cet::exception("Geometry") << "Given AuxDetGeo with index " << adGeoIdx
                                     << " does not correspond to any vector of sensitive volumes";
  }

  geo::SigType_t ChannelMapAlg::SignalType(raw::ChannelID_t const channel) const
  {
    return SignalTypeForChannelImpl(channel);
  }

  geo::SigType_t ChannelMapAlg::SignalType(readout::ROPID const& ropid) const
  {
    return SignalTypeForROPIDImpl(ropid);
  }

  geo::SigType_t ChannelMapAlg::SignalTypeForROPIDImpl(readout::ROPID const& ropid) const
  {
    return SignalType(FirstChannelInROP(ropid));
  }

  //......................................................................
  std::vector<raw::ChannelID_t> ChannelMapAlg::ChannelsInTPCs() const
  {
    std::vector<raw::ChannelID_t> channels;
    channels.reserve(Nchannels());

    for (auto const& ts : Iterate<readout::TPCsetID>()) {
      for (auto const t : TPCsetToTPCs(ts)) {
        for (auto const& wire : fGeom->Iterate<WireID>(t)) {
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
  unsigned int ChannelMapAlg::NOpChannels() const { return NOpChannels(fGeom->NOpDets()); }

  //......................................................................
  unsigned int ChannelMapAlg::MaxOpChannel() const { return MaxOpChannel(fGeom->NOpDets()); }

  //......................................................................
  bool ChannelMapAlg::IsValidOpChannel(int opChannel) const
  {
    return IsValidOpChannel(opChannel, fGeom->NOpDets());
  }

  //......................................................................
  SigType_t ChannelMapAlg::SignalType(PlaneID const& pid) const
  {
    // map wire plane -> readout plane -> first channel, then use SignalType(channel)

    auto const ropid = WirePlaneToROP(pid);
    if (!ropid.isValid) {
      throw cet::exception("ChannelMapAlg") << "SignalType(): Mapping of wire plane "
                                            << std::string(pid) << " to readout plane failed!\n";
    }
    return SignalType(ropid);
  }

  //......................................................................
  View_t ChannelMapAlg::View(readout::ROPID const& ropid) const
  {
    auto const pid = FirstWirePlaneInROP(ropid);
    return pid ? fGeom->Plane(pid).View() : kUnknown;
  }

  View_t ChannelMapAlg::View(raw::ChannelID_t const channel) const
  {
    return (channel == raw::InvalidChannelID) ? kUnknown : View(ChannelToROP(channel));
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t ChannelMapAlg::NearestChannel(Point_t const& worldPos,
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
    WireID const wireID = fGeom->Plane(planeid).NearestWireID(worldPos);
    return wireID ? PlaneWireToChannel(wireID) : raw::InvalidChannelID;
  }

  //......................................................................
  unsigned int ChannelMapAlg::FindAuxDetAtPosition(Point_t const& point, double tolerance) const
  {
    return NearestAuxDet(point, fGeom->AuxDets(), tolerance);
  }

  //......................................................................
  AuxDetGeo const& ChannelMapAlg::PositionToAuxDet(Point_t const& point,
                                                   unsigned int& ad,
                                                   double tolerance) const
  {
    // FIXME (KJK): Possibly remove this function?
    // locate the desired Auxiliary Detector
    ad = FindAuxDetAtPosition(point, tolerance);
    return fGeom->AuxDet(ad);
  }

  //......................................................................
  void ChannelMapAlg::FindAuxDetSensitiveAtPosition(Point_t const& point,
                                                    std::size_t& adg,
                                                    std::size_t& sv,
                                                    double tolerance) const
  {
    adg = FindAuxDetAtPosition(point, tolerance);
    sv = NearestSensitiveAuxDet(point, fGeom->AuxDets(), tolerance);
  }

  //......................................................................
  AuxDetSensitiveGeo const& ChannelMapAlg::PositionToAuxDetSensitive(Point_t const& point,
                                                                     size_t& ad,
                                                                     size_t& sv,
                                                                     double tolerance) const
  {
    // locate the desired Auxiliary Detector
    FindAuxDetSensitiveAtPosition(point, ad, sv, tolerance);
    return fGeom->AuxDet(ad).SensitiveVolume(sv);
  }

  //......................................................................
  AuxDetGeo const& ChannelMapAlg::ChannelToAuxDet(std::string const& auxDetName,
                                                  uint32_t const channel) const
  {
    size_t adIdx = ChannelToAuxDet(fGeom->AuxDets(), auxDetName, channel);
    return fGeom->AuxDet(adIdx);
  }

  //......................................................................
  AuxDetSensitiveGeo const& ChannelMapAlg::ChannelToAuxDetSensitive(std::string const& auxDetName,
                                                                    uint32_t const channel) const
  {
    auto idx = ChannelToSensitiveAuxDet(fGeom->AuxDets(), auxDetName, channel);
    return fGeom->AuxDet(idx.first).SensitiveVolume(idx.second);
  }

  //......................................................................
  std::optional<WireIDIntersection> ChannelMapAlg::ChannelsIntersect(
    raw::ChannelID_t const c1,
    raw::ChannelID_t const c2) const
  {
    // [GP] these errors should be exceptions, and this function is deprecated
    // because it violates interoperability
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

    return fGeom->WireIDsIntersect(chan1wires[0], chan2wires[0]);
  }

  readout::TPCsetID ChannelMapAlg::FindTPCsetAtPosition(Point_t const& worldLoc) const
  {
    return TPCtoTPCset(fGeom->FindTPCAtPosition(worldLoc));
  }

}
