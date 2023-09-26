////////////////////////////////////////////////////////////////////////
/// \file  WireReadoutStandardGeom.cxx
/// \brief Interface to algorithm class for the standar, simplest detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/WireReadoutStandardGeom.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/WireReadoutGeomBuilderStandard.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
  //----------------------------------------------------------------------------
  readout::TPCsetID ConvertTPCtoTPCset(geo::TPCID const& tpcid)
  {
    if (!tpcid.isValid) return {}; // invalid ID, default-constructed
    return {(readout::CryostatID::CryostatID_t)tpcid.Cryostat,
            (readout::TPCsetID::TPCsetID_t)tpcid.TPC};
  }

  //----------------------------------------------------------------------------
  geo::TPCID ConvertTPCsetToTPC(readout::TPCsetID const& tpcsetid)
  {
    if (!tpcsetid.isValid) return {};
    return {(geo::CryostatID::CryostatID_t)tpcsetid.Cryostat, (geo::TPCID::TPCID_t)tpcsetid.TPCset};
  }

  //----------------------------------------------------------------------------
  readout::ROPID ConvertWirePlaneToROP(geo::PlaneID const& planeid)
  {
    if (!planeid.isValid) return {}; // invalid ID, default-constructed
    return {(readout::CryostatID::CryostatID_t)planeid.Cryostat,
            (readout::TPCsetID::TPCsetID_t)planeid.TPC,
            (readout::ROPID::ROPID_t)planeid.Plane};
  }

  //----------------------------------------------------------------------------
  geo::PlaneID ConvertROPtoWirePlane(readout::ROPID const& ropid)
  {
    if (!ropid.isValid) return {};
    return {(geo::CryostatID::CryostatID_t)ropid.Cryostat,
            (geo::TPCID::TPCID_t)ropid.TPCset,
            (geo::PlaneID::PlaneID_t)ropid.ROP};
  }

}

namespace geo {

  //----------------------------------------------------------------------------
  WireReadoutStandardGeom::WireReadoutStandardGeom(fhicl::ParameterSet const& pset,
                                                   GeometryCore const* geom,
                                                   std::unique_ptr<WireReadoutSorter> sorter)
    : WireReadoutGeom{geom,
                      std::make_unique<WireReadoutGeomBuilderStandard>(
                        pset.get<fhicl::ParameterSet>("Builder", {})),
                      std::move(sorter)}
  {
    fNcryostat = geom->Ncryostats();

    mf::LogInfo("WireReadoutStandardGeom") << "Initializing standard wire-readout geometry...";

    fNTPC.resize(fNcryostat);
    fWireCounts.resize(fNcryostat);
    fNPlanes.resize(fNcryostat);
    fFirstWireProj.resize(fNcryostat);
    fOrthVectorsY.resize(fNcryostat);
    fOrthVectorsZ.resize(fNcryostat);
    fPlaneBaselines.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    fPlaneIDs.clear();
    fTopChannel = 0;

    int RunningTotal = 0;

    for (auto const& cryo : geom->Iterate<CryostatGeo>()) {
      auto const cs = cryo.ID().Cryostat;
      fNTPC[cs] = cryo.NTPC();

      // Size up all the vectors
      fWireCounts[cs].resize(fNTPC[cs]);
      fFirstWireProj[cs].resize(fNTPC[cs]);
      fOrthVectorsY[cs].resize(fNTPC[cs]);
      fOrthVectorsZ[cs].resize(fNTPC[cs]);
      fPlaneBaselines[cs].resize(fNTPC[cs]);
      fWiresPerPlane[cs].resize(fNTPC[cs]);
      fNPlanes[cs].resize(fNTPC[cs]);
      fFirstChannelInThisPlane[cs].resize(fNTPC[cs]);
      fFirstChannelInNextPlane[cs].resize(fNTPC[cs]);

      for (auto const& TPC : geom->Iterate<TPCGeo>(cryo.ID())) {
        auto const TPCCount = TPC.ID().TPC;
        unsigned int PlanesThisTPC = Nplanes(TPC.ID());
        fWireCounts[cs][TPCCount].resize(PlanesThisTPC);
        fFirstWireProj[cs][TPCCount].resize(PlanesThisTPC);
        fOrthVectorsY[cs][TPCCount].resize(PlanesThisTPC);
        fOrthVectorsZ[cs][TPCCount].resize(PlanesThisTPC);
        fNPlanes[cs][TPCCount] = PlanesThisTPC;
        for (auto const& plane : Iterate<PlaneGeo>(TPC.ID())) {
          fPlaneIDs.insert(plane.ID());
          auto const PlaneCount = plane.ID().Plane;

          double ThisWirePitch = plane.WirePitch();
          fWireCounts[cs][TPCCount][PlaneCount] = plane.Nwires();

          WireGeo const& firstWire = plane.Wire(0);
          double const sth = firstWire.SinThetaZ(), cth = firstWire.CosThetaZ();

          auto WireCenter1 = firstWire.GetCenter();
          auto WireCenter2 = plane.Wire(1).GetCenter();

          // figure out if we need to flip the orthogonal vector
          // (should point from wire n -> n+1)
          double OrthY = cth, OrthZ = -sth;
          if (((WireCenter2.Y() - WireCenter1.Y()) * OrthY +
               (WireCenter2.Z() - WireCenter1.Z()) * OrthZ) < 0) {
            OrthZ *= -1;
            OrthY *= -1;
          }

          // Overall we are trying to build an expression that looks like
          //  int NearestWireNumber = round((worldPos.OrthVector - FirstWire.OrthVector)/WirePitch);
          // That runs as fast as humanly possible.
          // We predivide everything by the wire pitch so we don't do this in the loop.
          //
          // Putting this together into the useful constants we will use later per plane and tpc:
          fOrthVectorsY[cs][TPCCount][PlaneCount] = OrthY / ThisWirePitch;
          fOrthVectorsZ[cs][TPCCount][PlaneCount] = OrthZ / ThisWirePitch;

          fFirstWireProj[cs][TPCCount][PlaneCount] =
            WireCenter1.Y() * OrthY + WireCenter1.Z() * OrthZ;
          fFirstWireProj[cs][TPCCount][PlaneCount] /= ThisWirePitch;

          // now to count up wires in each plane and get first channel in each plane
          int WiresThisPlane = plane.Nwires();
          fWiresPerPlane[cs].at(TPCCount).push_back(WiresThisPlane);
          fPlaneBaselines[cs].at(TPCCount).push_back(RunningTotal);

          RunningTotal += WiresThisPlane;

          fFirstChannelInThisPlane[cs].at(TPCCount).push_back(fTopChannel);
          fTopChannel += WiresThisPlane;
          fFirstChannelInNextPlane[cs].at(TPCCount).push_back(fTopChannel);

        } // end loop over planes
      }   // end loop over TPCs
    }     // end loop over cryostats

    // calculate the total number of channels in the detector
    fNchannels = fTopChannel;

    MF_LOG_DEBUG("ChannelMapStandard") << "# of channels is " << fNchannels;
  }

  //----------------------------------------------------------------------------
  std::vector<WireID> WireReadoutStandardGeom::ChannelToWire(raw::ChannelID_t channel) const
  {
    std::vector<WireID> AllSegments;
    unsigned int cstat = 0;
    unsigned int tpc = 0;
    unsigned int plane = 0;
    unsigned int wire = 0;

    // first check if this channel ID is legal
    if (channel > fTopChannel)
      throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel << "\n";

    // then go find which plane, tpc and cryostat it is in from the information we stored earlier
    bool foundWid(false);
    for (unsigned int csloop = 0; csloop != fNcryostat; ++csloop) {
      for (unsigned int tpcloop = 0; tpcloop != fNTPC[csloop]; ++tpcloop) {
        for (unsigned int planeloop = 0;
             planeloop != fFirstChannelInNextPlane[csloop][tpcloop].size();
             ++planeloop) {
          if (channel < fFirstChannelInNextPlane[csloop][tpcloop][planeloop]) {
            cstat = csloop;
            tpc = tpcloop;
            plane = planeloop;
            wire = channel - fFirstChannelInThisPlane[cstat][tpcloop][planeloop];
            foundWid = true;
            break;
          }
          if (foundWid) break;
        } // end plane loop
        if (foundWid) break;
      } // end tpc loop
      if (foundWid) break;
    } // end cryostat loop

    WireID CodeWire(cstat, tpc, plane, wire);

    AllSegments.push_back(CodeWire);

    return AllSegments;
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t WireReadoutStandardGeom::Nchannels() const { return fNchannels; }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutStandardGeom::Nchannels(readout::ROPID const& ropid) const
  {
    if (!HasROP(ropid)) return 0;
    // The number of channels matches the number of wires. Life is easy.
    return WireCount(FirstWirePlaneInROP(ropid));
  }

  //----------------------------------------------------------------------------
  double WireReadoutStandardGeom::WireCoordinate(double YPos,
                                                 double ZPos,
                                                 PlaneID const& planeID) const
  {
    return YPos * AccessElement(fOrthVectorsY, planeID) +
           ZPos * AccessElement(fOrthVectorsZ, planeID) - AccessElement(fFirstWireProj, planeID);
  }

  //----------------------------------------------------------------------------
  WireID WireReadoutStandardGeom::NearestWireID(Point_t const& worldPos,
                                                PlaneID const& planeID) const
  {

    // This part is the actual calculation of the nearest wire number, where we assume
    //  uniform wire pitch and angle within a wireplane

    // add 0.5 to have the correct rounding
    int NearestWireNumber = int(0.5 + WireCoordinate(worldPos.Y(), worldPos.Z(), planeID));

    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if (NearestWireNumber < 0 || (unsigned int)NearestWireNumber >= WireCount(planeID)) {
      int wireNumber = NearestWireNumber; // save for the output

      if (NearestWireNumber < 0)
        NearestWireNumber = 0;
      else
        NearestWireNumber = WireCount(planeID) - 1;

      throw InvalidWireError("Geometry", planeID, wireNumber, NearestWireNumber)
        << "Can't Find Nearest Wire for position (" << worldPos.X() << "," << worldPos.Y() << ","
        << worldPos.Z() << ")"
        << " in plane " << std::string(planeID) << " approx wire number # " << wireNumber
        << " (capped from " << NearestWireNumber << ")\n";
    }

    return WireID(planeID, (WireID::WireID_t)NearestWireNumber);
  }

  //----------------------------------------------------------------------------
  // This method returns the channel number, assuming the numbering scheme
  // is heirachical - that is, channel numbers run in order, for example:
  //                                             (Ben J Oct 2011)
  //                    Wire1     | 0
  //           Plane1 { Wire2     | 1
  //    TPC1 {          Wire3     | 2
  //           Plane2 { Wire1     | 3   increasing channel number
  //                    Wire2     | 4     (with no gaps)
  //    TPC2 { Plane1 { Wire1     | 5
  //           Plane2 { Wire1     | 6
  //                    Wire2     v 7
  //
  raw::ChannelID_t WireReadoutStandardGeom::PlaneWireToChannel(WireID const& wireID) const
  {
    unsigned int const* pBaseLine = GetElementPtr(fPlaneBaselines, wireID);
    // This is the actual lookup part - first make sure coordinates are legal
    if (pBaseLine) {
      // if the channel has legal coordinates, its ID is given by the wire
      // number above the number of wires in lower planes, tpcs and cryostats
      return *pBaseLine + wireID.Wire;
    }
    else {
      // if the coordinates were bad, throw an exception
      throw cet::exception("WireReadoutStandardGeom")
        << "NO CHANNEL FOUND for " << std::string(wireID);
    }

    // made it here, that shouldn't happen, return raw::InvalidChannelID
    mf::LogWarning("WireReadoutStandardGeom")
      << "should not be at the point in the function, returning "
      << "invalid channel";
    return raw::InvalidChannelID;
  }

  //----------------------------------------------------------------------------
  SigType_t WireReadoutStandardGeom::SignalTypeForChannelImpl(raw::ChannelID_t const channel) const
  {

    // still assume one cryostat for now -- faster
    unsigned int nChanPerTPC = fNchannels / fNTPC[0];
    // casting wil trunc towards 0 -- faster than floor
    unsigned int tpc = channel / nChanPerTPC;
    //need number of planes to know Collection
    unsigned int PlanesThisTPC = fNPlanes[0][tpc];

    SigType_t sigt = kMysteryType;
    if ((channel >= fFirstChannelInThisPlane[0][tpc][0]) &&
        (channel < fFirstChannelInNextPlane[0][tpc][PlanesThisTPC - 2])) {
      sigt = kInduction;
    }
    else if ((channel >= fFirstChannelInThisPlane[0][tpc][PlanesThisTPC - 1]) &&
             (channel < fFirstChannelInNextPlane[0][tpc][PlanesThisTPC - 1])) {
      sigt = kCollection;
    }
    else
      mf::LogWarning("BadChannelSignalType")
        << "Channel " << channel << " not given signal type." << std::endl;

    return sigt;
  }

  //----------------------------------------------------------------------------
  std::set<PlaneID> const& WireReadoutStandardGeom::PlaneIDs() const { return fPlaneIDs; }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutStandardGeom::NTPCsets(readout::CryostatID const& cryoid) const
  {
    // return the same number as the number of TPCs
    return (cryoid.isValid && cryoid.Cryostat < fNTPC.size()) ? fNTPC[cryoid.Cryostat] : 0;
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutStandardGeom::MaxTPCsets() const { return MaxTPCs(); }

  //----------------------------------------------------------------------------
  bool WireReadoutStandardGeom::HasTPCset(readout::TPCsetID const& tpcsetid) const
  {
    return tpcsetid.TPCset < NTPCsets(tpcsetid);
  }

  //----------------------------------------------------------------------------
  readout::TPCsetID WireReadoutStandardGeom::TPCtoTPCset(TPCID const& tpcid) const
  {
    return ConvertTPCtoTPCset(tpcid);
  }

  //----------------------------------------------------------------------------
  std::vector<TPCID> WireReadoutStandardGeom::TPCsetToTPCs(readout::TPCsetID const& tpcsetid) const
  {
    std::vector<TPCID> IDs;
    if (tpcsetid.isValid) IDs.emplace_back(ConvertTPCsetToTPC(tpcsetid));
    return IDs;
  }

  //----------------------------------------------------------------------------
  TPCID WireReadoutStandardGeom::FirstTPCinTPCset(readout::TPCsetID const& tpcsetid) const
  {
    return ConvertTPCsetToTPC(tpcsetid);
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutStandardGeom::MaxTPCs() const
  {
    unsigned int max = 0;
    for (unsigned int nTPCs : fNTPC)
      if (nTPCs > max) max = nTPCs;
    return max;
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutStandardGeom::NROPs(readout::TPCsetID const& tpcsetid) const
  {
    if (!HasTPCset(tpcsetid)) return 0;
    return AccessElement(fNPlanes, FirstTPCinTPCset(tpcsetid));
  }

  //----------------------------------------------------------------------------
  unsigned int WireReadoutStandardGeom::MaxROPs() const
  {
    unsigned int max = 0;
    for (auto const& cryo_tpc : fNPlanes)
      for (unsigned int nPlanes : cryo_tpc)
        if (nPlanes > max) max = nPlanes;
    return max;
  }

  //----------------------------------------------------------------------------
  bool WireReadoutStandardGeom::HasROP(readout::ROPID const& ropid) const
  {
    return ropid.ROP < NROPs(ropid);
  }

  //----------------------------------------------------------------------------
  readout::ROPID WireReadoutStandardGeom::WirePlaneToROP(PlaneID const& planeid) const
  {
    return ConvertWirePlaneToROP(planeid);
  }

  //----------------------------------------------------------------------------
  std::vector<PlaneID> WireReadoutStandardGeom::ROPtoWirePlanes(readout::ROPID const& ropid) const
  {
    std::vector<PlaneID> IDs;
    if (ropid.isValid) IDs.emplace_back(FirstWirePlaneInROP(ropid));
    return IDs;
  }

  //----------------------------------------------------------------------------
  std::vector<TPCID> WireReadoutStandardGeom::ROPtoTPCs(readout::ROPID const& ropid) const
  {
    std::vector<TPCID> IDs;
    // we take the TPC set of the ROP and convert it straight into a TPC ID
    if (ropid.isValid) IDs.emplace_back(ConvertTPCsetToTPC(ropid.asTPCsetID()));
    return IDs;
  }

  //----------------------------------------------------------------------------
  readout::ROPID WireReadoutStandardGeom::ChannelToROP(raw::ChannelID_t channel) const
  {
    if (!raw::isValidChannelID(channel)) return {}; // invalid ROP returned

    // which wires does the channel cover?
    std::vector<WireID> wires = ChannelToWire(channel);

    // - none:
    if (wires.empty()) return {}; // invalid ROP returned

    // - one: maps its plane ID into a ROP ID
    return WirePlaneToROP(wires[0]);
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t WireReadoutStandardGeom::FirstChannelInROP(readout::ROPID const& ropid) const
  {
    if (!ropid.isValid) return raw::InvalidChannelID;
    return (raw::ChannelID_t)AccessElement(fPlaneBaselines, ConvertROPtoWirePlane(ropid));
  }

  //----------------------------------------------------------------------------
  PlaneID WireReadoutStandardGeom::FirstWirePlaneInROP(readout::ROPID const& ropid) const
  {
    return ConvertROPtoWirePlane(ropid);
  }

} // namespace
