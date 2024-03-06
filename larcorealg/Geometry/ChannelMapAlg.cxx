////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapAlg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/AuxDetGeo.h"

namespace geo {

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
    // By defualt, 1 channel per optical detector
    return 1;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapAlg::OpChannel(unsigned int detNum, unsigned int /* channel */) const
  {
    return detNum;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapAlg::OpDetFromOpChannel(unsigned int opChannel) const
  {
    return opChannel;
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
    if (opChannel >= this->NOpChannels(NOpDets)) return false;

    // Check opdet number
    unsigned int opdet = this->OpDetFromOpChannel(opChannel);
    if (opdet >= NOpDets) return false;

    // Check hardware channel number
    unsigned int hChan = this->HardwareChannelFromOpChannel(opChannel);
    if (hChan >= this->NOpHardwareChannels(opdet)) return false;

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
    size_t adGeoIdx = this->ChannelToAuxDet(auxDets, detName, channel);

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

  geo::SigType_t ChannelMapAlg::SignalTypeForChannel(raw::ChannelID_t const channel) const
  {
    return SignalTypeForChannelImpl(channel);
  }

  geo::SigType_t ChannelMapAlg::SignalTypeForROPID(readout::ROPID const& ropid) const
  {
    return SignalTypeForROPIDImpl(ropid);
  }

  geo::SigType_t ChannelMapAlg::SignalTypeForROPIDImpl(readout::ROPID const& ropid) const
  {
    return SignalTypeForChannel(FirstChannelInROP(ropid));
  }
}
