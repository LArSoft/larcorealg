////////////////////////////////////////////////////////////////////////
/// @file  WireReadoutGeom.cxx
/// @brief Interface to algorithm class for a specific detector channel mapping
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/AuxDetReadoutGeom.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <limits>

namespace {
  constexpr auto invalid_index = std::numeric_limits<std::size_t>::max();
}

namespace geo {
  AuxDetReadoutInitializers AuxDetInitializer::init(std::vector<AuxDetGeo> const& ads) const
  {
    return initialize(ads);
  }

  //----------------------------------------------------------------------------
  AuxDetReadoutGeom::AuxDetReadoutGeom(AuxDetReadoutInitializers initializers)
    : fADGeoToName{std::move(initializers.ADGeoToName)}
    , fNameToADGeo{std::move(initializers.NameToADGeo)}
    , fADGeoToChannelAndSV{std::move(initializers.ADGeoToChannelAndSV)}
  {}

  //----------------------------------------------------------------------------
  std::size_t AuxDetReadoutGeom::NearestAuxDet(Point_t const& point,
                                               std::vector<AuxDetGeo> const& auxDets,
                                               double tolerance,
                                               bool throwIfAbsent) const
  {
    for (std::size_t a = 0; a < auxDets.size(); ++a) {
      auto const localPoint = auxDets[a].toLocalCoords(point);
      double const HalfCenterWidth = 0.5 * (auxDets[a].HalfWidth1() + auxDets[a].HalfWidth2());

      if (localPoint.Z() >= (-auxDets[a].Length() / 2 - tolerance) &&
          localPoint.Z() <= (auxDets[a].Length() / 2 + tolerance) &&
          localPoint.Y() >= (-auxDets[a].HalfHeight() - tolerance) &&
          localPoint.Y() <= (auxDets[a].HalfHeight() + tolerance) &&
          // if AuxDet a is a box, then HalfSmallWidth = HalfWidth
          localPoint.X() >= (-HalfCenterWidth +
                             localPoint.Z() * (HalfCenterWidth - auxDets[a].HalfWidth2()) /
                               (0.5 * auxDets[a].Length()) -
                             tolerance) &&
          localPoint.X() <= (HalfCenterWidth -
                             localPoint.Z() * (HalfCenterWidth - auxDets[a].HalfWidth2()) /
                               (0.5 * auxDets[a].Length()) +
                             tolerance))
        return a;
    } // for loop over AudDet a

    std::ostringstream msg;
    msg << "Can't find AuxDet for position " << point;

    if (throwIfAbsent) { throw cet::exception("AuxDetReadoutGeom") << msg.str(); }

    mf::LogDebug("AuxDetReadoutGeom") << msg.str();
    return invalid_index;
  }

  // ----------------------------------------------------------------------------
  std::size_t AuxDetReadoutGeom::NearestSensitiveAuxDet(Point_t const& point,
                                                        std::vector<AuxDetGeo> const& auxDets,
                                                        double tolerance,
                                                        bool throwIfAbsent) const
  {
    std::size_t const auxDetIdx = NearestAuxDet(point, auxDets, tolerance, throwIfAbsent);
    if (auxDetIdx == invalid_index) { return invalid_index; }

    AuxDetGeo const& adg = auxDets[auxDetIdx];
    for (std::size_t a = 0; a < adg.NSensitiveVolume(); ++a) {
      AuxDetSensitiveGeo const& adsg = adg.SensitiveVolume(a);
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

    std::ostringstream msg;
    msg << "Can't find AuxDetSensitive for position " << point;

    if (throwIfAbsent) { throw cet::exception("AuxDetReadoutGeom") << msg.str(); }

    mf::LogDebug("AuxDetReadoutGeom") << msg.str();
    return invalid_index;
  }

  //----------------------------------------------------------------------------
  // the first member of the pair is the index in the auxDets vector for the AuxDetGeo,
  // the second member is the index in the vector of AuxDetSensitiveGeos for that AuxDetGeo
  std::pair<std::size_t, std::size_t> AuxDetReadoutGeom::ChannelToSensitiveAuxDet(
    std::string const& detName,
    std::uint32_t channel) const
  {
    std::size_t adGeoIdx = DetNameToAuxDet(detName);

    // look for the index of the sensitive volume for the given channel
    auto it = fADGeoToChannelAndSV.find(adGeoIdx);
    if (it == fADGeoToChannelAndSV.cend()) {
      throw cet::exception("Geometry")
        << "Given AuxDetSensitive channel, " << channel
        << ", cannot be found in vector associated to AuxDetGeo index: " << adGeoIdx
        << ". Vector has size " << it->second.size();
    }

    // get the vector of channels to AuxDetSensitiveGeo index
    if (channel < it->second.size()) return {adGeoIdx, it->second[channel].second};

    throw cet::exception("Geometry") << "Given AuxDetGeo with index " << adGeoIdx
                                     << " does not correspond to any vector of sensitive volumes";
  }

  //----------------------------------------------------------------------------
  geo::Point_t AuxDetReadoutGeom::AuxDetChannelToPosition(
    std::uint32_t const channel,
    std::string const& auxDetName,
    std::vector<geo::AuxDetGeo> const& auxDets) const
  {
    // Figure out which detector we are in
    auto ad_it = fNameToADGeo.find(auxDetName);
    if (ad_it == fNameToADGeo.cend()) {
      throw cet::exception("CRTWireReadoutGeom") << "No AuxDetGeo with name " << auxDetName;
    }
    std::size_t const ad = ad_it->second;

    // Get the vector of channel and sensitive volume pairs
    auto it = fADGeoToChannelAndSV.find(ad);
    if (it == fADGeoToChannelAndSV.end()) {
      throw cet::exception("CRTWireReadoutGeom") << "No entry in channel and sensitive volume"
                                                 << " map for AuxDet index " << ad << " bail";
    }

    // Loop over the vector of channel and sensitive volumes to determine the sensitive
    // volume for this channel. Then get the origin of the sensitive volume in the world
    // coordinate system.
    for (auto const& [ch, svindex] : it->second) {
      if (ch == channel) { return auxDets[ad].SensitiveVolume(svindex).GetCenter(); }
    }

    return {};
  }

  //----------------------------------------------------------------------------
  std::size_t AuxDetReadoutGeom::DetNameToAuxDet(std::string const& detName) const
  {
    // loop over the map of AuxDet names to Geo object numbers to determine which auxdet
    // we have.  If no name in the map matches the provided string, throw an exception;
    // the list of AuxDetGeo passed as argument is ignored!  Note that fADGeoToName must
    // have been updated by a derived class.
    auto it = std::find_if(fADGeoToName.begin(), fADGeoToName.end(), [&detName](auto const& pr) {
      return pr.second == detName;
    });
    if (it == fADGeoToName.end()) {
      throw cet::exception("Geometry") << "No AuxDetGeo matching name: " << detName;
    }
    return it->first;
  }
}
