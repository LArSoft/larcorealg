////////////////////////////////////////////////////////////////////////
/// \file  WireReadoutGeom.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "cetlib_except/exception.h"

#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/AuxDetWireReadoutGeom.h"

namespace geo {

  //----------------------------------------------------------------------------
  size_t AuxDetWireReadoutGeom::NearestAuxDet(Point_t const& point,
                                              std::vector<AuxDetGeo> const& auxDets,
                                              double tolerance) const
  {
    for (size_t a = 0; a < auxDets.size(); ++a) {
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

    // throw an exception because we couldn't find the sensitive volume
    throw cet::exception("AuxDetWireReadoutGeom") << "Can't find AuxDet for position (" << point.X()
                                                  << "," << point.Y() << "," << point.Z() << ")\n";
  }

  //----------------------------------------------------------------------------
  size_t AuxDetWireReadoutGeom::NearestSensitiveAuxDet(Point_t const& point,
                                                       std::vector<AuxDetGeo> const& auxDets,
                                                       size_t& ad,
                                                       double tolerance) const
  {
    ad = NearestAuxDet(point, auxDets, tolerance);

    AuxDetGeo const& adg = auxDets[ad];

    for (size_t a = 0; a < adg.NSensitiveVolume(); ++a) {
      AuxDetSensitiveGeo const& adsg = adg.SensitiveVolume(a);
      auto const localPoint = adsg.toLocalCoords(point);

      double const HalfCenterWidth = 0.5 * (adsg.HalfWidth1() + adsg.HalfWidth2());

      if (localPoint.Z() >= -(adsg.Length() / 2 + tolerance) &&
          localPoint.Z() <= adsg.Length() / 2 + tolerance &&
          localPoint.Y() >= -(adsg.HalfHeight() + tolerance) &&
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
  size_t AuxDetWireReadoutGeom::ChannelToAuxDet(std::vector<AuxDetGeo> const& /* auxDets */,
                                                std::string const& detName,
                                                uint32_t /*channel*/) const
  {
    // loop over the map of AuxDet names to Geo object numbers to determine which auxdet
    // we have.  If no name in the map matches the provided string, throw an exception;
    // the list of AuxDetGeo passed as argument is ignored!
    // Note that fADGeoToName must have been updated by a derived class.
    for (auto itr : fADGeoToName)
      if (itr.second.compare(detName) == 0) return itr.first;

    throw cet::exception("Geometry") << "No AuxDetGeo matching name: " << detName;
  }

  //----------------------------------------------------------------------------
  // the first member of the pair is the index in the auxDets vector for the AuxDetGeo,
  // the second member is the index in the vector of AuxDetSensitiveGeos for that AuxDetGeo
  std::pair<size_t, size_t> AuxDetWireReadoutGeom::ChannelToSensitiveAuxDet(
    std::vector<AuxDetGeo> const& auxDets,
    std::string const& detName,
    uint32_t channel) const
  {
    size_t adGeoIdx = ChannelToAuxDet(auxDets, detName, channel);

    // look for the index of the sensitive volume for the given channel
    if (fADGeoToChannelAndSV.count(adGeoIdx) > 0) {

      auto itr = fADGeoToChannelAndSV.find(adGeoIdx);

      // get the vector of channels to AuxDetSensitiveGeo index
      if (channel < itr->second.size())
        return std::make_pair(adGeoIdx, itr->second[channel].second);

      throw cet::exception("Geometry")
        << "Given AuxDetSensitive channel, " << channel
        << ", cannot be found in vector associated to AuxDetGeo index: " << adGeoIdx
        << ". Vector has size " << itr->second.size();
    }

    throw cet::exception("Geometry") << "Given AuxDetGeo with index " << adGeoIdx
                                     << " does not correspond to any vector of sensitive volumes";
  }

}
