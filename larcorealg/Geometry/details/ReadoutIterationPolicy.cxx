#include "larcorealg/Geometry/details/ReadoutIterationPolicy.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"

namespace geo::details {

  ReadoutIterationPolicy::ReadoutIterationPolicy(GeometryCore const* geom,
                                                 ChannelMapAlg const* channelMapAlg)
    : fGeom{geom}, fChannelMapAlg{channelMapAlg}
  {}

  unsigned int ReadoutIterationPolicy::NSiblings(CryostatID const& id) const
  {
    return fGeom->NSiblingElements(id);
  }

  unsigned int ReadoutIterationPolicy::NSiblings(readout::TPCsetID const& id) const
  {
    return fChannelMapAlg->NTPCsets(id);
  }

  unsigned int ReadoutIterationPolicy::NSiblings(readout::ROPID const& id) const
  {
    return fChannelMapAlg->NROPs(id);
  }

  // CryostatID
  CryostatID ReadoutIterationPolicy::EndCryostatID() const
  {
    return CryostatID{fGeom->Ncryostats()};
  }

  // TPCsetID
  readout::TPCsetID ReadoutIterationPolicy::EndTPCsetID() const { return {EndCryostatID(), 0}; }

  readout::TPCsetID ReadoutIterationPolicy::EndTPCsetID(CryostatID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

  // ROPID
  readout::ROPID ReadoutIterationPolicy::EndROPID() const { return {EndTPCsetID(), 0}; }

  readout::ROPID ReadoutIterationPolicy::EndROPID(CryostatID const& id) const
  {
    return {EndTPCsetID(id), 0};
  }

  readout::ROPID ReadoutIterationPolicy::EndROPID(readout::TPCsetID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

}
