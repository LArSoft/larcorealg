#ifndef LARCOREALG_GEOMETRY_DETAILS_READOUTITERATIONPOLICY_H
#define LARCOREALG_GEOMETRY_DETAILS_READOUTITERATIONPOLICY_H

// LArSoft libraries
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

namespace geo::details {
  class ReadoutIterationPolicy {
  public:
    ReadoutIterationPolicy() = default; // Required for subranges
    ReadoutIterationPolicy(GeometryCore const* geom, ChannelMapAlg const* channelMapAlg);
    unsigned int NSiblings(CryostatID const& id) const;
    unsigned int NSiblings(readout::TPCsetID const& id) const;
    unsigned int NSiblings(readout::ROPID const& id) const;

    template <typename GeoID>
    GeoID GetEndID() const;

    template <typename GeoID, typename ContextID>
    GeoID GetEndID(ContextID const& id) const;

  private:
    CryostatID EndCryostatID() const;

    readout::TPCsetID EndTPCsetID() const;
    readout::TPCsetID EndTPCsetID(CryostatID const& id) const;

    readout::ROPID EndROPID() const;
    readout::ROPID EndROPID(CryostatID const& id) const;
    readout::ROPID EndROPID(readout::TPCsetID const& id) const;

    GeometryCore const* fGeom{nullptr};
    ChannelMapAlg const* fChannelMapAlg{nullptr};
  };

  // CryostatID
  template <>
  inline CryostatID ReadoutIterationPolicy::GetEndID<CryostatID>() const
  {
    return EndCryostatID();
  }

  // TPCsetID
  template <>
  inline readout::TPCsetID ReadoutIterationPolicy::GetEndID<readout::TPCsetID>() const
  {
    return EndTPCsetID();
  }

  template <>
  inline readout::TPCsetID ReadoutIterationPolicy::GetEndID<readout::TPCsetID, CryostatID>(
    CryostatID const& id) const
  {
    return EndTPCsetID(id);
  }

  // ROPID
  template <>
  inline readout::ROPID ReadoutIterationPolicy::GetEndID<readout::ROPID>() const
  {
    return EndROPID();
  }

  template <>
  inline readout::ROPID ReadoutIterationPolicy::GetEndID<readout::ROPID, CryostatID>(
    CryostatID const& id) const
  {
    return EndROPID(id);
  }

  template <>
  inline readout::ROPID ReadoutIterationPolicy::GetEndID<readout::ROPID, readout::TPCsetID>(
    readout::TPCsetID const& id) const
  {
    return EndROPID(id);
  }

}

#endif // LARCOREALG_GEOMETRY_DETAILS_READOUTITERATIONPOLICY_H
