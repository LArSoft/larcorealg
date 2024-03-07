#include "larcorealg/Geometry/details/ReadoutIterationPolicy.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

namespace geo::details {

  ReadoutIterationPolicy::ReadoutIterationPolicy(GeometryCore const* geom,
                                                 WireReadoutGeom const* wireReadoutGeom)
    : GeometryIterationPolicy{geom}, fGeom{geom}, fWireReadoutGeom{wireReadoutGeom}
  {}

  unsigned int ReadoutIterationPolicy::NSiblings(PlaneID const& id) const
  {
    return fWireReadoutGeom->NSiblingElements(id);
  }

  unsigned int ReadoutIterationPolicy::NSiblings(WireID const& id) const
  {
    return fWireReadoutGeom->NSiblingElements(id);
  }

  unsigned int ReadoutIterationPolicy::NSiblings(readout::TPCsetID const& id) const
  {
    return fWireReadoutGeom->NTPCsets(id);
  }

  unsigned int ReadoutIterationPolicy::NSiblings(readout::ROPID const& id) const
  {
    return fWireReadoutGeom->NROPs(id);
  }

  // PlaneID
  PlaneID ReadoutIterationPolicy::EndPlaneID() const
  {
    auto id = PlaneID::first();
    if (fWireReadoutGeom->MaxPlanes() != 0) { id.Cryostat = fGeom->Ncryostats(); }
    return id;
  }

  PlaneID ReadoutIterationPolicy::EndPlaneID(CryostatID const& id) const
  {
    return {EndTPCID(id), 0};
  }

  PlaneID ReadoutIterationPolicy::EndPlaneID(TPCID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

  // WireID
  WireID ReadoutIterationPolicy::EndWireID() const
  {
    auto id = WireID::first();
    if (fWireReadoutGeom->MaxWires() != 0) { id.Cryostat = fGeom->Ncryostats(); }
    return id;
  }

  WireID ReadoutIterationPolicy::EndWireID(CryostatID const& id) const
  {
    return {EndPlaneID(id), 0};
  }

  WireID ReadoutIterationPolicy::EndWireID(TPCID const& id) const
  {
    return {EndPlaneID(id), 0};
  }

  WireID ReadoutIterationPolicy::EndWireID(PlaneID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

  // TPCsetID
  readout::TPCsetID ReadoutIterationPolicy::EndTPCsetID() const
  {
    return {EndCryostatID(), 0};
  }

  readout::TPCsetID ReadoutIterationPolicy::EndTPCsetID(CryostatID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

  // ROPID
  readout::ROPID ReadoutIterationPolicy::EndROPID() const
  {
    return {EndTPCsetID(), 0};
  }

  readout::ROPID ReadoutIterationPolicy::EndROPID(CryostatID const& id) const
  {
    return {EndTPCsetID(id), 0};
  }

  readout::ROPID ReadoutIterationPolicy::EndROPID(readout::TPCsetID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

  // ......
  PlaneGeo const* getElementPtr(WireReadoutGeom const* wireReadoutGeom, PlaneID const& id)
  {
    assert(wireReadoutGeom);
    return wireReadoutGeom->PlanePtr(id);
  }

  WireGeo const* getElementPtr(WireReadoutGeom const* wireReadoutGeom, WireID const& id)
  {
    assert(wireReadoutGeom);
    return wireReadoutGeom->WirePtr(id);
  }

  bool validElement(WireReadoutGeom const* wireReadoutGeom, PlaneID const& id)
  {
    return wireReadoutGeom && getElementPtr(wireReadoutGeom, id) != nullptr;
  }

  bool validElement(WireReadoutGeom const* wireReadoutGeom, WireID const& id)
  {
    return wireReadoutGeom && getElementPtr(wireReadoutGeom, id) != nullptr;
  }

}
