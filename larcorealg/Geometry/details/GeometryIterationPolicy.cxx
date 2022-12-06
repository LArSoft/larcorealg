#include "larcorealg/Geometry/details/GeometryIterationPolicy.h"
#include "larcorealg/Geometry/GeometryCore.h"

namespace geo::details {

  GeometryIterationPolicy::GeometryIterationPolicy(GeometryCore const* geom) : fGeom{geom} {}

  unsigned int GeometryIterationPolicy::NSiblings(CryostatID const& id) const
  {
    return fGeom->NSiblingElements(id);
  }

  unsigned int GeometryIterationPolicy::NSiblings(TPCID const& id) const
  {
    return fGeom->NSiblingElements(id);
  }

  unsigned int GeometryIterationPolicy::NSiblings(PlaneID const& id) const
  {
    return fGeom->NSiblingElements(id);
  }

  unsigned int GeometryIterationPolicy::NSiblings(WireID const& id) const
  {
    return fGeom->NSiblingElements(id);
  }

  // CryostatID
  CryostatID GeometryIterationPolicy::EndCryostatID() const
  {
    return CryostatID{fGeom->Ncryostats()};
  }

  // TPCID
  TPCID GeometryIterationPolicy::EndTPCID() const
  {
    auto id = TPCID::first();
    if (fGeom->MaxTPCs() != 0) { id.Cryostat = fGeom->Ncryostats(); }
    return id;
  }

  TPCID GeometryIterationPolicy::EndTPCID(CryostatID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

  // PlaneID
  PlaneID GeometryIterationPolicy::EndPlaneID() const
  {
    auto id = PlaneID::first();
    if (fGeom->MaxPlanes() != 0) { id.Cryostat = fGeom->Ncryostats(); }
    return id;
  }

  PlaneID GeometryIterationPolicy::EndPlaneID(CryostatID const& id) const
  {
    return {EndTPCID(id), 0};
  }

  PlaneID GeometryIterationPolicy::EndPlaneID(TPCID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

  // WireID
  WireID GeometryIterationPolicy::EndWireID() const
  {
    auto id = WireID::first();
    if (fGeom->MaxWires() != 0) { id.Cryostat = fGeom->Ncryostats(); }
    return id;
  }

  WireID GeometryIterationPolicy::EndWireID(CryostatID const& id) const
  {
    return {EndPlaneID(id), 0};
  }

  WireID GeometryIterationPolicy::EndWireID(TPCID const& id) const { return {EndPlaneID(id), 0}; }

  WireID GeometryIterationPolicy::EndWireID(PlaneID const& id) const
  {
    return {GetNextID(id, *this), 0};
  }

  // ......
  CryostatGeo const* getElementPtr(GeometryCore const* geom, CryostatID const& id)
  {
    assert(geom);
    return geom->GetElementPtr(id);
  }

  TPCGeo const* getElementPtr(GeometryCore const* geom, TPCID const& id)
  {
    assert(geom);
    return geom->GetElementPtr(id);
  }

  PlaneGeo const* getElementPtr(GeometryCore const* geom, PlaneID const& id)
  {
    assert(geom);
    return geom->GetElementPtr(id);
  }

  WireGeo const* getElementPtr(GeometryCore const* geom, WireID const& id)
  {
    assert(geom);
    return geom->GetElementPtr(id);
  }

  bool validElement(GeometryCore const* geom, CryostatID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

  bool validElement(GeometryCore const* geom, TPCID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

  bool validElement(GeometryCore const* geom, PlaneID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

  bool validElement(GeometryCore const* geom, WireID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

}
