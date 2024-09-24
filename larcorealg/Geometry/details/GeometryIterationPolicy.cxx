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

  bool validElement(GeometryCore const* geom, CryostatID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

  bool validElement(GeometryCore const* geom, TPCID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

}
