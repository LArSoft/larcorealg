#ifndef LARCOREALG_GEOMETRY_DETAILS_GEOMETRYITERATIONPOLICY_H
#define LARCOREALG_GEOMETRY_DETAILS_GEOMETRYITERATIONPOLICY_H

// LArSoft libraries
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace geo::details {
  class GeometryIterationPolicy {
  public:
    GeometryIterationPolicy() = default; // Required for subranges
    explicit GeometryIterationPolicy(GeometryCore const* geom);
    unsigned int NSiblings(CryostatID const& id) const;
    unsigned int NSiblings(TPCID const& id) const;
    unsigned int NSiblings(PlaneID const& id) const;
    unsigned int NSiblings(WireID const& id) const;

    template <typename GeoID>
    GeoID GetEndID() const;

    template <typename GeoID, typename ContextID>
    GeoID GetEndID(ContextID const& id) const;

  private:
    CryostatID EndCryostatID() const;

    TPCID EndTPCID() const;
    TPCID EndTPCID(CryostatID const& id) const;

    PlaneID EndPlaneID() const;
    PlaneID EndPlaneID(CryostatID const& id) const;
    PlaneID EndPlaneID(TPCID const& id) const;

    WireID EndWireID() const;
    WireID EndWireID(CryostatID const& id) const;
    WireID EndWireID(TPCID const& id) const;
    WireID EndWireID(PlaneID const& id) const;

    GeometryCore const* fGeom{nullptr};
  };

  /// Initializes the specified ID with the invalid ID after the last cryostat
  template <>
  inline CryostatID GeometryIterationPolicy::GetEndID<CryostatID>() const
  {
    return EndCryostatID();
  }

  // TPCID
  template <>
  inline TPCID GeometryIterationPolicy::GetEndID<TPCID>() const
  {
    return EndTPCID();
  }

  template <>
  inline TPCID GeometryIterationPolicy::GetEndID<TPCID, CryostatID>(CryostatID const& id) const
  {
    return EndTPCID(id);
  }

  // PlaneID
  template <>
  inline PlaneID GeometryIterationPolicy::GetEndID<PlaneID>() const
  {
    return EndPlaneID();
  }

  template <>
  inline PlaneID GeometryIterationPolicy::GetEndID<PlaneID, CryostatID>(CryostatID const& id) const
  {
    return EndPlaneID(id);
  }

  template <>
  inline PlaneID GeometryIterationPolicy::GetEndID<PlaneID, TPCID>(TPCID const& id) const
  {
    return EndPlaneID(id);
  }

  // WireID
  template <>
  inline WireID GeometryIterationPolicy::GetEndID<WireID>() const
  {
    return EndWireID();
  }

  template <>
  inline WireID GeometryIterationPolicy::GetEndID<WireID, CryostatID>(CryostatID const& id) const
  {
    return EndWireID(id);
  }

  template <>
  inline WireID GeometryIterationPolicy::GetEndID<WireID, TPCID>(TPCID const& id) const
  {
    return EndWireID(id);
  }

  template <>
  inline WireID GeometryIterationPolicy::GetEndID<WireID, PlaneID>(PlaneID const& id) const
  {
    return EndWireID(id);
  }

  CryostatGeo const* getElementPtr(GeometryCore const* geom, CryostatID const& id);
  TPCGeo const* getElementPtr(GeometryCore const* geom, TPCID const& id);
  PlaneGeo const* getElementPtr(GeometryCore const* geom, PlaneID const& id);
  WireGeo const* getElementPtr(GeometryCore const* geom, WireID const& id);

  bool validElement(GeometryCore const* geom, CryostatID const& id);
  bool validElement(GeometryCore const* geom, TPCID const& id);
  bool validElement(GeometryCore const* geom, PlaneID const& id);
  bool validElement(GeometryCore const* geom, WireID const& id);
}

#endif // LARCOREALG_GEOMETRY_DETAILS_GEOMETRYITERATIONPOLICY_H
