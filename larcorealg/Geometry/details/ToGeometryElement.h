#ifndef LARCOREALG_GEOMETRY_DETAILS_TOGEOMETRYELEMENT_H
#define LARCOREALG_GEOMETRY_DETAILS_TOGEOMETRYELEMENT_H

// LArSoft libraries
#include "larcorealg/Geometry/details/element_iterators.h"
#include "larcorealg/Geometry/fwd.h"

namespace geo::details {
  class ToGeometryElement {
  public:
    explicit ToGeometryElement(GeometryCore const* geom) : fGeom{geom} {}
    template <typename T, typename Iterator>
    auto transform(Iterator const& iterator) const
    {
      return geometry_element_iterator<GeometryCore, T, Iterator>{fGeom, iterator};
    }

  private:
    GeometryCore const* fGeom;
  };

  class ToReadoutGeometryElement {
  public:
    explicit ToReadoutGeometryElement(WireReadoutGeom const* wireReadoutGeom)
      : fWireReadoutGeom{wireReadoutGeom}
    {}

    template <typename T, typename Iterator>
    auto transform(Iterator const& iterator) const
    {
      return geometry_element_iterator<WireReadoutGeom, T, Iterator>{fWireReadoutGeom, iterator};
    }

  private:
    WireReadoutGeom const* fWireReadoutGeom;
  };
}

#endif // LARCOREALG_GEOMETRY_DETAILS_TOGEOMETRYELEMENT_H
