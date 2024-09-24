/**
 * \file GeoDirectedLine.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class DirectedLine
 *
 * @author David Caratelli
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEODIRECTEDLINE_H
#define BASICTOOL_GEODIRECTEDLINE_H

#include "larcorealg/GeoAlgo/GeoLine.h"
#include "larcorealg/GeoAlgo/GeoVector.h"

namespace geoalgo {

  class HalfLine;

  /**
     \class DirectedLine
     @brief Representation of a 3D infinite line.
     Defines an infinite 3D line with a point and a direction.
     Line points are constructed like this:
     (pt, dir) -> (pt, pt+dir)
     It hides the point attributes from users for protecting the dimensionality.
  */
  class DirectedLine : public Line {

  public:
    /// Default ctor
    DirectedLine();

    /// Alternative ctor (1)
    DirectedLine(double const x,
                 double const y,
                 double const z,
                 double const dirx,
                 double const diry,
                 double const dirz);

    /// Altenartive ctor (2)
    DirectedLine(Point_t const& pt, Vector_t const& dir);

    /// Alternative ctor (3)
    DirectedLine(HalfLine const& l);

    /// Alternative ctor using template (3)
    template <class T, class U>
    DirectedLine(T const& pt, U const& dir) : Line(Point_t(pt), Point_t(pt + dir))
    {}

    Vector_t Dir() const;
  };

  typedef DirectedLine DirectedLine_t;
}
#endif
/** @} */ // end of doxygen group
