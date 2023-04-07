#include "larcorealg/GeoAlgo/GeoDirectedLine.h"
#include "larcorealg/GeoAlgo/GeoHalfLine.h"

namespace geoalgo {

  DirectedLine::DirectedLine() : Line() {}

  DirectedLine::DirectedLine(double const x,
                             double const y,
                             double const z,
                             double const dirx,
                             double const diry,
                             double const dirz)
    : Line(x, y, z, x + dirx, y + diry, z + dirz)
  {
    check_and_raise(_pt1, _pt2);
  }

  DirectedLine::DirectedLine(Point_t const& pt, Vector_t const& dir) : Line(pt, pt + dir)
  {
    check_and_raise(_pt1, _pt2);
  }

  DirectedLine::DirectedLine(HalfLine const& l) : Line(l.Start(), l.Start() + l.Dir())
  {
    check_and_raise(_pt1, _pt2);
  }

  Vector_t DirectedLine::Dir() const { return _pt2 - _pt1; }

}
