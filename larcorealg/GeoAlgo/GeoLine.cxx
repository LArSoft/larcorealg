#include "larcorealg/GeoAlgo/GeoLine.h"
#include "larcorealg/GeoAlgo/GeoAlgoException.h"

namespace geoalgo {

  Line::Line() : _pt1(3), _pt2(3) {}

  Line::Line(double const x1,
             double const y1,
             double const z1,
             double const x2,
             double const y2,
             double const z2)
    : _pt1(x1, y1, z1), _pt2(x2, y2, z2)
  {
    check_and_raise(_pt1, _pt2);
  }

  Line::Line(Point_t const& pt1, Point_t const& pt2) : _pt1(pt1), _pt2(pt2)
  {
    check_and_raise(pt1, pt2);
  }

  Point_t const& Line::Pt1() const { return _pt1; }
  Point_t const& Line::Pt2() const { return _pt2; }

  void Line::Pt1(double const x, double const y, double const z)
  {
    _pt1[0] = x;
    _pt1[1] = y;
    _pt1[2] = z;
    check_and_raise(_pt1, _pt2);
  }

  void Line::Pt2(double const x, double const y, double const z)
  {
    _pt2[0] = x;
    _pt2[1] = y;
    _pt2[2] = z;
    check_and_raise(_pt1, _pt2);
  }

  void Line::check_and_raise(Point_t const& p1, Point_t const& p2) const
  {
    if (p1.size() != 3)
      throw GeoAlgoException("<<check_and_raise>> Pt1 is not 3 dimensional point!");
    if (p2.size() != 3)
      throw GeoAlgoException("<<check_and_raise>> Pt2 is not 3 dimensional point!");
    if (p1 == p2)
      throw GeoAlgoException("<<check_and_raise>> Two identical points not allowed for Line ctor!");
  }

}
