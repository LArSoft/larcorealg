#include "larcorealg/GeoAlgo/GeoLineSegment.h"
#include "larcorealg/GeoAlgo/GeoAlgoException.h"

namespace geoalgo {

  LineSegment::LineSegment() : _start(3), _end(3), _dir(3)
  {
    DirReset();
  }

  LineSegment::LineSegment(double const start_x,
                           double const start_y,
                           double const start_z,
                           double const end_x,
                           double const end_y,
                           double const end_z)
    : _start(start_x, start_y, start_z), _end(end_x, end_y, end_z), _dir(3)
  {
    DirReset();
  }

  LineSegment::LineSegment(Point_t const& start, Point_t const& end)
    : _start(start), _end(end), _dir(3)
  {
    if (start.size() != 3 || end.size() != 3)
      throw GeoAlgoException("LineSegment ctor accepts only 3D Point!");
    DirReset();
  }

  Point_t const& LineSegment::Start() const
  {
    return _start;
  }

  Point_t const& LineSegment::End() const
  {
    return _end;
  }

  Vector_t const LineSegment::Dir() const
  {
    return _dir;
  }

  void LineSegment::Start(double const x, double const y, double const z)
  {
    _start[0] = x;
    _start[1] = y;
    _start[2] = z;
    DirReset();
  }

  void LineSegment::End(double const x, double const y, double const z)
  {
    _end[0] = x;
    _end[1] = y;
    _end[2] = z;
    DirReset();
  }

  void LineSegment::DirReset()
  {
    _dir = _end - _start;
  }

}
