#include "larcorealg/GeoAlgo/GeoCone.h"
#include "larcorealg/GeoAlgo/GeoAlgoException.h" // for GeoAlgoException

#include <math.h>
#include <sstream>

namespace geoalgo {

  Cone::Cone() : HalfLine()
  {
    _length = 1;
    _radius = 1;
    _angle = atan(_radius / _length);
  }

  Cone::Cone(double const x,
             double const y,
             double const z,
             double const dirx,
             double const diry,
             double const dirz,
             double const length,
             double const radius)
    : HalfLine(x, y, z, dirx, diry, dirz)
  {
    if (length == 0) {
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
          << " Cone Length cannot be 0." << std::endl;
      throw GeoAlgoException(msg.str());
    }
    _length = length;
    _radius = radius;
    _angle = atan(_radius / _length);
  }

  Cone::Cone(Point_t const& start, Vector_t const& dir, double const length, double const radius)
    : HalfLine(start, dir)
  {
    if (length == 0) {
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
          << " Cone Length cannot be 0." << std::endl;
      throw GeoAlgoException(msg.str());
    }
    _length = length;
    _radius = radius;
    _angle = atan(_radius / _length);
  }

  double Cone::Length() const { return _length; }

  double Cone::Radius() const { return _radius; }

  double Cone::Angle() const { return _angle; }

  void Cone::Length(double const l)
  {
    if (l == 0) {
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
          << " Cone Length cannot be 0." << std::endl;
      throw GeoAlgoException(msg.str());
    }
    _length = l;
    _angle = atan(_radius / _length);
  }

  void Cone::Radius(double const r)
  {
    _radius = r;
    _angle = atan(_radius / _length);
  }
}
