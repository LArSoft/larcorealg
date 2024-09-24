#include "larcorealg/GeoAlgo/GeoCylinder.h"
#include "larcorealg/GeoAlgo/GeoAlgoException.h"

namespace geoalgo {

  Cylinder::Cylinder() : Line(), _radius(0.) {}

  Cylinder::Cylinder(double const x_min,
                     double const y_min,
                     double const z_min,
                     double const x_max,
                     double const y_max,
                     double const z_max,
                     double const radius)
    : Line(x_min, y_min, z_min, x_max, y_max, z_max), _radius(radius)
  {}

  Cylinder::Cylinder(Point_t const& min, Vector_t const& max, double const radius)
    : Line(min, max), _radius(radius)
  {
    if (min.size() != 3 || max.size() != 3)
      throw GeoAlgoException("Cylinder ctor accepts only 3D Point!");
  }

  bool Cylinder::Contain(Point_t const& pt) const
  {

    // get a vector that defines the axis of the cylinder
    Vector_t axis = _pt1 - _pt2;
    Vector_t dirpt = pt - _pt2;

    // angle of point w.r.t. the axis
    double angleMin = axis.Angle(dirpt);

    // if the angle is > 90 -> outside -> return
    if (angleMin > 0.5 * 3.14) return false;

    // revert the axis direction
    axis = _pt2 - _pt1;
    dirpt = pt - _pt1;
    angleMin = axis.Angle(dirpt);

    // if the angle is > 90 -> outside -> return
    if (angleMin > 0.5 * 3.14) return false;

    // if still here, all that is left to verify is
    // that the point isn't more than a radius
    // away from the cylinder axis
    // 1) make a line corresponding to the axis
    // 2) get the distance between the point and the line
    double radial_dist_sq = _geoAlgo.SqDist(*this, pt);

    if (radial_dist_sq > _radius * _radius) return false;

    return true;
  }
}
