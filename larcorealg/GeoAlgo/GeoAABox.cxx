#include "larcorealg/GeoAlgo/GeoAABox.h"
#include "larcorealg/GeoAlgo/GeoAlgoException.h"

namespace geoalgo {

  AABox::AABox() : _min(3), _max(3) {}

  AABox::AABox(double const x_min,
               double const y_min,
               double const z_min,
               double const x_max,
               double const y_max,
               double const z_max)
    : _min(x_min, y_min, z_min), _max(x_max, y_max, z_max)
  {}

  AABox::AABox(Point_t const& min, Vector_t const& max) : _min(min), _max(max)
  {
    if (min.size() != 3 || max.size() != 3)
      throw GeoAlgoException("AABox ctor accepts only 3D Point!");
  }

  Point_t const& AABox::Min() const { return _min; }
  Point_t const& AABox::Max() const { return _max; }

  void AABox::Min(double const x, double const y, double const z)
  {
    _min[0] = x;
    _min[1] = y;
    _min[2] = z;
  }
  void AABox::Max(double const x, double const y, double const z)
  {
    _max[0] = x;
    _max[1] = y;
    _max[2] = z;
  }

  bool AABox::Contain(Point_t const& pt) const
  {
    return !((pt[0] < _min[0] || _max[0] < pt[0]) || // point is outside X boundaries OR
             (pt[1] < _min[1] || _max[1] < pt[1]) || // point is outside Y boundaries OR
             (pt[2] < _min[2] || _max[2] < pt[2])    // point is outside Z boundaries
    );
  }
}
