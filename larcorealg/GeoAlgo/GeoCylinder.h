/**
 * \file GeoCylinder.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class Cylinder
 *
 * @author david caratelli
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOCYLINDER_H
#define BASICTOOL_GEOCYLINDER_H

#include "larcorealg/GeoAlgo/GeoAlgo.h"
#include "larcorealg/GeoAlgo/GeoLine.h"
#include "larcorealg/GeoAlgo/GeoVector.h"

namespace geoalgo {
  /**
     \class Cylinder
     @brief Representation of a 3D Cylinder volume.
     A Cylinder object inherits from a geoalgo::Line
     @remark input:
     * 2 points, which define the line representing
     the central axis of the cylinder
     * a radius, defining the radius of the cylinder
  */
  class Cylinder : public Line {

  public:
    /// Default constructor
    Cylinder();

    /// Default destructor
    virtual ~Cylinder(){};

    /// Alternative ctor (0)
    Cylinder(double const x_min,
             double const y_min,
             double const z_min,
             double const x_max,
             double const y_max,
             double const z_max,
             double const radius);

    /// Altenartive ctor (1)
    Cylinder(Point_t const& min, Vector_t const& max, double const radius);

    /// Containment evaluation
    bool Contain(Point_t const& pt) const; ///< Test if a point is contained within the box

    /// Getters
    double GetRadius() { return _radius; }
    /// Setters
    void SetRadius(double r) { _radius = r; }

  protected:
    double _radius; ///< Radius of the cylinder

    // geoalgo utility
    GeoAlgo _geoAlgo;
  };

  typedef Cylinder Cylinder_t;
}
#endif
/** @} */ // end of doxygen group
