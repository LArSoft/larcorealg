/**
 * \file GeoLineSegment.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class LineSegment
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOLINESEGMENT_H
#define BASICTOOL_GEOLINESEGMENT_H

#include "larcorealg/GeoAlgo/GeoVector.h"

namespace geoalgo {
  /**
     \class LineSegment
     @brief Representation of a simple 3D line segment
     Defines a finite 3D straight line by having the start and end position (Point_t). \n
  */
  class LineSegment {

  public:
    /// Default constructor
    LineSegment();

    /// Default destructor
    virtual ~LineSegment() {}

    /// Alternative ctor (1)
    LineSegment(double const start_x,
                double const start_y,
                double const start_z,
                double const end_x,
                double const end_y,
                double const end_z);

    /// Altenartive ctor (2)
    LineSegment(Point_t const& start, Point_t const& end);

    //
    // Getters
    //
    Point_t const& Start() const; ///< Start getter
    Point_t const& End() const;   ///< End getter
    Vector_t const Dir() const;   ///< Direction getter

    //
    // Setters
    //
    void Start(double const x, double const y, double const z); ///< Start setter
    void End(double const x, double const y, double const z);   ///< End setter

  protected:
    void DirReset(); ///< Internal function to reset direction
    Point_t _start;  ///< Start position of a line
    Point_t _end;    ///< End position of a line
    Vector_t _dir;   ///< Direction

  public:
    //
    // Template
    //
    /// Alternative ctor using template (3)
    template <class T, class U>
    LineSegment(T const& start, U const& end) : LineSegment(Point_t(start), Point_t(end))
    {}
  };

  typedef LineSegment LineSegment_t;
}

#endif
/** @} */ // end of doxygen group
