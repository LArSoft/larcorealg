/**
 * \file GeoLine.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class Line
 *
 * @author David Caratelli
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOLINE_H
#define BASICTOOL_GEOLINE_H

#include "larcorealg/GeoAlgo/GeoVector.h"

namespace geoalgo {
  /**
     \class Line
     @brief Representation of a 3D infinite line.
     Defines an infinite 3D line by having 2 points which completely determine the line
     along which the line extends. It hides the point attributes from users for   \n
     protecting the dimensionality.
  */
  class Line {

  public:
    /// Default constructor
    Line();

    /// Default destructor
    virtual ~Line() {}

    /// Alternative ctor (1)
    Line(double const x1,
         double const y1,
         double const z1,
         double const x2,
         double const y2,
         double const z2);

    /// Altenartive ctor (2)
    Line(Point_t const& pt1, Point_t const& pt2);

    //
    // Getters
    //
    Point_t const& Pt1() const; ///< Start getter
    Point_t const& Pt2() const; ///< Direction getter

    //
    // Setters
    //
    void Pt1(double const x, double const y, double const z); ///< Pt1 setter
    void Pt2(double const x, double const y, double const z); ///< Pt2 setter

  protected:
    /// Compatibility check
    void check_and_raise(Point_t const& p1, Point_t const& p2) const;

    Point_t _pt1;  ///< First point denoting infinite line
    Vector_t _pt2; ///< Second point denoting infinite line

  public:
    //
    // Template
    //
    /// Alternative ctor using template (3)
    template <class T, class U>
    Line(T const& pt1, U const& pt2) : Line(Point_t(pt1), Point_t(pt2))
    {}

    /// Pt1 setter template
    template <class T>
    void Pt1(T const& pt1)
    {
      _pt1 = Point_t(pt1);
      check_and_raise(_pt1, _pt2);
    }

    /// Pt2 setter template
    template <class T>
    void Pt2(T const& pt2)
    {
      _pt2 = Vector_t(pt2);
      check_and_raise(_pt1, _pt2);
    }
  };

  typedef Line Line_t;
}
#endif
/** @} */ // end of doxygen group
