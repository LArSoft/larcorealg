/**
 * \file GeoAABox.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class AABox
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOAABOX_H
#define BASICTOOL_GEOAABOX_H

#include "larcorealg/GeoAlgo/GeoVector.h" // for Point_t, Vector_t

namespace geoalgo {
  /**
     \class AABox
     @brief Representation of a 3D rectangular box which sides are aligned w/ coordinate axis.
     A representation of an Axis-Aligned-Boundary-Box, a simple & popular representation of   \n
     3D boundary box for collision detection. The concept was taken from the reference,       \n
     Real-Time-Collision-Detection (RTCD), and in particular Ch. 4.2 (page 77):               \n

     Ref: http://realtimecollisiondetection.net

     This class uses one of the simplest representation for AABox: "min-max" representation.   \n
     Though this method requires storing 6 floating point values, class attributes (i.e.      \n
     "min" and "max" points) store intuitive values for most UB analyzers. Also it simplifies \n
     utility function implementations.
  */
  class AABox {

  public:
    /// Default constructor
    AABox();

    /// Default destructor
    virtual ~AABox(){};

    /// Alternative ctor (0)
    AABox(double const x_min,
          double const y_min,
          double const z_min,
          double const x_max,
          double const y_max,
          double const z_max);

    /// Altenartive ctor (1)
    AABox(Point_t const& min, Vector_t const& max);

    //
    // Attribute accessor
    //
    Point_t const& Min() const;                               ///< Minimum point getter
    Point_t const& Max() const;                               ///< Maximum point getter
    void Min(double const x, double const y, double const z); ///< Minimum point setter
    void Max(double const x, double const y, double const z); ///< Maximum point setter
    bool Contain(Point_t const& pt) const; ///< Test if a point is contained within the box

  protected:
    Point_t _min; ///< Minimum point
    Point_t _max; ///< Maximum point

  public:
    //
    // Template
    //
    /// Alternative ctor using template (3)
    template <class T, class U>
    AABox(T const& min, U const& max) : AABox(Point_t(min), Point_t(max))
    {}
  };

  typedef AABox AABox_t;
}
#endif
/** @} */ // end of doxygen group
