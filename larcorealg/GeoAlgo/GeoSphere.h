/**
 * \file GeoSphere.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class HalfLine
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOSPHERE_H
#define BASICTOOL_GEOSPHERE_H

#include "larcorealg/GeoAlgo/GeoVector.h"

#include <vector>

namespace geoalgo {
  /**
     \class Spehere
     @brief Representation of a 3D sphere
     Defines a 3D Sphere having an center (Point_t) and a radius (double)
  */
  class Sphere {

  public:
    Sphere();            ///< Default ctor
    virtual ~Sphere() {} ///< Default dtor

    /// Alternative ctor (0)
    Sphere(double const& x, double const& y, double const& z, double const& r);

    /// Altenartive ctor (1) - 1 Point
    Sphere(Point_t const& center, double const r = 0);

    /// Alternative ctor (2) - 2 Points
    Sphere(Point_t const& pt1, Point_t const& pt2);

    /// Alternative ctor (3) - 3 Points
    Sphere(Point_t const& A, Point_t const& B, Point_t const& C);

    //  Alternative ctor (4) - 4 Points
    Sphere(Point_t const& A, Point_t const& B, Point_t const& C, Point_t const& D);

    // Alternative ctor (5) - Set of points
    Sphere(const std::vector<::geoalgo::Point_t>& pts);

    //
    // Getters
    //
    Point_t const& Center() const; ///< Center getter
    double Radius() const;         ///< Radius getter

    //
    // Setters
    //
    void Center(double const x, double const y, double const z); ///< Center setter
    void Center(Point_t const& pt);                              ///< Center setter
    void Radius(double const& r);                                ///< Radius setter

    //
    // Utilities
    //
    bool Contain(Point_t const& p) const; ///< Judge if a point is contained within a sphere

  protected:
    void compat(Point_t const& p, double const r = 0) const; ///< 3D point compatibility check
    void compat(double const& r) const; ///< Positive radius compatibility check

    /// Center of Sphere
    Point_t _center;

    /// Radius of Sphere
    double _radius;

  public:
    //
    // Templates
    //
    /*
#ifndef __CINT__ // Not sure why but CINT has a problem with this ctor. FIXME
    template <class T> Sphere(T const& center, double const r=0)
      : Sphere(Point_t(center),r)
    {}
#endif
    */
    template <class T>
    Sphere(T const& pt1, T const& pt2) : Sphere(Point_t(pt1), Point_t(pt2))
    {}

    template <class T>
    Sphere(T const& A, T const& B, T const& C) : Sphere(Point_t(A), Point_t(B), Point_t(C))
    {}

    template <class T>
    Sphere(T const& A, T const& B, T const& C, T const& D)
      : Sphere(Point_t(A), Point_t(B), Point_t(C), Point_t(D))
    {}

    template <class T>
    Sphere(const std::vector<T>& pts)
    {
      std::vector<::geoalgo::Vector> geo_pts;
      geo_pts.reserve(pts);
      for (auto const& p : pts)
        geo_pts.emplace_back(p);
      (*this) = Sphere(geo_pts);
    }

    template <class T>
    void Center(T const& pt)
    {
      Center(Point_t(pt));
    }

    template <class T>
    bool Contain(T const& p) const
    {
      return Contain(Point_t(p));
    }
  };

  typedef Sphere Sphere_t;
}
#endif
/** @} */ // end of doxygen group
