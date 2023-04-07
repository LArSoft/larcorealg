/**
 * \file GeoCone.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class HalfLine
 *
 * @author David Caratelli
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOCONE_H
#define BASICTOOL_GEOCONE_H

#include "larcorealg/GeoAlgo/GeoHalfLine.h"
#include "larcorealg/GeoAlgo/GeoVector.h"
namespace geoalgo {
  /**
     \class Cone
     @brief Representation of a 3D semi-infinite line.
     Defines a 3D cone with the following properties:                                           \n
     Start point (or vertex), Direction, Length (or Length), Radius, opening angle              \n
     When 2 of Length, Radius, opening angle are defined the third is automatically set
  */
  class Cone : public HalfLine {

  public:
    /// Default constructor
    Cone();

    /// Default destructor
    virtual ~Cone(){};

    /// Alternative ctor (1)
    Cone(double const x,
         double const y,
         double const z,
         double const dirx,
         double const diry,
         double const dirz,
         double const length,
         double const radius);

    /// Alternative ctor (2)
    Cone(Point_t const& start, Vector_t const& dir, double const length, double const radius);

    //
    // Getters
    //
    double Length() const; ///< Length getter
    double Radius() const; ///< Length getter
    double Angle() const;  ///< Angle getter

    //
    // Setters
    //
    void Length(double const l); ///< Length setter
    void Radius(double const r); ///< Radius setter

  protected:
    double _length; ///< Helight (length) of the cone
    double _radius; ///< Radius of the cone at the base
    double _angle;  ///< Opening Angle

  public:
    //
    // Template
    //
    /// Alternative ctor using template (3)
    template <class T, class U>
    Cone(T const& start, U const& dir) : Cone(Point_t(start), Vector_t(dir))
    {}
  };

  typedef Cone Cone_t;
}
#endif
/** @} */ // end of doxygen group
