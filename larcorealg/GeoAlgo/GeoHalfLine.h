/**
 * \file GeoHalfLine.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class HalfLine
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOHALFLINE_H
#define BASICTOOL_GEOHALFLINE_H

#include "larcorealg/GeoAlgo/GeoAlgoException.h"
#include "larcorealg/GeoAlgo/GeoVector.h"

#include "TVector3.h"

namespace geoalgo {
  /**
     \class HalfLine
     @brief Representation of a 3D semi-infinite line.
     Defines a semi-infinite 3D line by having a start point (Point_t) and direction (Vector_t) \n
     along which the line extends. It hides the start and direction attributes from users for   \n
     protecting the dimensionality.
  */
  class HalfLine {

  public:
    /// Default constructor
    HalfLine();

    /// Default destructor
    virtual ~HalfLine(){};

    /// Alternative ctor (1)
    HalfLine(double const x,
             double const y,
             double const z,
             double const dirx,
             double const diry,
             double const dirz);

    /// Altenartive ctor (2)
    HalfLine(Point_t const& start, Vector_t const& dir);

    Point_t const& Start() const; ///< Start getter
    Vector_t const& Dir() const;  ///< Direction getter

    void Start(double const x, double const y, double const z); ///< Start setter
    void Dir(double const x, double const y, double const z);   ///< Dir setter

    void Start(TVector3 const& pt); ///< Start setter
    void Dir(TVector3 const& dir);  ///< Dir setter

  protected:
    void Normalize(); ///< Normalize direction
    Point_t _start;   ///< Beginning of the half line
    Vector_t _dir;    ///< Direction of the half line from _start

  public:
    //
    // Template
    //

    /// Alternative ctor using template (3)
    template <class T, class U>
    HalfLine(T const& start, U const& dir) : HalfLine(Point_t(start), Vector_t(dir))
    {}

    /// Start setter template
    template <class T>
    void Start(T const& pos)
    {
      _start = Point_t(pos);
      if (_start.size() != 3)
        throw GeoAlgoException("<<Start>> Only 3 dimensional start point allowed!");
    }

    /// Dir setter template
    template <class T>
    void Dir(T const& dir)
    {
      _dir = Vector_t(dir);
      if (_dir.size() != 3)
        throw GeoAlgoException("<<Start>> Only 3 dimensional start point allowed!");
      Normalize();
    }
  };

  typedef HalfLine HalfLine_t;
}
#endif
/** @} */ // end of doxygen group
