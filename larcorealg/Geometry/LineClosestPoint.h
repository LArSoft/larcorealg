/**
 * @file   larcorealg/Geometry/LineClosestPoint.h
 * @brief  Utility for intersection of two 3D lines.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    larcorealg/Geometry/LineClosestPoint.tcc,
 *         larcorealg/Geometry/WireGeo.h
 *
 * This utility is used as implementation of detector wires, but it can stand on
 * its own.
 * 
 * This library is header-only and (likely) no additional linkage.
 */


#ifndef LARCOREALG_GEOMETRY_LINECLOSESTPOINT_H
#define LARCOREALG_GEOMETRY_LINECLOSESTPOINT_H


// C++ standard library
#include <utility> // std::pair<>


// -----------------------------------------------------------------------------
namespace geo {

  /**
   * @brief Returns the point of a line that is closest to a second line.
   * @tparam Point a type describing a point
   * @tparam Vector a type describing a direction (displacement vector)
   * @param refA a reference point on the first line
   * @param dirA the direction of the first line
   * @param refB a reference point on the second line
   * @param dirB the direction of the second line
   * @param[out] locOnLines pointer to additional output (see description)
   * @return the point of `A` closest to `B`
   *
   * The point of line `A` that is closest to line `B` is returned.
   * 
   * The two lines are _assumed_ not to be parallel, and when this prerequisite
   * is not met the behaviour is undefined.
   * 
   * @note This formulation is valid for lines in a Euclidean space of any
   *       dimension; the minimized distance is the Euclidean one.
   * 
   * If `locOnLines` is specified, a pair is returned with the distance of the
   * closest point from the reference points on line A (`first`) and B
   * (`second`), in units of `dirA` and `dirB` respectively.
   * 
   * 
   * Requirements
   * -------------
   * 
   * The following operations between points, vectors and scalars must be
   * defined:
   * 
   * * `Vector operator- (Point, Point)`: the difference of two points always
   *   exists and it returns a `Vector`;
   * * `Point operator+ (Point, Vector)`: translation of a point by a
   *   displacement vector;
   * * `Vector operator* (Vector, double)`: vector scaling by a real factor;
   * * `double dot(Vector, Vector)`: scalar (inner) product of two vectors.
   * 
   */
  template <typename Point, typename Vector>
  Point LineClosestPoint(
    Point const& startA, Vector const& dirA,
    Point const& startB, Vector const& dirB,
    std::pair<double, double>* locOnLines = nullptr
    );

  
  /**
   * @brief Returns the point of a line that is closest to a second line.
   * @tparam Point a type describing a point
   * @tparam UnitVector a type describing a direction (unit vector)
   * @param refA a reference point on the first line
   * @param dirA the direction of the first line (unity-normed)
   * @param refB a reference point on the second line
   * @param dirB the direction of the second line (unity-normed)
   * @param[out] locOnLines pointer to additional output (see description)
   * @return the point of `A` closest to `B`
   *
   * The point of line `A` that is closest to line `B` is returned.
   * 
   * The two lines are _assumed_ not to be parallel, and when this prerequisite
   * is not met the behaviour is undefined.
   * 
   * The two directions are _required_ to have norm `1`. While formally if this
   * prerequisite is not met the result is undefined, the result will be still
   * mostly correct if their norm departs from unity only by a rounding error.
   * The more the vectors deviate from that condition, the larger the error
   * in the result.
   *
   * @note This formulation is valid for lines in a Euclidean space of any
   *       dimension; the minimized distance is the Euclidean one.
   * 
   * If `locOnLines` is specified, a pair is returned with the distance of the
   * closest point from the reference points on line A (`first`) and B
   * (`second`), in the same space units as the points and vectors.
   * 
   * 
   * Requirements
   * -------------
   * 
   * The following operations between points, vectors and scalars must be
   * defined:
   * 
   * * `Vector operator* (UnitVector, double)`: scaling of a unit vector by a
   *   real factor to produce a non-unit vector;
   * * `Vector operator- (Point, Point)`: the difference of two points always
   *   exists and it returns a `Vector`;
   * * `Point operator+ (Point, Vector)`: translation of a point by a
   *   displacement vector;
   * * `double dot(Vector, UnitVector)`, `double dot(UnitVector, UnitVector)`:
   *   scalar (inner) product of one unit vector and a vector, or two unit
   *   vectors.
   * 
   */
  template <typename Point, typename UnitVector>
  Point LineClosestPointWithUnitVectors(
    Point const& startA, UnitVector const& dirA,
    Point const& startB, UnitVector const& dirB,
    std::pair<double, double>* locOnLines = nullptr
    );
  
  
} // namespace geo


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
#include "larcorealg/Geometry/LineClosestPoint.tcc"

// -----------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_LINECLOSESTPOINT_H
