/**
 * @file   larcorealg/Geometry/LineClosestPoint.h
 * @brief  Utility for intersection of two 3D lines.
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

  /// Data structure for return values of `LineClosestPointAndOffsets()`.
  template <typename Point>
  struct IntersectionPointAndOffsets {

    Point point;    ///< Intersection point.
    double offset1; ///< Distance from reference point of first line.
    double offset2; ///< Distance from reference point of second line.

    /// Helper to assign to `std::tie()`.
    operator std::tuple<Point&, double&, double&>() noexcept { return {point, offset1, offset2}; }

  }; // IntersectionPointAndOffsets<>

  /**
   * @brief Returns the point of a line that is closest to a second line.
   * @tparam Point a type describing a point
   * @tparam Vector a type describing a direction (displacement vector)
   * @param refA a reference point on the first line
   * @param dirA the direction of the first line
   * @param refB a reference point on the second line
   * @param dirB the direction of the second line
   * @return a data structure with three fields:
   *         `point`: the point of `A` closest to `B`,
   *         `offset1`: its offset on `A` in units of `dirA`,
   *         `offset2`: its offset on `B` in units of `dirB`
   * @see `LineClosestPointWithUnitVectors()`
   * @see `LineClosestPointAndOffsets()`
   *
   * The point of line `A` that is closest to line `B` is returned.
   *
   * This function is equivalent to `LineClosestPoint()`, but it
   * returns in addition the offsets of the intersection point from the
   * reference points of the two lines, in the direction specified by
   * `dirA`/`dirB`.
   *
   * The return value is a data structure of type
   * `geo::IntersectionPointAndOffsets`, which is most easily unpacked
   * immediately:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto [ point, offsetA, offsetB ] = geo::LineClosestPointAndOffsets(
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0.0,   0.5, 0.0 },
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 0.866, 0.0, 0.0 }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will set `point` to `geo::Point{ 2, 1, 1 }`, `offsetA` to `2` and `offsetB`
   * to `2.309...`.
   * To reassign the variables after they have been defined, though, a temporary
   * structure is needed:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const xsectAndOfs = geo::LineClosestPointAndOffsets(
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 0.866, 0.0, 0.0 },
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0.0,   0.5, 0.0 }
   *   );
   * point = xsectAndOfs.point;
   * offsetA = xsectAndOfs.offset1;
   * offsetB = xsectAndOfs.offset2;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`point` to `geo::Point{ 2, 1, 0 }`, `offsetA` to `2.039...` and `offsetB`
   * to `2`, because the intersection point is always on the first line).
   *
   */
  template <typename Point, typename Vector>
  IntersectionPointAndOffsets<Point> LineClosestPointAndOffsets(Point const& startA,
                                                                Vector const& dirA,
                                                                Point const& startB,
                                                                Vector const& dirB);

  /**
   * @brief Returns the point of a line that is closest to a second line.
   * @tparam Point a type describing a point
   * @tparam Vector a type describing a direction (displacement vector)
   * @param refA a reference point on the first line
   * @param dirA the direction of the first line
   * @param refB a reference point on the second line
   * @param dirB the direction of the second line
   * @return the point of `A` closest to `B`
   * @see LineClosestPointAndOffsets(), LineClosestPointWithUnitVectors()
   *
   * The point of line `A` that is closest to line `B` is returned.
   *
   * The two lines are _assumed_ not to be parallel, and when this prerequisite
   * is not met the behaviour is undefined.
   *
   * @note This formulation is valid for lines in a Euclidean space of any
   *       dimension; the minimized distance is the Euclidean one.
   *
   * A separate function, `LineClosestPointAndOffsets()`,
   * also returns the offset of the intersection from the two reference points.
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
  Point LineClosestPoint(Point const& startA,
                         Vector const& dirA,
                         Point const& startB,
                         Vector const& dirB);

  /**
   * @brief Returns the point of a line that is closest to a second line.
   * @tparam Point a type describing a point
   * @tparam UnitVector a type describing a direction (unit vector)
   * @param refA a reference point on the first line
   * @param dirA the direction of the first line (unity-normed)
   * @param refB a reference point on the second line
   * @param dirB the direction of the second line (unity-normed)
   * @return a data structure with three fields:
   *         `point`: the point of `A` closest to `B`,
   *         `offset1`: its offset on `A` in units of `dirA` (i.e. unity),
   *         `offset2`: its offset on `B` in units of `dirB` (i.e. unity)
   * @see `LineClosestPointWithUnitVectors()`
   * @see `LineClosestPointAndOffsets()`
   *
   * The point of line `A` that is closest to line `B` is returned.
   *
   * The return value is a data structure of type
   * `geo::IntersectionPointAndOffsets`, which is most easily unpacked
   * immediately:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto [ point, offsetA, offsetB ] = geo::LineClosestPointAndOffsetsWithUnitVectors(
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0, 1, 0 },
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 1, 0, 0 }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will set `point` to `geo::Point{ 2, 1, 1 }`, `offsetA` to `1` and `offsetB`
   * to `2`.
   * To reassign the variables after they have been defined, instead:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const xsectAndOfs = geo::LineClosestPointAndOffsetsWithUnitVectors(
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 1, 0, 0 },
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0, 1, 0 }
   *   );
   * point = xsectAndOfs.point;
   * offsetA = xsectAndOfs.offset1;
   * offsetB = xsectAndOfs.offset2;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`point` to `geo::Point{ 2, 1, 0 }`, `offsetA` to `2` and `offsetB` to `1`,
   * because the intersection point is always on the first line).
   *
   */
  template <typename Point, typename UnitVector>
  IntersectionPointAndOffsets<Point> LineClosestPointAndOffsetsWithUnitVectors(
    Point const& startA,
    UnitVector const& dirA,
    Point const& startB,
    UnitVector const& dirB);

  /**
   * @brief Returns the point of a line that is closest to a second line.
   * @tparam Point a type describing a point
   * @tparam UnitVector a type describing a direction (unit vector)
   * @param refA a reference point on the first line
   * @param dirA the direction of the first line (unity-normed)
   * @param refB a reference point on the second line
   * @param dirB the direction of the second line (unity-normed)
   * @return the point of `A` closest to `B`
   * @see LineClosestPointAndOffsetsWithUnitVectors(), LineClosestPoint()
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
   * A separate function, `LineClosestPointAndOffsetsWithUnitVectors()`,
   * also returne the offset of the intersection from the two reference points.
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
  Point LineClosestPointWithUnitVectors(Point const& startA,
                                        UnitVector const& dirA,
                                        Point const& startB,
                                        UnitVector const& dirB);

} // namespace geo

// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
#include "larcorealg/Geometry/LineClosestPoint.tcc"

// -----------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_LINECLOSESTPOINT_H
