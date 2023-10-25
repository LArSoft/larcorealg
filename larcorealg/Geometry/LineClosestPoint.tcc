/**
 * @file   larcorealg/Geometry/LineClosestPoint.tcc
 * @brief  Utility for intersection of two 3D lines: template implementation.
 * @see    larcorealg/Geometry/LineClosestPoint.h
 */

#ifndef LARCOREALG_GEOMETRY_LINECLOSESTPOINT_H
#error                                                                                             \
  "LineClosestPoint.tcc should not be included directly: #include \"larcorealg/Geometry/LineClosestPoint.h\" instead."
#endif

// LArSoft and framework libraries
#include "cetlib/pow.h" // cet::square

// C/C++ libraries
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::dot()

// C++ standard library
#include <cassert>
#include <cmath> // std::abs()
#include <tuple> // std::tie()

// -----------------------------------------------------------------------------
// ---  implementation details
// -----------------------------------------------------------------------------
namespace geo::details {

  template <typename Point, typename Vector>
  Point LineClosestPointImpl(Point const& startA,
                             Vector const& dirA,
                             Point const& startB,
                             Vector const& dirB,
                             std::pair<double, double>* locOnLines /* = nullptr */
  )
  {
    /*
     * The point on the first line ("A"):
     *
     *     p1(t) = c1 + t w1 (c1 the starting point, w1 its direction)
     *
     * has the minimal distance from the other line (c2 + u w2, with the same
     * notation) at
     *
     *     t = [-(dc,w1)(w2,w2) + (dc,w2)(w1,w2)] / [ (w1,w2)^2 - (w1,w1)(w2,w2) ]
     *     u = [ (dc,w2)(w2,w2) - (dc,w1)(w1,w2)] / [ (w1,w2)^2 - (w1,w1)(w2,w2) ]
     *
     * (where (a,b) is a scalar product and dc = (c1 - c2) ).
     * If w1 and w2 are unit vectors, t and u are in fact the distance of the
     * point from the "start" of the respective lines in "standard" geometry
     * units.
     *
     * If `locOnLines` is non-null, it is filled with { t, u }.
     */

    // aliases for quick notation
    auto const& [c1, w1] = std::tie(startA, dirA);
    auto const& [c2, w2] = std::tie(startB, dirB);

    auto const dc = c2 - c1;

    using geo::vect::dot;
    double const dcw1 = dot(dc, w1);
    double const dcw2 = dot(dc, w2);
    double const w1w2 = dot(w1, w2);                 // this is cos(angle), angle between lines
    assert(std::abs(std::abs(w1w2) - 1.0) >= 1e-10); // prerequisite: not parallel

    using geo::vect::mag2;
    double const w1w1 = mag2(w1);
    double const w2w2 = mag2(w2);
    double const inv_den = 1.0 / (cet::square(w1w2) - (w1w1 * w2w2));
    double const t = ((dcw2 * w1w2) - (dcw1 * w2w2)) * inv_den;

    if (locOnLines) {
      double const u = ((dcw2 * w1w1) - (dcw1 * w1w2)) * inv_den;
      *locOnLines = {t, u};
    }

    return c1 + w1 * t;

  } // geo::details::LineClosestPointImpl()

  // ---------------------------------------------------------------------------
  template <typename Point, typename UnitVector>
  Point LineClosestPointWithUnitVectorsImpl(Point const& startA,
                                            UnitVector const& dirA,
                                            Point const& startB,
                                            UnitVector const& dirB,
                                            std::pair<double, double>* locOnLines /* = nullptr */
  )
  {
    /*
     * The implementation is the same as in `LineClosestPoint()`,
     * with some computation skipped while relying on the assumption of unit
     * vectors.
     */

    using geo::vect::mag2;
    assert(std::abs(mag2(dirA) - 1.0) < 1e-10); // prerequisite: dirA unit vector
    assert(std::abs(mag2(dirB) - 1.0) < 1e-10); // prerequisite: dirB unit vector

    // aliases for quick notation
    Point const& c1 = startA;
    UnitVector const& w1 = dirA;

    Point const& c2 = startB;
    UnitVector const& w2 = dirB;

    auto const dc = c2 - c1;

    using geo::vect::dot;
    double const dcw1 = dot(dc, w1);
    double const dcw2 = dot(dc, w2);
    double const w1w2 = dot(w1, w2);                 // this is cos(angle), angle between lines
    assert(std::abs(std::abs(w1w2) - 1.0) >= 1e-10); // prerequisite: not parallel

    // we keep the generic math here, but freeze the modulus parameters
    constexpr double const w1w1{1.0};
    constexpr double const w2w2{1.0};
    double const inv_den = 1.0 / (cet::square(w1w2) - (w1w1 * w2w2));
    double const t = ((dcw2 * w1w2) - (dcw1 * w2w2)) * inv_den;

    if (locOnLines) {
      double const u = ((dcw2 * w1w1) - (dcw1 * w1w2)) * inv_den;
      *locOnLines = {t, u};
    }

    return c1 + w1 * t;

  } // geo::details::LineClosestPointWithUnitVectorsImpl()

  // ---------------------------------------------------------------------------

} // namespace geo::details

// -----------------------------------------------------------------------------
// ---  interface
// -----------------------------------------------------------------------------
template <typename Point, typename Vector>
geo::IntersectionPointAndOffsets<Point> geo::LineClosestPointAndOffsets(Point const& startA,
                                                                        Vector const& dirA,
                                                                        Point const& startB,
                                                                        Vector const& dirB)
{
  std::pair<double, double> locOnLines;
  auto const point = details::LineClosestPointImpl(startA, dirA, startB, dirB, &locOnLines);
  return {point, locOnLines.first, locOnLines.second};
} // geo::LineClosestPointAndOffsets()

// -----------------------------------------------------------------------------
template <typename Point, typename Vector>
Point geo::LineClosestPoint(Point const& startA,
                            Vector const& dirA,
                            Point const& startB,
                            Vector const& dirB)
{
  return details::LineClosestPointImpl(startA, dirA, startB, dirB, nullptr);
} // geo::LineClosestPoint()

// -----------------------------------------------------------------------------
template <typename Point, typename UnitVector>
geo::IntersectionPointAndOffsets<Point> geo::LineClosestPointAndOffsetsWithUnitVectors(
  Point const& startA,
  UnitVector const& dirA,
  Point const& startB,
  UnitVector const& dirB)
{

  std::pair<double, double> locOnLines;
  auto const point =
    details::LineClosestPointWithUnitVectorsImpl(startA, dirA, startB, dirB, &locOnLines);
  return {point, locOnLines.first, locOnLines.second};
}

// -----------------------------------------------------------------------------
template <typename Point, typename UnitVector>
Point geo::LineClosestPointWithUnitVectors(Point const& startA,
                                           UnitVector const& dirA,
                                           Point const& startB,
                                           UnitVector const& dirB)
{
  return details::LineClosestPointWithUnitVectorsImpl(startA, dirA, startB, dirB, nullptr);
}

// -----------------------------------------------------------------------------
