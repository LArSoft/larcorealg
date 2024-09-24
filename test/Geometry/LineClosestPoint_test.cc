/**
 * @file   LineClosestPoint_test.cc
 * @brief  Test of `LineClosestPoint.h` utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 29, 2021
 * @see    `larcorealg/Geometry/LineClosestPoint.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE LineClosestPoint_test
#include <boost/test/unit_test.hpp> // BOOST_AUTO_TEST_CASE(), BOOST_TEST()

// LArSoft libraries
#include "larcorealg/Geometry/LineClosestPoint.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// C++ standard library
#include <cmath>       // std::sqrt()
#include <type_traits> // std::is_same_v<>
#include <utility>     // std::pair<>

using namespace geo;

namespace {

  void LineClosestPointSimple_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    Point_t const p = LineClosestPoint(
      origin() - Zaxis() - 3.0 * Xaxis(), Xaxis(), origin() + Zaxis() + 2.0 * Yaxis(), Yaxis());

    BOOST_TEST(p.X() == 0.0, tol);
    BOOST_TEST(p.Y() == 0.0, tol);
    BOOST_TEST(p.Z() == -1.0, tol);
  }

  // -----------------------------------------------------------------------------
  void LineClosestPointSimple45_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    Point_t const p =
      LineClosestPoint(origin(), Xaxis(), origin() + Yaxis(), (Xaxis() + Zaxis()) / std::sqrt(2.0));

    BOOST_TEST(p.X() == 0.0, tol);
    BOOST_TEST(p.Y() == 0.0, tol);
    BOOST_TEST(p.Z() == 0.0, tol);
  }

  // -----------------------------------------------------------------------------
  void LineClosestPointAndOffsets_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    auto const [p, ofsA, ofsB] = LineClosestPointAndOffsets(
      origin() - 3.0 * Xaxis(), Xaxis(), origin() + Zaxis() + 2.0 * Yaxis(), Yaxis());

    BOOST_TEST(p.X() == 0.0, tol);
    BOOST_TEST(p.Y() == 0.0, tol);
    BOOST_TEST(p.Z() == 0.0, tol);
    BOOST_TEST(ofsA == +3.0, tol);
    BOOST_TEST(ofsB == -2.0, tol);
  }

  // -----------------------------------------------------------------------------
  void LineClosestPointWithScaledDirs_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    auto const [p, ofsA, ofsB] = LineClosestPointAndOffsets(
      origin() - 3.0 * Xaxis(), 2.0 * Xaxis(), origin() + Zaxis() + 2.0 * Yaxis(), -2.0 * Yaxis());

    BOOST_TEST(p.X() == 0.0, tol);
    BOOST_TEST(p.Y() == 0.0, tol);
    BOOST_TEST(p.Z() == 0.0, tol);
    BOOST_TEST(ofsA == (+3.0 / 2.0), tol);
    BOOST_TEST(ofsB == (-2.0 / -2.0), tol);
  }

  // -----------------------------------------------------------------------------
  void LineClosestPointWithNonHomogeneousDirs_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    auto const [p, ofsA, ofsB] = LineClosestPointAndOffsets(
      origin() - 3.0 * Xaxis(), 1.5 * Xaxis(), origin() + Zaxis() + 2.0 * Yaxis(), -2.0 * Yaxis());

    BOOST_TEST(p.X() == 0.0, tol);
    BOOST_TEST(p.Y() == 0.0, tol);
    BOOST_TEST(p.Z() == 0.0, tol);
    BOOST_TEST(ofsA == (+3.0 / 1.5), tol);
    BOOST_TEST(ofsB == (-2.0 / -2.0), tol);
  }

  // -----------------------------------------------------------------------------
  void LineClosestPointAndOffsetsDocumentation_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    /*
   * The promise:
   *
   * --- 8< --------------------------------------------------------------------
   * The return value is a triplet, which is most easily unpacked immediately:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto [ point, offsetA, offsetB ] = LineClosestPointAndOffsets(
   *   Point_t{ 2, 0, 1 }, Vector_t{ 0.0,   0.5, 0.0 },
   *   Point_t{ 0, 1, 0 }, Vector_t{ 0.866, 0.0, 0.0 }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will set `point` to `Point{ 2, 1, 1 }`, `offsetA` to `2` and `offsetB`
   * to `2.309...`.
   * To reassign the variables after they have been defined, instead:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const xsectAndOfs = LineClosestPointAndOffsets(
   *   Point_t{ 0, 1, 0 }, Vector_t{ 0.866, 0.0, 0.0 },
   *   Point_t{ 2, 0, 1 }, Vector_t{ 0.0,   0.5, 0.0 }
   *   );
   * point = xsectAndOfs.point;
   * offsetA = xsectAndOfs.offset1;
   * offsetB = xsectAndOfs.offset2;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`point` to `Point{ 2, 1, 0 }`, `offsetA` to `2.039...` and `offsetB`
   * to `2`, because the intersection point is always on the first line).
   * --- 8< --------------------------------------------------------------------
   */

    auto [point, offsetA, offsetB] = LineClosestPointAndOffsets(
      Point_t{2, 0, 1}, Vector_t{0.0, 0.5, 0.0}, Point_t{0, 1, 0}, Vector_t{0.866, 0.0, 0.0});

    // a way to check we did not mess with the assignment above too much
    static_assert(std::is_same_v<decltype(point), Point_t>, "Unexpected point type");
    static_assert(std::is_same_v<decltype(offsetA), double>, "Unexpected first offset type");
    static_assert(std::is_same_v<decltype(offsetB), double>, "Unexpected second offset type");

    BOOST_TEST(point.X() == 2.0, tol);
    BOOST_TEST(point.Y() == 1.0, tol);
    BOOST_TEST(point.Z() == 1.0, tol);
    BOOST_TEST(offsetA == 2.0, tol);
    BOOST_TEST(offsetB == (2.0 / 0.866), tol);

    auto const xsectAndOfs = LineClosestPointAndOffsets(
      Point_t{0, 1, 0}, Vector_t{0.866, 0.0, 0.0}, Point_t{2, 0, 1}, Vector_t{0.0, 0.5, 0.0});
    point = xsectAndOfs.point;
    offsetA = xsectAndOfs.offset1;
    offsetB = xsectAndOfs.offset2;

    BOOST_TEST(point.X() == 2.0, tol);
    BOOST_TEST(point.Y() == 1.0, tol);
    BOOST_TEST(point.Z() == 0.0, tol);
    BOOST_TEST(offsetA == 2.0 / 0.866, tol);
    BOOST_TEST(offsetB == 2.0, tol);

    // actually we _can_ assign with `std::tie()`:
    std::tie(point, offsetA, offsetB) = LineClosestPointAndOffsets(
      Point_t{2, 0, 1}, Vector_t{0.0, 0.5, 0.0}, Point_t{0, 1, 0}, Vector_t{0.866, 0.0, 0.0});

    BOOST_TEST(point.X() == 2.0, tol);
    BOOST_TEST(point.Y() == 1.0, tol);
    BOOST_TEST(point.Z() == 1.0, tol);
    BOOST_TEST(offsetA == 2.0, tol);
    BOOST_TEST(offsetB == (2.0 / 0.866), tol);
  }

  void LineClosestPointWithUnitVectorsSimple_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    Point_t const p = LineClosestPointWithUnitVectors(
      origin() - Zaxis() - 3.0 * Xaxis(), Xaxis(), origin() + Zaxis() + 2.0 * Yaxis(), Yaxis());

    BOOST_TEST(p.X() == 0.0, tol);
    BOOST_TEST(p.Y() == 0.0, tol);
    BOOST_TEST(p.Z() == -1.0, tol);
  }

  // -----------------------------------------------------------------------------
  void LineClosestPointWithUnitVectorsSimple45_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    Point_t const p = LineClosestPointWithUnitVectors(
      origin(), Xaxis(), origin() + Yaxis(), (Xaxis() + Zaxis()) / std::sqrt(2.0));

    BOOST_TEST(p.X() == 0.0, tol);
    BOOST_TEST(p.Y() == 0.0, tol);
    BOOST_TEST(p.Z() == 0.0, tol);
  }

  // -----------------------------------------------------------------------------
  void LineClosestPointWithUnitVectorsAndOffsets_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    auto const [p, ofsA, ofsB] = LineClosestPointAndOffsetsWithUnitVectors(
      origin() - 3.0 * Xaxis(), Xaxis(), origin() + Zaxis() + 2.0 * Yaxis(), Yaxis());

    BOOST_TEST(p.X() == 0.0, tol);
    BOOST_TEST(p.Y() == 0.0, tol);
    BOOST_TEST(p.Z() == 0.0, tol);
    BOOST_TEST(ofsA == +3.0, tol);
    BOOST_TEST(ofsB == -2.0, tol);
  }

  // -----------------------------------------------------------------------------
  void LineClosestPointAndOffsetsWithUnitVectorsDocumentation_test()
  {
    auto const tol = boost::test_tools::tolerance(0.001);

    /*
   * The promise:
   *
   * --- 8< --------------------------------------------------------------------
   * The return value is a triplet, which is most easily unpacked immediately:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto [ point, offsetA, offsetB ] = LineClosestPointAndOffsetsWithUnitVectors(
   *   Point_t{ 2, 0, 1 }, Vector_t{ 0, 1, 0 },
   *   Point_t{ 0, 1, 0 }, Vector_t{ 1, 0, 0 }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will set `point` to `Point{ 2, 1, 1 }`, `offsetA` to `1` and `offsetB`
   * to `2`.
   * To reassign the variables after they have been defined, instead:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const xsectAndOfs = LineClosestPointAndOffsetsWithUnitVectors(
   *   Point_t{ 0, 1, 0 }, Vector_t{ 1, 0, 0 },
   *   Point_t{ 2, 0, 1 }, Vector_t{ 0, 1, 0 }
   *   );
   * point = xsectAndOfs.point;
   * offsetA = xsectAndOfs.offset1;
   * offsetB = xsectAndOfs.offset2;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`point` to `Point{ 2, 1, 0 }`, `offsetA` to `2` and `offsetB` to `1`,
   * because the intersection point is always on the first line).
   * --- 8< --------------------------------------------------------------------
   */

    auto [point, offsetA, offsetB] = LineClosestPointAndOffsetsWithUnitVectors(
      Point_t{2, 0, 1}, Vector_t{0, 1, 0}, Point_t{0, 1, 0}, Vector_t{1, 0, 0});

    // a way to check we did not mess with the assignment above too much
    static_assert(std::is_same_v<decltype(point), Point_t>, "Unexpected point type");
    static_assert(std::is_same_v<decltype(offsetA), double>, "Unexpected first offset type");
    static_assert(std::is_same_v<decltype(offsetB), double>, "Unexpected second offset type");

    BOOST_TEST(point.X() == 2.0, tol);
    BOOST_TEST(point.Y() == 1.0, tol);
    BOOST_TEST(point.Z() == 1.0, tol);
    BOOST_TEST(offsetA == 1.0, tol);
    BOOST_TEST(offsetB == 2.0, tol);

    auto const xsectAndOfs = LineClosestPointAndOffsetsWithUnitVectors(
      Point_t{0, 1, 0}, Vector_t{1, 0, 0}, Point_t{2, 0, 1}, Vector_t{0, 1, 0});
    point = xsectAndOfs.point;
    offsetA = xsectAndOfs.offset1;
    offsetB = xsectAndOfs.offset2;

    BOOST_TEST(point.X() == 2.0, tol);
    BOOST_TEST(point.Y() == 1.0, tol);
    BOOST_TEST(point.Z() == 0.0, tol);
    BOOST_TEST(offsetA == 2.0, tol);
    BOOST_TEST(offsetB == 1.0, tol);

    // actually we _can_ assign with `std::tie()`:
    std::tie(point, offsetA, offsetB) = LineClosestPointAndOffsetsWithUnitVectors(
      Point_t{2, 0, 1}, Vector_t{0, 1, 0}, Point_t{0, 1, 0}, Vector_t{1, 0, 0});

    BOOST_TEST(point.X() == 2.0, tol);
    BOOST_TEST(point.Y() == 1.0, tol);
    BOOST_TEST(point.Z() == 1.0, tol);
    BOOST_TEST(offsetA == 1.0, tol);
    BOOST_TEST(offsetB == 2.0, tol);
  }
} // anon. namespace

// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(LineClosestPointTestCase)
{
  LineClosestPointSimple_test();
  LineClosestPointSimple45_test();
  LineClosestPointAndOffsets_test();
  LineClosestPointWithScaledDirs_test();
  LineClosestPointWithNonHomogeneousDirs_test();
  LineClosestPointAndOffsetsDocumentation_test();
}

// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(LineClosestPointWithUnitVectorsTestCase)
{
  LineClosestPointWithUnitVectorsSimple_test();
  LineClosestPointWithUnitVectorsSimple45_test();
  LineClosestPointWithUnitVectorsAndOffsets_test();
  LineClosestPointAndOffsetsWithUnitVectorsDocumentation_test();
}
