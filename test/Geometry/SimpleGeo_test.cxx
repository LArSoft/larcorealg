/**
 * @file   SimpleGeo_test.cxx
 * @brief  Unit test for SimpleGeo library
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 11, 2016
 *
 * Usage: just run the executable.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate
// a main() function; Boost is pulled in by boost_unit_test_base.h
#define BOOST_TEST_MODULE SimpleGeoTest

// LArSoft libraries
#include "larcorealg/Geometry/SimpleGeo.h"

// Boost libraries
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <array>

//------------------------------------------------------------------------------
namespace {
  void Point2DTest()
  {
    using Point_t = lar::util::simple_geo::Point2D<float>;

    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    std::array<Point_t::Data_t, 2U> buffer = {{1., 2.}};

    // default constructor
    Point_t p1;
    BOOST_TEST(p1.x == 0.);
    BOOST_TEST(p1.y == 0.);

    // element constructor
    Point_t p2(1., 2.);
    BOOST_TEST(p2.x == 1.);
    BOOST_TEST(p2.y == 2.);

    // data buffer constructor
    Point_t p3(buffer.data());
    BOOST_TEST(p3.x == 1.);
    BOOST_TEST(p3.y == 2.);

    // implicit comparisons
    BOOST_TEST(p2 == p3);
    BOOST_TEST(p1 != p2);

    // arithmetic operations
    Point_t p4{2., 4.};
    BOOST_TEST(p1 + p2 == p2);
    BOOST_TEST(p2 + p1 == p2);
    BOOST_TEST(p1 * 2. == p1);
    BOOST_TEST(p2 * 2.0 == p4);
    BOOST_TEST(p2 / 0.5 == p4);

  } // Point2DTest()

  void Point3DTest()
  {

    using Point_t = lar::util::simple_geo::Point3D<float>;

    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    std::array<Point_t::Data_t, 3U> buffer = {{1., 2., 4.}};

    // default constructor
    Point_t p1;
    BOOST_TEST(p1.x == 0.);
    BOOST_TEST(p1.y == 0.);
    BOOST_TEST(p1.z == 0.);

    // element constructor
    Point_t p2(1., 2., 4.);
    BOOST_TEST(p2.x == 1.);
    BOOST_TEST(p2.y == 2.);
    BOOST_TEST(p2.z == 4.);

    // data buffer constructor
    Point_t p3(buffer.data());
    BOOST_TEST(p3.x == 1.);
    BOOST_TEST(p3.y == 2.);
    BOOST_TEST(p2.z == 4.);

    // implicit comparisons
    BOOST_TEST(p2 == p3);
    BOOST_TEST(p1 != p2);

    // arithmetic operations
    Point_t p4{2., 4., 8.};
    BOOST_TEST(p1 + p2 == p2);
    BOOST_TEST(p2 + p1 == p2);
    BOOST_TEST(p1 * 2. == p1);
    BOOST_TEST(p2 * 2.0 == p4);
    BOOST_TEST(p2 / 0.5 == p4);

  } // Point3DTest()

  void AreaTest()
  {
    using Point_t = lar::util::simple_geo::Point2D<float>;
    using Area_t = lar::util::simple_geo::Area<Point_t>;

    // default constructor
    Area_t A1;
    BOOST_TEST(A1.Min() == Area_t::Point_t{});
    BOOST_TEST(A1.Max() == Area_t::Point_t{});

    // point constructor
    Area_t A2({1., 6.}, {4., 2.});
    BOOST_TEST(A2.Min() == Area_t::Point_t(1., 2.));
    BOOST_TEST(A2.Max() == Area_t::Point_t(4., 6.));

    Area_t A3({1., 2.}, {2., 2.});

    // queries
    BOOST_TEST(A2.Center() == Area_t::Point_t(2.5, 4.0));

    BOOST_TEST(A1.DeltaX() == 0.0);
    BOOST_TEST(A2.DeltaX() == 3.0);

    BOOST_TEST(!A1.isNullX());
    BOOST_TEST(!A2.isNullX());
    BOOST_TEST(!A3.isNullX());

    BOOST_TEST(A1.isEmptyX());
    BOOST_TEST(!A2.isEmptyX());
    BOOST_TEST(!A3.isEmptyX());

    BOOST_TEST(A1.nonEmptyDims() == 0U);
    BOOST_TEST(A2.nonEmptyDims() == 2U);
    BOOST_TEST(A3.nonEmptyDims() == 1U);

    BOOST_TEST(!A1.isNull());
    BOOST_TEST(!A2.isNull());
    BOOST_TEST(!A3.isNull());

    BOOST_TEST(A1.isEmpty());
    BOOST_TEST(!A2.isEmpty());
    BOOST_TEST(!A3.isEmpty());

    BOOST_TEST(!A1.isLine());
    BOOST_TEST(!A2.isLine());
    BOOST_TEST(A3.isLine());

    BOOST_TEST(A1.thinnestSize() == 0.0);
    BOOST_TEST(A2.thinnestSize() == 3.0);
    BOOST_TEST(A3.thinnestSize() == 0.0);

    BOOST_TEST(((A1.thinnestSide() == 0) || (A1.thinnestSide() == 1)));
    BOOST_TEST(A2.thinnestSide() == 0U);
    BOOST_TEST(A3.thinnestSide() == 1U);

    BOOST_TEST(A1.DeltaY() == 0.0);
    BOOST_TEST(A2.DeltaY() == 4.0);
    BOOST_TEST(A3.DeltaY() == 0.0);

    BOOST_TEST(!A1.isNullY());
    BOOST_TEST(!A2.isNullY());
    BOOST_TEST(!A3.isNullY());

    BOOST_TEST(A1.isEmptyY());
    BOOST_TEST(!A2.isEmptyY());
    BOOST_TEST(A3.isEmptyY());

    BOOST_TEST(!A1.isPlane());
    BOOST_TEST(A2.isPlane());
    BOOST_TEST(!A3.isPlane());

    // modifications
    BOOST_TEST(A3 == Area_t({1.0, 2.0}, {2.0, 2.0}));
    BOOST_TEST(A3 != A2);
    A3.Intersect(A2);
    BOOST_TEST(A3 == Area_t({1.0, 2.0}, {2.0, 2.0}));

    A3.IncludePoint({1.5, 6.0});
    BOOST_TEST(A3 == Area_t({1.0, 2.0}, {2.0, 6.0}));
    A3.IncludePoint({4.0, 5.0});
    BOOST_TEST(A3 == Area_t({1.0, 2.0}, {4.0, 6.0}));
    BOOST_TEST(A3 == A2);

    A3.Intersect(A2);
    BOOST_TEST(A3 == Area_t({1.0, 2.0}, {4.0, 6.0}));

    A3.IncludePoint({0.0, 4.0});
    BOOST_TEST(A3 == Area_t({0.0, 2.0}, {4.0, 6.0}));

    A3.Include(Area_t({1.0, -1.0}, {2.0, 7.0}));
    BOOST_TEST(A3 == Area_t({0.0, -1.0}, {4.0, 7.0}));

    A3.Include(Area_t({-1.0, 0.0}, {3.0, 1.0}));
    BOOST_TEST(A3 == Area_t({-1.0, -1.0}, {4.0, 7.0}));

    // intersection result is empty
    BOOST_CHECK_THROW(A3.Intersect(Area_t({8.0, 9.0}, {8.0, 9.0})), Area_t::NullIntersection);
    BOOST_TEST(A3.isNull());

  } // AreaTest()

  void VolumeTest()
  {
    using Point_t = lar::util::simple_geo::Point3D<float>;
    using Volume_t = lar::util::simple_geo::Volume<Point_t>;

    // default constructor
    Volume_t A1;
    BOOST_TEST(A1.Min() == Volume_t::Point_t{});
    BOOST_TEST(A1.Max() == Volume_t::Point_t{});

    // point constructor
    Volume_t A2({1., 6., 4.}, {4., 2., 8.});
    BOOST_TEST(A2.Min() == Volume_t::Point_t(1., 2., 4.));
    BOOST_TEST(A2.Max() == Volume_t::Point_t(4., 6., 8.));

    Volume_t A3({1., 2., 3.}, {2., 2., 6.});
    Volume_t A4({2., 3., 3.}, {2., 2., 3.});

    // queries
    BOOST_TEST(A2.Center() == Volume_t::Point_t(2.5, 4.0, 6.0));

    BOOST_TEST(A1.DeltaX() == 0.0);
    BOOST_TEST(A2.DeltaX() == 3.0);
    BOOST_TEST(A3.DeltaX() == 1.0);
    BOOST_TEST(A4.DeltaX() == 0.0);

    BOOST_TEST(!A1.isNullX());
    BOOST_TEST(!A2.isNullX());
    BOOST_TEST(!A3.isNullX());
    BOOST_TEST(!A4.isNullX());

    BOOST_TEST(A1.isEmptyX());
    BOOST_TEST(!A2.isEmptyX());
    BOOST_TEST(!A3.isEmptyX());
    BOOST_TEST(A4.isEmptyX());

    BOOST_TEST(A1.nonEmptyDims() == 0U);
    BOOST_TEST(A2.nonEmptyDims() == 3U);
    BOOST_TEST(A3.nonEmptyDims() == 2U);
    BOOST_TEST(A4.nonEmptyDims() == 1U);

    BOOST_TEST(!A1.isNull());
    BOOST_TEST(!A2.isNull());
    BOOST_TEST(!A3.isNull());
    BOOST_TEST(!A4.isNull());

    BOOST_TEST(A1.isEmpty());
    BOOST_TEST(!A2.isEmpty());
    BOOST_TEST(!A3.isEmpty());
    BOOST_TEST(!A4.isEmpty());

    BOOST_TEST(!A1.isLine());
    BOOST_TEST(!A2.isLine());
    BOOST_TEST(!A3.isLine());
    BOOST_TEST(A4.isLine());

    BOOST_TEST(A1.thinnestSize() == 0.0);
    BOOST_TEST(A2.thinnestSize() == 3.0);
    BOOST_TEST(A3.thinnestSize() == 0.0);
    BOOST_TEST(A4.thinnestSize() == 0.0);

    BOOST_TEST(((A1.thinnestSide() == 0) || (A1.thinnestSide() == 1) || (A1.thinnestSide() == 2)));
    BOOST_TEST(A2.thinnestSide() == 0U);
    BOOST_TEST(A3.thinnestSide() == 1U);
    BOOST_TEST(((A4.thinnestSide() == 0) || (A4.thinnestSide() == 2)));

    BOOST_TEST(A1.DeltaY() == 0.0);
    BOOST_TEST(A2.DeltaY() == 4.0);
    BOOST_TEST(A3.DeltaY() == 0.0);
    BOOST_TEST(A4.DeltaY() == 1.0);

    BOOST_TEST(!A1.isNullY());
    BOOST_TEST(!A2.isNullY());
    BOOST_TEST(!A3.isNullY());
    BOOST_TEST(!A4.isNullY());

    BOOST_TEST(A1.isEmptyY());
    BOOST_TEST(!A2.isEmptyY());
    BOOST_TEST(A3.isEmptyY());
    BOOST_TEST(!A4.isEmptyY());

    BOOST_TEST(!A1.isPlane());
    BOOST_TEST(!A2.isPlane());
    BOOST_TEST(A3.isPlane());
    BOOST_TEST(!A4.isPlane());

    BOOST_TEST(A1.DeltaZ() == 0.0);
    BOOST_TEST(A2.DeltaZ() == 4.0);
    BOOST_TEST(A3.DeltaZ() == 3.0);
    BOOST_TEST(A4.DeltaZ() == 0.0);

    BOOST_TEST(!A1.isNullZ());
    BOOST_TEST(!A2.isNullZ());
    BOOST_TEST(!A3.isNullZ());
    BOOST_TEST(!A4.isNullZ());

    BOOST_TEST(A1.isEmptyZ());
    BOOST_TEST(!A2.isEmptyZ());
    BOOST_TEST(!A3.isEmptyZ());
    BOOST_TEST(A4.isEmptyZ());

    BOOST_TEST(!A1.isVolume());
    BOOST_TEST(A2.isVolume());
    BOOST_TEST(!A3.isVolume());
    BOOST_TEST(!A4.isVolume());

    // modifications
    BOOST_TEST(A3 == Volume_t({1.0, 2.0, 3.0}, {2.0, 2.0, 6.0}));
    BOOST_TEST(A2 == Volume_t({1.0, 2.0, 4.0}, {4.0, 6.0, 8.0}));
    BOOST_TEST(A3 != A2);
    A3.Intersect(A2);
    BOOST_TEST(A3 == Volume_t({1.0, 2.0, 4.0}, {2.0, 2.0, 6.0}));

    A3.IncludePoint({1.5, 6.0, 4.0});
    BOOST_TEST(A3 == Volume_t({1.0, 2.0, 4.0}, {2.0, 6.0, 6.0}));
    A3.IncludePoint({4.0, 5.0, 8.0});
    BOOST_TEST(A3 == Volume_t({1.0, 2.0, 4.0}, {4.0, 6.0, 8.0}));
    BOOST_TEST(A3 == A2);

    A3.Intersect(A2);
    BOOST_TEST(A3 == Volume_t({1.0, 2.0, 4.0}, {4.0, 6.0, 8.0}));

    A3.IncludePoint({0.0, 4.0, 4.0});
    BOOST_TEST(A3 == Volume_t({0.0, 2.0, 4.0}, {4.0, 6.0, 8.0}));

    A3.Include(Volume_t({1.0, -1.0, 3.0}, {2.0, 7.0, 7.0}));
    BOOST_TEST(A3 == Volume_t({0.0, -1.0, 3.0}, {4.0, 7.0, 8.0}));

    A3.Include(Volume_t({-1.0, 0.0, 2.0}, {3.0, 1.0, 6.0}));
    BOOST_TEST(A3 == Volume_t({-1.0, -1.0, 2.0}, {4.0, 7.0, 8.0}));

    // intersection result is empty
    BOOST_CHECK_THROW(A3.Intersect(Volume_t({8.0, 9.0, 2.0}, {8.0, 9.0, 3.0})),
                      Volume_t::NullIntersection);
    BOOST_TEST(A3.isNull());
  } // VolumeTest()
} // anon. namespace

// -----------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(SimpleGeo_test)

BOOST_AUTO_TEST_CASE(PointTestCase)
{
  Point2DTest();
  Point3DTest();
}

BOOST_AUTO_TEST_CASE(AreaTestCase)
{
  AreaTest();
  VolumeTest();
}

BOOST_AUTO_TEST_SUITE_END()
