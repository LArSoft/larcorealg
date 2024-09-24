/**
 * @file   geo_vectors_utils_test.cc
 * @brief  Test of geo_vectors_utils.h utilities.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 23, 2017
 */

// Boost libraries
#define BOOST_TEST_MODULE (geo_vectors_test)
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/MetaUtils.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h"

// ROOT libraries
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"

// C/C++ standard library
#include <array>
#include <cmath> // std::nan
#include <iostream>
#include <numeric> // std::iota()
#include <sstream>
#include <stdexcept> // std::runtime_error
#include <string>
#include <type_traits> // std::is_same, std::decay_t
#include <vector>

using boost::test_tools::tolerance;

//------------------------------------------------------------------------------
template <typename PointA, typename PointB>
void CheckPoint(PointA const& test, PointB const& ref, std::string tag = "")
{
  auto const tol = 0.001 % tolerance();

  if (!tag.empty()) BOOST_TEST_CHECKPOINT(tag);

  for (auto ic : geo::vect::indices(test)) {
    BOOST_TEST_CHECKPOINT("  coordinate #" << ic);
    BOOST_TEST(geo::vect::coord(test, ic)() == geo::vect::coord(ref, ic)(), tol);
  }

} // CheckPoint()

//------------------------------------------------------------------------------
void test_MiddlePointAccumulator_defaultConstructor()
{
  geo::Point_t expected(2., 4., 6.);

  std::vector<geo::Point_t> points{
    geo::Point_t(1., 2., 3.), geo::Point_t(2., 4., 6.), geo::Point_t(3., 6., 9.)};
  TVector3 another(expected.X(), expected.Y(), expected.Z());

  //
  // default construction, then bulk addition
  //
  geo::vect::MiddlePointAccumulator acc;
  BOOST_TEST(acc.empty());
  BOOST_TEST(acc.weight() == 0.0, 0.001 % tolerance());
  // add a single point
  acc.add(another);
  BOOST_TEST(!acc.empty());
  BOOST_TEST(acc.weight() == 1.0, 0.001 % tolerance());
  CheckPoint(acc.middlePoint(), expected, "Single add");
  // add many points
  acc.add(points.begin(), points.end());
  BOOST_TEST(!acc.empty());
  BOOST_TEST(acc.weight() == 1.0 + points.size(), 0.001 % tolerance());
  CheckPoint(acc.middlePoint(), expected, "Single add plus sequence");

  //
  // clear test
  //
  acc.clear();
  BOOST_TEST(acc.empty());
  acc.add(geo::Point_t{expected.X() + 1.0, expected.Z(), expected.Y()});
  CheckPoint(
    acc.middlePoint(), geo::Point_t{expected.X() + 1.0, expected.Z(), expected.Y()}, "clear test");

  //
  // start over (same accumulator)
  //
  acc.clear();
  // add many points
  acc.add(points.begin(), points.end());
  BOOST_TEST(!acc.empty());
  CheckPoint(acc.middlePoint(), expected, "Sequence add");
  // add another one
  acc.add(another);
  BOOST_TEST(!acc.empty());
  CheckPoint(acc.middlePoint(), expected, "Sequence add plus single point");

} // test_MiddlePointAccumulator_defaultConstructor()

void test_MiddlePointAccumulator_sequenceConstructor()
{

  geo::Point_t expected(2., 4., 6.);

  std::vector<geo::Point_t> points{
    geo::Point_t(1., 2., 3.), geo::Point_t(2., 4., 6.), geo::Point_t(3., 6., 9.)};
  TVector3 another(expected.X(), expected.Y(), expected.Z());

  //
  // sequence constructor
  //
  geo::vect::MiddlePointAccumulator acc(points.begin(), points.end());
  BOOST_TEST(!acc.empty());
  CheckPoint(acc.middlePoint(), expected, "Sequence construction");
  // add another one
  acc.add(another);
  BOOST_TEST(!acc.empty());
  CheckPoint(acc.middlePoint(), expected, "Sequence construction plus single");

} // test_MiddlePointAccumulator_sequenceConstructor()

//------------------------------------------------------------------------------
template <typename Point>
void test_MiddlePointAccumulator_generic()
{

  constexpr unsigned int Dim = geo::vect::dimension<Point>();
  using Scalar_t = geo::vect::coordinate_t<Point>;

  // prepare the input from larger dimension input data
  constexpr unsigned int MaxDim = 4U;

  static_assert(Dim < MaxDim, "This test supports only up to dimension 4");
  using GenType = std::array<Scalar_t, MaxDim>;
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  GenType const genExpected{{2., 4., 6., 8.}};

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::vector<GenType> genPoints{
    GenType{{1., 2., 3., 4.}}, GenType{{2., 4., 6., 8.}}, GenType{{3., 6., 9., 12.}}};

  std::vector<Point> points;
  for (auto const& genPoint : genPoints)
    points.push_back(geo::vect::makeFromCoords<Point>(genPoint.data()));
  auto const expected = geo::vect::makeFromCoords<Point>(genExpected.data());
  auto const another = expected;

  //
  // default construction, then bulk addition
  //
  geo::vect::MiddlePointAccumulator acc;
  BOOST_TEST(acc.empty());
  BOOST_TEST(acc.weight() == 0.0, 0.001 % tolerance());
  // add a single point
  acc.add(another);
  BOOST_TEST(!acc.empty());
  BOOST_TEST(acc.weight() == 1.0, 0.001 % tolerance());
  CheckPoint(acc.middlePoint(), expected, "Single add");
  // add many points
  acc.add(points.begin(), points.end());
  BOOST_TEST(!acc.empty());
  BOOST_TEST(acc.weight() == 1.0 + points.size(), 0.001 % tolerance());
  CheckPoint(acc.middlePoint(), expected, "Single add plus sequence");

  //
  // clear test
  //
  acc.clear();
  BOOST_TEST(acc.empty());
  acc.add(geo::Point_t{expected.X() + 1.0, expected.Z(), expected.Y()});
  CheckPoint(
    acc.middlePoint(), geo::Point_t{expected.X() + 1.0, expected.Z(), expected.Y()}, "clear test");

  //
  // start over (same accumulator)
  //
  acc.clear();
  // add many points
  acc.add(points.begin(), points.end());
  BOOST_TEST(!acc.empty());
  CheckPoint(acc.middlePoint(), expected, "Sequence add");
  // add another one
  acc.add(another);
  BOOST_TEST(!acc.empty());
  CheckPoint(acc.middlePoint(), expected, "Sequence add plus single point");

} // test_MiddlePointAccumulator_generic()

void test_MiddlePointAccumulator_documentation_class()
{

  geo::Point_t expected{0.0, 1.0, 0.0};
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::array<geo::Point_t, 4> const points = {
   *   geo::Point_t{ 0.0,  1.0,  2.0 },
   *   geo::Point_t{ 0.0, -1.0,  2.0 },
   *   geo::Point_t{ 0.0,  1.0, -2.0 },
   *   geo::Point_t{ 0.0, -1.0, -2.0 }
   * };
   *
   * geo::vect::MiddlePointAccumulator pointsAboveGround;
   * for (auto const& point: points)
   *   if (point.Y() > 0.0) pointsAboveGround.add(point);
   *
   * if (pointsAboveGround.empty())
   *   throw std::runtime_error("No point above ground!");
   *
   * auto middleAboveGround = pointsAboveGround.middlePoint();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<geo::Point_t, 4> const points = {{geo::Point_t{0.0, 1.0, 2.0},
                                               geo::Point_t{0.0, -1.0, 2.0},
                                               geo::Point_t{0.0, 1.0, -2.0},
                                               geo::Point_t{0.0, -1.0, -2.0}}};

  geo::vect::MiddlePointAccumulator pointsAboveGround;
  for (auto const& point : points)
    if (point.Y() > 0.0) pointsAboveGround.add(point);

  if (pointsAboveGround.empty()) throw std::runtime_error("No point above ground!");

  auto middleAboveGround = pointsAboveGround.middlePoint();

  static_assert(std::is_same<decltype(middleAboveGround), geo::Point_t>::value,
                "unexpected return type for geo::vect::MiddlePointAccumulator::middlePoint()");
  CheckPoint(middleAboveGround, expected, "MiddlePointAccumulator::middlePoint()");

} // test_MiddlePointAccumulator_documentation_middlePointAs()

void test_MiddlePointAccumulator_documentation_middlePointAs()
{

  //
  // middlePointAs()
  //
  geo::vect::MiddlePointAccumulator accumulator;
  accumulator.add(geo::Point_t());

  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto mp = accumulator.middlePointAs<TVector3>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  auto mp = accumulator.middlePointAs<TVector3>();

  CheckPoint(mp, geo::Point_t(), "MiddlePointAccumulator::middlePointAs()");

} // test_MiddlePointAccumulator_documentation_middlePointAs()

void test_MiddlePointAccumulator_documentation()
{

  test_MiddlePointAccumulator_documentation_class();
  test_MiddlePointAccumulator_documentation_middlePointAs();

} // test_MiddlePointAccumulator_documentation()

//------------------------------------------------------------------------------
void test_middlePoint()
{

  geo::Point_t expected(2., 4., 6.);

  std::vector<geo::Point_t> points{
    geo::Point_t(1., 2., 3.), geo::Point_t(2., 4., 6.), geo::Point_t(3., 6., 9.)};

  //
  // sequence
  //
  CheckPoint(geo::vect::middlePoint(points.begin(), points.end()), expected, "iterators");

  //
  // points (initializer list)
  //
  CheckPoint(
    geo::vect::middlePoint({points[0], points[1], points[2]}), expected, "initializer list");

  //
  // middlePointAs() (sequence)
  //
  auto const mp3 = geo::vect::middlePointAs<TVector3>(points.begin(), points.end());
  static_assert(std::is_same<std::decay_t<decltype(mp3)>, TVector3>::value,
                "geo::vect::middlePointAs<TVector3> does not return a TVector3!");
  CheckPoint(mp3, expected, "geo::vect::middlePointAs(sequence)");

} // test_middlePoint()

void test_middlePointAs_documentation()
{

  /*
   * std::vector<geo::Point_t> points {
   *   geo::Point_t(1., 2., 3.),
   *   geo::Point_t(2., 4., 6.),
   *   geo::Point_t(3., 6., 9.)
   *   };
   * auto mp = geo::vect::middlePointAs<geo::Vector_t>(points.begin(), points.end());
   */
  std::vector<geo::Point_t> points{
    geo::Point_t(1., 2., 3.), geo::Point_t(2., 4., 6.), geo::Point_t(3., 6., 9.)};
  auto mp = geo::vect::middlePointAs<geo::Vector_t>(points.begin(), points.end());

  static_assert(std::is_same<std::decay_t<decltype(mp)>, geo::Vector_t>::value,
                "geo::vect::middlePointAs<geo::Vector_t> result is not geo::Vector_t");
  CheckPoint(mp, geo::Vector_t(2., 4., 6.));

} // test_middlePointAs_documentation()

void test_middlePoint_iterators_documentation()
{

  /*
   * std::vector<geo::Point_t> points {
   *   geo::Point_t(1., 2., 3.),
   *   geo::Point_t(2., 4., 6.),
   *   geo::Point_t(3., 6., 9.)
   *   };
   *
   * auto mp = geo::vect::middlePoint(points.begin(), points.end());
   */

  std::vector<geo::Point_t> points{
    geo::Point_t(1., 2., 3.), geo::Point_t(2., 4., 6.), geo::Point_t(3., 6., 9.)};

  auto mp = geo::vect::middlePoint(points.begin(), points.end());

  static_assert(std::is_same<std::decay_t<decltype(mp)>, geo::Point_t>::value,
                "geo::vect::middlePoint() result is not geo::Point_t");
  CheckPoint(mp, geo::Point_t(2., 4., 6.));

} // test_middlePoint_iterators_documentation()

void test_middlePoint_initlist_documentation()
{

  /*
   * auto mp = geo::vect::middlePoint
   *   ({ geo::Point_t(1., 2., 3.), geo::Point_t(3., 6., 9.) });
   *
   */

  auto mp = geo::vect::middlePoint({geo::Point_t(1., 2., 3.), geo::Point_t(3., 6., 9.)});

  static_assert(std::is_same<std::decay_t<decltype(mp)>, geo::Point_t>::value,
                "geo::vect::middlePoint() result is not geo::Point_t");
  CheckPoint(mp, geo::Point_t(2., 4., 6.));

} // test_middlePoint_initlist_documentation()

//------------------------------------------------------------------------------
template <typename Vector>
struct test_vectorAccess {

  test_vectorAccess()
  {
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    std::array<double, 3U> const coords = {{1.0, 5.0, 9.0}};
    Vector v{coords[0], coords[1], coords[2]};

    unsigned int iCoord = 0;
    for (auto coordMan : geo::vect::coordManagers<Vector>()) {
      auto const expected = coords[iCoord++];
      auto mc = geo::vect::bindCoord(v, coordMan);
      BOOST_TEST(mc() == expected);
      BOOST_TEST(mc == expected);
    } // for
    BOOST_TEST(iCoord == 3U);

    auto x = geo::vect::Xcoord(v);
    auto c0 = geo::vect::coord(v, 0U);
    auto mx = geo::vect::bindCoord(v, geo::vect::XcoordManager<Vector>);
    auto mc0 = geo::vect::bindCoord(v, geo::vect::coordManager<Vector>(0U));
    BOOST_TEST(x() == 1.0);
    BOOST_TEST(x == 1.0);
    BOOST_TEST(x() == v.X());
    BOOST_TEST(c0() == v.X());
    BOOST_TEST(mx() == 1.0);
    BOOST_TEST(mx == 1.0);
    BOOST_TEST(mx() == v.X());
    BOOST_TEST(mc0() == v.X());

    x = 2.0;
    BOOST_TEST(x() == 2.0);
    BOOST_TEST(x == 2.0);
    BOOST_TEST(x() == v.X());
    BOOST_TEST(c0() == v.X());
    BOOST_TEST(mx() == 2.0);
    BOOST_TEST(mx == 2.0);
    BOOST_TEST(mx() == v.X());
    BOOST_TEST(mc0() == v.X());

    x += 2.0;
    BOOST_TEST(x() == 4.0);
    BOOST_TEST(x == 4.0);
    BOOST_TEST(x() == v.X());
    BOOST_TEST(c0() == v.X());
    BOOST_TEST(mx() == 4.0);
    BOOST_TEST(mx == 4.0);
    BOOST_TEST(mx() == v.X());
    BOOST_TEST(mc0() == v.X());

    x -= 2.0;
    BOOST_TEST(x() == 2.0);
    BOOST_TEST(x == 2.0);
    BOOST_TEST(x() == v.X());
    BOOST_TEST(c0() == v.X());
    BOOST_TEST(mx() == 2.0);
    BOOST_TEST(mx == 2.0);
    BOOST_TEST(mx() == v.X());
    BOOST_TEST(mc0() == v.X());

    x *= 4.0;
    BOOST_TEST(x() == 8.0);
    BOOST_TEST(x == 8.0);
    BOOST_TEST(x() == v.X());
    BOOST_TEST(c0() == v.X());
    BOOST_TEST(mx() == 8.0);
    BOOST_TEST(mx == 8.0);
    BOOST_TEST(mx() == v.X());
    BOOST_TEST(mc0() == v.X());

    x /= 4.0;
    BOOST_TEST(x() == 2.0);
    BOOST_TEST(x == 2.0);
    BOOST_TEST(x() == v.X());
    BOOST_TEST(c0() == v.X());
    BOOST_TEST(mx() == 2.0);
    BOOST_TEST(mx == 2.0);
    BOOST_TEST(mx() == v.X());
    BOOST_TEST(mc0() == v.X());

    auto y = geo::vect::Ycoord(v);
    auto c1 = geo::vect::coord(v, 1U);
    auto my = geo::vect::bindCoord(v, geo::vect::YcoordManager<Vector>);
    auto mc1 = geo::vect::bindCoord(v, geo::vect::coordManager<Vector>(1U));
    BOOST_TEST(y() == 5.0);
    BOOST_TEST(y == 5.0);
    BOOST_TEST(y() == v.Y());
    BOOST_TEST(c1() == v.Y());
    BOOST_TEST(my() == 5.0);
    BOOST_TEST(my == 5.0);
    BOOST_TEST(my() == v.Y());
    BOOST_TEST(mc1() == v.Y());

    y = 2.0;
    BOOST_TEST(y() == 2.0);
    BOOST_TEST(y == 2.0);
    BOOST_TEST(y() == v.Y());
    BOOST_TEST(c1() == v.Y());
    BOOST_TEST(my() == 2.0);
    BOOST_TEST(my == 2.0);
    BOOST_TEST(my() == v.Y());
    BOOST_TEST(mc1() == v.Y());

    y += 2.0;
    BOOST_TEST(y() == 4.0);
    BOOST_TEST(y == 4.0);
    BOOST_TEST(y() == v.Y());
    BOOST_TEST(c1() == v.Y());
    BOOST_TEST(my() == 4.0);
    BOOST_TEST(my == 4.0);
    BOOST_TEST(my() == v.Y());
    BOOST_TEST(mc1() == v.Y());

    y -= 2.0;
    BOOST_TEST(y() == 2.0);
    BOOST_TEST(y == 2.0);
    BOOST_TEST(y() == v.Y());
    BOOST_TEST(c1() == v.Y());
    BOOST_TEST(my() == 2.0);
    BOOST_TEST(my == 2.0);
    BOOST_TEST(my() == v.Y());
    BOOST_TEST(mc1() == v.Y());

    y *= 4.0;
    BOOST_TEST(y() == 8.0);
    BOOST_TEST(y == 8.0);
    BOOST_TEST(y() == v.Y());
    BOOST_TEST(c1() == v.Y());
    BOOST_TEST(my() == 8.0);
    BOOST_TEST(my == 8.0);
    BOOST_TEST(my() == v.Y());
    BOOST_TEST(mc1() == v.Y());

    y /= 4.0;
    BOOST_TEST(y() == 2.0);
    BOOST_TEST(y == 2.0);
    BOOST_TEST(y() == v.Y());
    BOOST_TEST(c1() == v.Y());
    BOOST_TEST(my() == 2.0);
    BOOST_TEST(my == 2.0);
    BOOST_TEST(my() == v.Y());
    BOOST_TEST(mc1() == v.Y());

    auto z = geo::vect::Zcoord(v);
    auto c2 = geo::vect::coord(v, 2U);
    auto mz = geo::vect::bindCoord(v, geo::vect::ZcoordManager<Vector>);
    auto mc2 = geo::vect::bindCoord(v, geo::vect::coordManager<Vector>(2U));
    BOOST_TEST(z() == 9.0);
    BOOST_TEST(z == 9.0);
    BOOST_TEST(z() == v.Z());
    BOOST_TEST(c2() == v.Z());
    BOOST_TEST(mz() == 9.0);
    BOOST_TEST(mz == 9.0);
    BOOST_TEST(mz() == v.Z());
    BOOST_TEST(mc2() == v.Z());

    z = 2.0;
    BOOST_TEST(z() == 2.0);
    BOOST_TEST(z == 2.0);
    BOOST_TEST(z() == v.Z());
    BOOST_TEST(c2() == v.Z());
    BOOST_TEST(mz() == 2.0);
    BOOST_TEST(mz == 2.0);
    BOOST_TEST(mz() == v.Z());
    BOOST_TEST(mc2() == v.Z());

    z += 2.0;
    BOOST_TEST(z() == 4.0);
    BOOST_TEST(z == 4.0);
    BOOST_TEST(z() == v.Z());
    BOOST_TEST(c2() == v.Z());
    BOOST_TEST(mz() == 4.0);
    BOOST_TEST(mz == 4.0);
    BOOST_TEST(mz() == v.Z());
    BOOST_TEST(mc2() == v.Z());

    z -= 2.0;
    BOOST_TEST(z() == 2.0);
    BOOST_TEST(z == 2.0);
    BOOST_TEST(z() == v.Z());
    BOOST_TEST(c2() == v.Z());
    BOOST_TEST(mz() == 2.0);
    BOOST_TEST(mz == 2.0);
    BOOST_TEST(mz() == v.Z());
    BOOST_TEST(mc2() == v.Z());

    z *= 4.0;
    BOOST_TEST(z() == 8.0);
    BOOST_TEST(z == 8.0);
    BOOST_TEST(z() == v.Z());
    BOOST_TEST(c2() == v.Z());
    BOOST_TEST(mz() == 8.0);
    BOOST_TEST(mz == 8.0);
    BOOST_TEST(mz() == v.Z());
    BOOST_TEST(mc2() == v.Z());

    z /= 4.0;
    BOOST_TEST(z() == 2.0);
    BOOST_TEST(z == 2.0);
    BOOST_TEST(z() == v.Z());
    BOOST_TEST(c2() == v.Z());
    BOOST_TEST(mz() == 2.0);
    BOOST_TEST(mz == 2.0);
    BOOST_TEST(mz() == v.Z());
    BOOST_TEST(mc2() == v.Z());

  } // test_vectorAccess()

}; // struct test_vectorAccess

template <typename Vector>
struct test_vectorAccess<Vector const> {

  test_vectorAccess()
  {

    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    std::array<double, 3U> const coords = {{1.0, 5.0, 9.0}};
    Vector const v{coords[0], coords[1], coords[2]};

    unsigned int iCoord = 0;
    for (auto coordMan : geo::vect::coordReaders<Vector>()) {
      auto const expected = coords[iCoord++];
      auto mc = geo::vect::bindCoord(v, coordMan);
      BOOST_TEST(mc() == expected);
      BOOST_TEST(mc == expected);
    } // for
    BOOST_TEST(iCoord == 3U);

    auto x = geo::vect::Xcoord(v);
    auto c0 = geo::vect::coord(v, 0U);
    auto mx = geo::vect::bindCoord(v, geo::vect::XcoordManager<Vector>);
    auto mc0 = geo::vect::bindCoord(v, geo::vect::coordManager<Vector>(0U));
    BOOST_TEST(x() == 1.0);
    BOOST_TEST(x == 1.0);
    BOOST_TEST(x() == v.X());
    BOOST_TEST(c0() == v.X());
    BOOST_TEST(mx() == 1.0);
    BOOST_TEST(mx == 1.0);
    BOOST_TEST(mx() == v.X());
    BOOST_TEST(mc0() == v.X());

    auto y = geo::vect::Ycoord(v);
    auto c1 = geo::vect::coord(v, 1U);
    auto my = geo::vect::bindCoord(v, geo::vect::YcoordManager<Vector>);
    auto mc1 = geo::vect::bindCoord(v, geo::vect::coordManager<Vector>(1U));
    BOOST_TEST(y() == 5.0);
    BOOST_TEST(y == 5.0);
    BOOST_TEST(y() == v.Y());
    BOOST_TEST(c1() == v.Y());
    BOOST_TEST(my() == 5.0);
    BOOST_TEST(my == 5.0);
    BOOST_TEST(my() == v.Y());
    BOOST_TEST(mc1() == v.Y());

    auto z = geo::vect::Zcoord(v);
    auto c2 = geo::vect::coord(v, 2U);
    auto mz = geo::vect::bindCoord(v, geo::vect::ZcoordManager<Vector>);
    auto mc2 = geo::vect::bindCoord(v, geo::vect::coordManager<Vector>(2U));
    BOOST_TEST(z() == 9.0);
    BOOST_TEST(z == 9.0);
    BOOST_TEST(z() == v.Z());
    BOOST_TEST(c2() == v.Z());
    BOOST_TEST(mz() == 9.0);
    BOOST_TEST(mz == 9.0);
    BOOST_TEST(mz() == v.Z());
    BOOST_TEST(mc2() == v.Z());

  } // test_vectorAccess()
};  // struct test_vectorAccess<const>

//------------------------------------------------------------------------------
template <typename Vector, unsigned int Dim = geo::vect::dimension<Vector>()>
struct IsfiniteTester;

template <typename Vector>
struct IsfiniteTester<Vector, 4U> {
  IsfiniteTester()
  {
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 2.0, 3.0, 4.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{0.0, 2.0, 3.0, 4.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 0.0, 3.0, 4.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 2.0, 0.0, 4.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 2.0, 3.0, 0.0}));

    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), 2.0, 3.0, 4.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, std::nan(""), 3.0, 4.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, 2.0, std::nan(""), 4.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, 2.0, 3.0, std::nan("")}));

    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), std::nan(""), 3.0, 4.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), 2.0, std::nan(""), 4.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), 2.0, 3.0, std::nan("")}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, std::nan(""), std::nan(""), 4.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, std::nan(""), 3.0, std::nan("")}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, 2.0, std::nan(""), std::nan("")}));

    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, std::nan(""), std::nan(""), std::nan("")}));
    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), 2.0, std::nan(""), std::nan("")}));
    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), std::nan(""), 3.0, std::nan("")}));
    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), std::nan(""), std::nan(""), 4.0}));

    BOOST_TEST(
      !geo::vect::isfinite(Vector{std::nan(""), std::nan(""), std::nan(""), std::nan("")}));
  }
};

template <typename Vector>
struct IsfiniteTester<Vector, 3U> {
  IsfiniteTester()
  {
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 2.0, 3.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{0.0, 2.0, 3.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 0.0, 3.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 2.0, 0.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 2.0, 3.0}));

    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), 2.0, 3.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, std::nan(""), 3.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, 2.0, std::nan("")}));

    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, std::nan(""), std::nan("")}));
    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), 2.0, std::nan("")}));
    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), std::nan(""), 3.0}));

    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), std::nan(""), std::nan("")}));
  }
};

template <typename Vector>
struct IsfiniteTester<Vector, 2U> {
  IsfiniteTester()
  {
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 2.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{0.0, 2.0}));
    BOOST_TEST(geo::vect::isfinite(Vector{1.0, 0.0}));

    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), 2.0}));
    BOOST_TEST(!geo::vect::isfinite(Vector{1.0, std::nan("")}));

    BOOST_TEST(!geo::vect::isfinite(Vector{std::nan(""), std::nan("")}));
  }
};

template <typename Vector>
void test_vectorProcessing()
{

  (void)IsfiniteTester<Vector>();

} // test_vectorProcessing()

//------------------------------------------------------------------------------
template <typename Source, typename Dest>
void test_vector2Dconvert()
{

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<double, 4U> srcData{{1.0, 5.0, 9.0, 16.0}};

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  // Source const src{ srcData[0], srcData[1] };
  Source srcForClangBug;
  if constexpr (util::is_STLarray_v<Source>)
    srcForClangBug = Source{{srcData[0], srcData[1]}};
  else
    srcForClangBug = Source{srcData[0], srcData[1]};
  Source const src{srcForClangBug};

  auto dest = geo::vect::convertTo<Dest>(src);

  static_assert(std::is_same<decltype(dest), Dest>(), "Unexpected return type!");

  BOOST_TEST(geo::vect::Xcoord(dest) == srcData[0]);
  BOOST_TEST(geo::vect::Ycoord(dest) == srcData[1]);

} // test_vector2Dconvert()

//------------------------------------------------------------------------------
template <typename Source, typename Dest>
void test_vector3Dconvert()
{

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<double, 4U> srcData{{1.0, 5.0, 9.0, 16.0}};

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  // Source const src{ srcData[0], srcData[1], srcData[2] };
  Source srcForClangBug;
  if constexpr (util::is_STLarray_v<Source>)
    srcForClangBug = Source{{srcData[0], srcData[1], srcData[2]}};
  else
    srcForClangBug = Source{srcData[0], srcData[1], srcData[2]};
  Source const src{srcForClangBug};

  auto dest = geo::vect::convertTo<Dest>(src);

  static_assert(std::is_same<decltype(dest), Dest>(), "Unexpected return type!");

  BOOST_TEST(geo::vect::Xcoord(dest) == srcData[0]);
  BOOST_TEST(geo::vect::Ycoord(dest) == srcData[1]);
  BOOST_TEST(geo::vect::Zcoord(dest) == srcData[2]);

} // test_vector3Dconvert()

//------------------------------------------------------------------------------
template <typename Source, typename Dest>
void test_vector4Dconvert()
{

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<double, 4U> srcData{{1.0, 5.0, 9.0, 16.0}};

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  // Source const src{ srcData[0], srcData[1], srcData[2], srcData[3] };
  Source srcForClangBug;
  if constexpr (util::is_STLarray_v<Source>)
    srcForClangBug = Source{{srcData[0], srcData[1], srcData[2], srcData[3]}};
  else
    srcForClangBug = Source{srcData[0], srcData[1], srcData[2], srcData[3]};
  Source const src{srcForClangBug};

  auto dest = geo::vect::convertTo<Dest>(src);

  static_assert(std::is_same<decltype(dest), Dest>(), "Unexpected return type!");

  BOOST_TEST(geo::vect::Xcoord(dest) == srcData[0]);
  BOOST_TEST(geo::vect::Ycoord(dest) == srcData[1]);
  BOOST_TEST(geo::vect::Zcoord(dest) == srcData[2]);
  BOOST_TEST(geo::vect::Tcoord(dest) == srcData[3]);

} // test_vector4Dconvert()

//------------------------------------------------------------------------------
void test_makeFromCoords_documentation()
{

  /*
   * constexpr std::array<float, 5U> data { 2.0, 5.0, 7.0, 11.0, 15.5 };
   * constexpr auto p = geo::vect::makeFromCoords<geo::Point_t>(data);
   * auto v = geo::vect::makeFromCoords<geo::Vector_t>(data.data() + 1);
   */
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  constexpr std::array<float, 5U> data{{2.0, 5.0, 7.0, 11.0, 15.5}};
  auto const p = geo::vect::makeFromCoords<geo::Point_t>(data);
  auto const v = geo::vect::makeFromCoords<geo::Vector_t>(data.data() + 1);

  BOOST_TEST(p == (geo::Point_t{2.0, 5.0, 7.0}));
  BOOST_TEST(v == (geo::Vector_t{5.0, 7.0, 11.0}));

} // test_makeFromCoords_documentation()

//------------------------------------------------------------------------------
template <typename Vector>
void test_transform()
{

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<float, 5U> data{{4.0, 5.0, 7.0, 11.0, 15.5}};
  auto const v = geo::vect::makeFromCoords<Vector>(data);

  auto const neg_v = geo::vect::transformCoords(v, [](auto c) { return -c; });
  static_assert(std::is_same<decltype(neg_v), decltype(v)>(), "Unexpected return type");
  BOOST_TEST(neg_v == -v);

} // test_transform()

//------------------------------------------------------------------------------
void test_XcoordManager_documentation()
{
  std::ostringstream expected, out;
  /*
   * // constant vectors get a "reader" (read-only manager):
   * geo::Vector_t v { 1.0, 2.0, 3.0 };
   *
   * auto vx = geo::vect::bindCoord(v, geo::vect::XcoordManager<geo::Vector_t const>);
   * std::cout << v << " has x=" << vx() << std::endl;
   */

  // constant vectors get a "reader" (read-only manager):
  geo::Vector_t v{1.0, 2.0, 3.0};

  auto vx = geo::vect::bindCoord(v, geo::vect::XcoordManager<geo::Vector_t const>);
  std::cout << v << " has x=" << vx() << std::endl;
  out << v << " has x=" << vx();
  expected << v << " has x=" << v.X();

  BOOST_TEST(out.str() == expected.str());

  out.str("");
  expected.str("");
  /*
   * // mutable vectors get a full-featured "manager":
   * geo::Point_t p { 1.0, 2.0, 3.0 };
   * auto px = geo::vect::bindCoord(p, geo::vect::XcoordManager<geo::Point_t>);
   * px *= 5.0;
   * std::cout << p << " has now x=" << px() << std::endl;
   */

  // mutable vectors get a full-featured "manager":
  geo::Point_t p{1.0, 2.0, 3.0};
  auto px = geo::vect::bindCoord(p, geo::vect::XcoordManager<geo::Point_t>);
  px *= 5.0;
  std::cout << p << " has now x=" << px() << std::endl;
  out << p << " has now x=" << px();
  expected << p << " has now x=" << p.X();

  BOOST_TEST(out.str() == expected.str());

} // test_XcoordManager_documentation()

//------------------------------------------------------------------------------
struct VectorTraitsTester {
  // artificial ter vector types
  template <typename C>
  struct Vector0D {
    static constexpr unsigned int Dim = 0U;
    using Scalar = C;
  };
  template <typename C>
  struct Vector1D : public Vector0D<C> {
    static constexpr unsigned int Dim = 1U;
    C X() const;
  };
  template <typename C>
  struct Vector2D : public Vector1D<C> {
    static constexpr unsigned int Dim = 2U;
    C Y() const;
  };
  template <typename C>
  struct Vector3D : public Vector2D<C> {
    static constexpr unsigned int Dim = 3U;
    C Z() const;
  };
  template <typename C>
  struct Vector4D : public Vector3D<C> {
    static constexpr unsigned int Dim = 4U;
    C T() const;
  };
  template <typename C>
  struct Vector5D : public Vector4D<C> {
    static constexpr unsigned int Dim = 5U;
    C U() const;
  };

  // the test on details is here just to facilitate debugging;
  // failure of details tests is acceptable (and failing tests must be removed)
  static_assert(!geo::vect::details::HasX<Vector0D<double>>(), "Unexpected 0D::X()");
  static_assert(geo::vect::details::HasX<Vector1D<double>>(), "Unexpected 1D::X()");
  static_assert(geo::vect::details::HasX<Vector2D<double>>(), "Unexpected 2D::X()");
  static_assert(geo::vect::details::HasX<Vector3D<double>>(), "Unexpected 3D::X()");
  static_assert(geo::vect::details::HasX<Vector4D<double>>(), "Unexpected 4D::X()");
  static_assert(geo::vect::details::HasX<Vector5D<double>>(), "Unexpected 5D::X()");

  static_assert(!geo::vect::details::HasY<Vector0D<double>>(), "Unexpected 0D::Y()");
  static_assert(!geo::vect::details::HasY<Vector1D<double>>(), "Unexpected 1D::Y()");
  static_assert(geo::vect::details::HasY<Vector2D<double>>(), "Unexpected 2D::Y()");
  static_assert(geo::vect::details::HasY<Vector3D<double>>(), "Unexpected 3D::Y()");
  static_assert(geo::vect::details::HasY<Vector4D<double>>(), "Unexpected 4D::Y()");
  static_assert(geo::vect::details::HasY<Vector5D<double>>(), "Unexpected 5D::Y()");

  static_assert(!geo::vect::details::HasZ<Vector0D<double>>(), "Unexpected 0D::Z()");
  static_assert(!geo::vect::details::HasZ<Vector1D<double>>(), "Unexpected 1D::Z()");
  static_assert(!geo::vect::details::HasZ<Vector2D<double>>(), "Unexpected 2D::Z()");
  static_assert(geo::vect::details::HasZ<Vector3D<double>>(), "Unexpected 3D::Z()");
  static_assert(geo::vect::details::HasZ<Vector4D<double>>(), "Unexpected 4D::Z()");
  static_assert(geo::vect::details::HasZ<Vector5D<double>>(), "Unexpected 5D::Z()");

  static_assert(!geo::vect::details::HasT<Vector0D<double>>(), "Unexpected 0D::T()");
  static_assert(!geo::vect::details::HasT<Vector1D<double>>(), "Unexpected 1D::T()");
  static_assert(!geo::vect::details::HasT<Vector2D<double>>(), "Unexpected 2D::T()");
  static_assert(!geo::vect::details::HasT<Vector3D<double>>(), "Unexpected 3D::T()");
  static_assert(geo::vect::details::HasT<Vector4D<double>>(), "Unexpected 4D::T()");
  static_assert(geo::vect::details::HasT<Vector5D<double>>(), "Unexpected 5D::T()");

  static_assert(geo::vect::dimension<Vector0D<double>>() == 0U, "Unexpected 0D dimension");
  static_assert(geo::vect::dimension<Vector1D<double>>() == 1U, "Unexpected 1D dimension");
  static_assert(geo::vect::dimension<Vector2D<double>>() == 2U, "Unexpected 2D dimension");
  static_assert(geo::vect::dimension<Vector3D<double>>() == 3U, "Unexpected 3D dimension");
  static_assert(geo::vect::dimension<Vector4D<double>>() == 4U, "Unexpected 4D dimension");
  static_assert(geo::vect::dimension<Vector5D<double>>() == 4U, "Unexpected 5D dimension");

  //
  // TVector2
  //
  static_assert(geo::vect::details::HasX<TVector2>(), "Unexpected TVector2::X()");
  static_assert(geo::vect::details::HasY<TVector2>(), "Unexpected TVector2::Y()");
  static_assert(!geo::vect::details::HasZ<TVector2>(), "Unexpected TVector2::Z()");
  static_assert(!geo::vect::details::HasT<TVector2>(), "Unexpected TVector2::T()");
  static_assert(geo::vect::dimension<TVector2>() == 2U, "Unexpected TVector2 dimension");
  static_assert(std::is_same<geo::vect::coordinate_t<TVector2>, double>(),
                "Unexpected vector type");

  //
  // TVector3
  //
  static_assert(geo::vect::details::HasX<TVector3>(), "Unexpected TVector3::X()");
  static_assert(geo::vect::details::HasY<TVector3>(), "Unexpected TVector3::Y()");
  static_assert(geo::vect::details::HasZ<TVector3>(), "Unexpected TVector3::Z()");
  static_assert(!geo::vect::details::HasT<TVector3>(), "Unexpected TVector3::T()");
  static_assert(geo::vect::dimension<TVector3>() == 3U, "Unexpected TVector3 dimension");
  static_assert(std::is_same<geo::vect::coordinate_t<TVector3>, double>(),
                "Unexpected vector type");

  //
  // TLorentzVector
  //
  static_assert(geo::vect::details::HasX<TLorentzVector>(), "Unexpected TLorentzVector::X()");
  static_assert(geo::vect::details::HasY<TLorentzVector>(), "Unexpected TLorentzVector::Y()");
  static_assert(geo::vect::details::HasZ<TLorentzVector>(), "Unexpected TLorentzVector::Z()");
  static_assert(geo::vect::details::HasT<TLorentzVector>(), "Unexpected TLorentzVector::T()");
  static_assert(geo::vect::dimension<TLorentzVector>() == 4U,
                "Unexpected TLorentzVector dimension");
  static_assert(std::is_same<geo::vect::coordinate_t<TLorentzVector>, double>(),
                "Unexpected TLorentzVector coordinate type");

  //
  // geo::Vector_t
  //
  static_assert(geo::vect::details::HasX<geo::Vector_t>(), "Unexpected geo::Vector_t::X()");
  static_assert(geo::vect::details::HasY<geo::Vector_t>(), "Unexpected geo::Vector_t::Y()");
  static_assert(geo::vect::details::HasZ<geo::Vector_t>(), "Unexpected geo::Vector_t::Z()");
  static_assert(!geo::vect::details::HasT<geo::Vector_t>(), "Unexpected geo::Vector_t::T()");
  static_assert(geo::vect::dimension<geo::Vector_t>() == 3U, "Unexpected geo::Vector_t dimension");
  static_assert(std::is_same<geo::vect::coordinate_t<geo::Vector_t>, double>(),
                "Unexpected vector type");

  //
  // geo::Point_t
  //
  static_assert(geo::vect::details::HasX<geo::Point_t>(), "Unexpected geo::Point_t::X()");
  static_assert(geo::vect::details::HasY<geo::Point_t>(), "Unexpected geo::Point_t::Y()");
  static_assert(geo::vect::details::HasZ<geo::Point_t>(), "Unexpected geo::Point_t::Z()");
  static_assert(!geo::vect::details::HasT<geo::Point_t>(), "Unexpected geo::Point_t::T()");
  static_assert(geo::vect::dimension<geo::Point_t>() == 3U, "Unexpected geo::Point_t dimension");
  static_assert(std::is_same<geo::vect::coordinate_t<geo::Point_t>, double>(),
                "Unexpected vector type");

}; // struct VectorTraitsTester

//------------------------------------------------------------------------------
template <typename Vector>
void test_CoordConstIterator()
{

  using Vector_t = Vector;
  using Coord_t = geo::vect::coordinate_t<Vector_t>;

  std::array<Coord_t, geo::vect::dimension<Vector_t>()> expected;
  std::iota(expected.begin(), expected.end(), Coord_t(1));
  auto const v = geo::vect::makeFromCoords<Vector_t>(expected);

  unsigned int index = 0;
  for (Coord_t c : v) {
    BOOST_TEST(c == expected[index]);
    ++index;
  } // for

} // test_CoordConstIterator()

//------------------------------------------------------------------------------
template <typename Vector>
void test_fillCoords()
{

  using Vector_t = Vector;
  using Coord_t = geo::vect::coordinate_t<Vector_t>;

  std::array<Coord_t, geo::vect::dimension<Vector_t>()> expected;
  std::iota(expected.begin(), expected.end(), Coord_t(1));
  auto const v = geo::vect::makeFromCoords<Vector_t>(expected);

  Coord_t coords[geo::vect::dimension<Vector_t>()];
  auto const dim = geo::vect::fillCoords(coords, v);

  BOOST_TEST(dim == expected.size());

  for (unsigned int index = 0; index < dim; ++index)
    BOOST_TEST(coords[index] == expected[index]);

} // test_fillCoords()

BOOST_AUTO_TEST_SUITE(geo_vectors_utils_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(MiddlePointAccumulator_test)
{
  test_MiddlePointAccumulator_defaultConstructor();
  test_MiddlePointAccumulator_sequenceConstructor();
  test_MiddlePointAccumulator_documentation();
  test_MiddlePointAccumulator_generic<TVector3>();
  test_MiddlePointAccumulator_generic<geo::Vector_t>();
  test_MiddlePointAccumulator_generic<geo::Point_t>();
} // BOOST_AUTO_TEST_CASE(MiddlePointAccumulator_test)

BOOST_AUTO_TEST_CASE(middlePoint_test)
{
  test_middlePoint();
}

BOOST_AUTO_TEST_CASE(middlePoint_documentation_test)
{
  test_middlePointAs_documentation();
  test_middlePoint_iterators_documentation();
  test_middlePoint_initlist_documentation();
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(vectorAccess_test)
{
  test_vectorAccess<TVector3>();
  test_vectorAccess<geo::Point_t>();
  test_vectorAccess<geo::Vector_t>();
  test_vectorAccess<TVector3 const>();
  test_vectorAccess<geo::Point_t const>();
  test_vectorAccess<geo::Vector_t const>();
} // BOOST_AUTO_TEST_CASE(vectorAccess_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(vectorUtilDocumentation_test)
{
  test_XcoordManager_documentation();
  test_makeFromCoords_documentation();
} // BOOST_AUTO_TEST_CASE(vectorUtilDocumentation_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(vectorProperties_test)
{
  (void)VectorTraitsTester();
} // BOOST_AUTO_TEST_CASE(vectorConversion_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(vectorProcessing_test)
{
  test_vectorProcessing<TVector2>();
  test_vectorProcessing<geo::Point_t>();
  test_vectorProcessing<geo::Vector_t>();
  test_vectorProcessing<TVector3>();
  test_vectorProcessing<TLorentzVector>();
} // BOOST_AUTO_TEST_CASE(vectorProcessing_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(vectorConversion_test)
{
  // BUG until clang bug 21629 is fixed, the amount of work to make the *test*
  //     work on double[N] vectors is not worth; HERE are some things to restore
  //     when the bug is fixed:
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)

  // 2D
  test_vector2Dconvert<TVector2, TVector2>();
  test_vector2Dconvert<std::array<double, 2U>, TVector2>();
  // HERE
  // test_vector2Dconvert<double[2U]            , TVector2              >();

  // 3D
  test_vector3Dconvert<std::array<double, 3U>, TVector3>();

  // HERE
  // test_vector3Dconvert<double[3U]            , TVector3     >();
  test_vector3Dconvert<TVector3, TVector3>();
  test_vector3Dconvert<geo::Point_t, TVector3>();
  test_vector3Dconvert<geo::Vector_t, TVector3>();
  test_vector3Dconvert<std::array<double, 3U>, geo::Point_t>();
  // HERE
  // test_vector3Dconvert<double[3U]            , geo::Point_t >();
  test_vector3Dconvert<TVector3, geo::Point_t>();
  test_vector3Dconvert<geo::Point_t, geo::Point_t>();
  test_vector3Dconvert<geo::Vector_t, geo::Point_t>();
  test_vector3Dconvert<std::array<double, 3U>, geo::Vector_t>();
  // HERE
  // test_vector3Dconvert<double[3U]            , geo::Vector_t>();
  test_vector3Dconvert<TVector3, geo::Vector_t>();
  test_vector3Dconvert<geo::Point_t, geo::Vector_t>();
  test_vector3Dconvert<geo::Vector_t, geo::Vector_t>();

  test_transform<TVector3>();

  // 4D
  test_vector4Dconvert<TLorentzVector, TLorentzVector>();
  test_vector4Dconvert<std::array<double, 4U>, TLorentzVector>();
  // HERE
  // test_vector4Dconvert<double[4U]            , TLorentzVector>();

} // BOOST_AUTO_TEST_CASE(vectorAccess_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(vectorCoordinateIteration_test)
{

  test_CoordConstIterator<TVector2>();
  test_CoordConstIterator<TVector3>();
  test_CoordConstIterator<geo::Point_t>();
  test_CoordConstIterator<geo::Vector_t>();
  test_CoordConstIterator<TLorentzVector>();

  test_fillCoords<TVector2>();
  test_fillCoords<TVector3>();
  test_fillCoords<geo::Point_t>();
  test_fillCoords<geo::Vector_t>();
  test_fillCoords<TLorentzVector>();

} // BOOST_AUTO_TEST_CASE(vectorCoordinateIteration_test)

//------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
