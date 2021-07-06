/**
 * @file   LineClosestPoint_test.cc
 * @brief  Test of `LineClosestPoint.h` utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 29, 2021
 * @see    `larcorealg/Geometry/LineClosestPoint.h`
 */


// Boost libraries
#define BOOST_TEST_MODULE ( LineClosestPoint_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_CLOSE()

// LArSoft libraries
#include "larcorealg/Geometry/LineClosestPoint.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// C++ standard library
#include <utility> // std::pair<>
#include <cmath> // std::sqrt()


// =============================================================================
void LineClosestPointSimple_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  geo::Point_t const p = geo::LineClosestPoint(
    geo::origin() - geo::Zaxis() - 3.0 * geo::Xaxis(), geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), geo::Yaxis()
    );
  
  BOOST_CHECK_SMALL(p.X(),       tol);
  BOOST_CHECK_SMALL(p.Y(),       tol);
  BOOST_CHECK_CLOSE(p.Z(), -1.0, tol);
  
} // LineClosestPointSimple_test()


// -----------------------------------------------------------------------------
void LineClosestPointSimple45_test() {
  
  constexpr double tol = 0.001; // absolute
  
  geo::Point_t const p = geo::LineClosestPoint(
    geo::origin(),                geo::Xaxis(),
    geo::origin() + geo::Yaxis(), (geo::Xaxis() + geo::Zaxis()) / std::sqrt(2.0)
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  
} // LineClosestPointSimple45_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithLocOnLines_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  std::pair<double, double> locOnLines { 0.0, 0.0 };
  geo::Point_t const p = geo::LineClosestPoint(
    geo::origin()                - 3.0 * geo::Xaxis(), geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), geo::Yaxis(),
    &locOnLines
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  BOOST_CHECK_CLOSE(locOnLines.first,  +3.0, tol);
  BOOST_CHECK_CLOSE(locOnLines.second, -2.0, tol);
  
} // LineClosestPointWithLocOnLines_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithScaledDirs_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  std::pair<double, double> locOnLines { 0.0, 0.0 };
  geo::Point_t const p = geo::LineClosestPoint(
    geo::origin()                - 3.0 * geo::Xaxis(),  2.0 * geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), -2.0 * geo::Yaxis(),
    &locOnLines
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  BOOST_CHECK_CLOSE(locOnLines.first,  +3.0 /  2.0, tol);
  BOOST_CHECK_CLOSE(locOnLines.second, -2.0 / -2.0, tol);
  
} // LineClosestPointWithScaledDirs_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithNonHomogeneousDirs_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  std::pair<double, double> locOnLines { 0.0, 0.0 };
  geo::Point_t const p = geo::LineClosestPoint(
    geo::origin()                - 3.0 * geo::Xaxis(),  1.5 * geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), -2.0 * geo::Yaxis(),
    &locOnLines
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  BOOST_CHECK_CLOSE(locOnLines.first,  +3.0 /  1.5, tol);
  BOOST_CHECK_CLOSE(locOnLines.second, -2.0 / -2.0, tol);
  
} // LineClosestPointWithNonHomogeneousDirs_test()


// =============================================================================
void LineClosestPointWithUnitVectorsSimple_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  geo::Point_t const p = geo::LineClosestPointWithUnitVectors(
    geo::origin() - geo::Zaxis() - 3.0 * geo::Xaxis(), geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), geo::Yaxis()
    );
  
  BOOST_CHECK_SMALL(p.X(),       tol);
  BOOST_CHECK_SMALL(p.Y(),       tol);
  BOOST_CHECK_CLOSE(p.Z(), -1.0, tol);
  
} // LineClosestPointWithUnitVectorsSimple_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithUnitVectorsSimple45_test() {
  
  constexpr double tol = 0.001; // absolute
  
  geo::Point_t const p = geo::LineClosestPointWithUnitVectors(
    geo::origin(),                geo::Xaxis(),
    geo::origin() + geo::Yaxis(), (geo::Xaxis() + geo::Zaxis()) / std::sqrt(2.0)
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  
} // LineClosestPointWithUnitVectorsSimple45_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithUnitVectorsWithLocOnLines_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  std::pair<double, double> locOnLines { 0.0, 0.0 };
  geo::Point_t const p = geo::LineClosestPointWithUnitVectors(
    geo::origin()                - 3.0 * geo::Xaxis(), geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), geo::Yaxis(),
    &locOnLines
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  BOOST_CHECK_CLOSE(locOnLines.first,  +3.0, tol);
  BOOST_CHECK_CLOSE(locOnLines.second, -2.0, tol);
  
} // LineClosestPointWithUnitVectorsWithLocOnLines_test()


// =============================================================================
BOOST_AUTO_TEST_CASE(LineClosestPointTestCase) {
  
  LineClosestPointSimple_test();
  LineClosestPointSimple45_test();
  LineClosestPointWithLocOnLines_test();
  LineClosestPointWithScaledDirs_test();
  LineClosestPointWithNonHomogeneousDirs_test();
  
} // BOOST_AUTO_TEST_CASE(LineClosestPointTestCase)


// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(LineClosestPointWithUnitVectorsTestCase) {
  
  LineClosestPointWithUnitVectorsSimple_test();
  LineClosestPointWithUnitVectorsSimple45_test();
  LineClosestPointWithUnitVectorsWithLocOnLines_test();
  
} // BOOST_AUTO_TEST_CASE(LineClosestPointWithUnitVectorsTestCase)


// -----------------------------------------------------------------------------
