/**
 * @file    RealComparisons_test.cc
 * @brief   Unit test for `RealComparisons` class
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    November 4th, 2016
 * @see     larcorealg/CoreUtils/RealComparisons.h
 */

// Boost libraries
#define BOOST_TEST_MODULE ( RealComparisons_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_TEST()

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"

// C/C++ standard libraries
#include <ostream>
#include <cmath>


BOOST_AUTO_TEST_CASE(test_RealComparisons) {

  // instantiate a RealComparisons with two classes
  constexpr lar::util::RealComparisons check(1e-5f);

  constexpr double const sqrt2 = std::sqrt(2);

  // check zero()
  constexpr double const epsilon = 2.0 - (sqrt2 * sqrt2);
  BOOST_TEST( check.zero(epsilon));
  BOOST_TEST( check.zero(0.0));
  BOOST_TEST( check.zero(1e-5));
  BOOST_TEST( check.zero(-1e-5));
  BOOST_TEST(!check.zero(1.01e-5));
  BOOST_TEST(!check.zero(-1.01e-5 * 1.01));

  // check nonZero()
  BOOST_TEST(!check.nonZero(epsilon));
  BOOST_TEST(!check.nonZero(0.0));
  BOOST_TEST(!check.nonZero(1e-5));
  BOOST_TEST(!check.nonZero(-1e-5));
  BOOST_TEST( check.nonZero(1.01e-5));
  BOOST_TEST( check.nonZero(-1.01e-5 * 1.01));

  // check equal()
  BOOST_TEST(!check.equal(sqrt2, 1.4142));
  BOOST_TEST( check.nonEqual(sqrt2, 1.4142));
  BOOST_TEST( check.equal(sqrt2, 1.414213));
  BOOST_TEST(!check.nonEqual(sqrt2, 1.414213));

  // check strictlyNegative()
  BOOST_TEST(!check.strictlyNegative(+1e-5 + 1e-7)); // outside tolerance
  BOOST_TEST(!check.strictlyNegative(+1e-5 - 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyNegative(0.0));          // zero
  BOOST_TEST(!check.strictlyNegative(-1e-5 + 1e-7)); // within tolerance
  BOOST_TEST( check.strictlyNegative(-1e-5 - 1e-7)); // outside tolerance

  // check strictlyPositive()
  BOOST_TEST( check.strictlyPositive(+1e-5 + 1e-7)); // outside tolerance
  BOOST_TEST(!check.strictlyPositive(+1e-5 - 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyPositive(0.0));          // zero
  BOOST_TEST(!check.strictlyPositive(-1e-5 + 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyPositive(-1e-5 - 1e-7)); // outside tolerance

  // check nonNegative()
  BOOST_TEST( check.nonNegative(+1e-5 + 1e-7)); // outside tolerance
  BOOST_TEST( check.nonNegative(+1e-5 - 1e-7)); // within tolerance
  BOOST_TEST( check.nonNegative(0.0));          // zero
  BOOST_TEST( check.nonNegative(-1e-5 + 1e-7)); // within tolerance
  BOOST_TEST(!check.nonNegative(-1e-5 - 1e-7)); // outside tolerance

  // check nonPositive()
  BOOST_TEST(!check.nonPositive(+1e-5 + 1e-7)); // outside tolerance
  BOOST_TEST( check.nonPositive(+1e-5 - 1e-7)); // within tolerance
  BOOST_TEST( check.nonPositive(0.0));          // zero
  BOOST_TEST( check.nonPositive(-1e-5 + 1e-7)); // within tolerance
  BOOST_TEST( check.nonPositive(-1e-5 - 1e-7)); // outside tolerance

  // check strictlySmaller()
  BOOST_TEST(!check.strictlySmaller(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_TEST(!check.strictlySmaller(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlySmaller(1.0, 1.0));        // equal
  BOOST_TEST(!check.strictlySmaller(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_TEST( check.strictlySmaller(1.0, 1.0 + 1e-4)); // outside tolerance

  // check nonSmaller()
  BOOST_TEST( check.nonSmaller(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_TEST( check.nonSmaller(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_TEST( check.nonSmaller(1.0, 1.0));        // equal
  BOOST_TEST( check.nonSmaller(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_TEST(!check.nonSmaller(1.0, 1.0 + 1e-4)); // outside tolerance

  // check strictlyGreater()
  BOOST_TEST( check.strictlyGreater(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_TEST(!check.strictlyGreater(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyGreater(1.0, 1.0));        // equal
  BOOST_TEST(!check.strictlyGreater(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyGreater(1.0, 1.0 + 1e-4)); // outside tolerance

  // check nonGreater()
  BOOST_TEST(!check.nonGreater(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_TEST( check.nonGreater(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_TEST( check.nonGreater(1.0, 1.0));        // equal
  BOOST_TEST( check.nonGreater(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_TEST( check.nonGreater(1.0, 1.0 + 1e-4)); // outside tolerance

  // check within()
  BOOST_TEST(!check.within(sqrt2, 0., 1.41420));
  BOOST_TEST( check.within(sqrt2, 0., 1.41421));
  BOOST_TEST( check.within(sqrt2, 1.41422, 2.));
  BOOST_TEST(!check.within(sqrt2, 1.41423, 2.));

  // check inverted limits
  BOOST_TEST(!check.within(sqrt2, 1.41421, 0.));
  BOOST_TEST(check.withinSorted(sqrt2, 1.41421, 0.));

} // BOOST_AUTO_TEST_CASE(test_RealComparisons)
