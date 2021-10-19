/**
 * @file   enumerate_test.cc
 * @brief  Test of `util::enumerate` and support utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 14, 2019
 *
 */


// testing library
#include "larcorealg/CoreUtils/enumerate.h"

// Boost libraries
#define BOOST_TEST_MODULE ( enumerate_test )
#include <boost/test/unit_test.hpp>

// C/C++ libraries
#include <vector>
#include <array>
#include <cstddef> // std::size_t
#include <cassert>


// -----------------------------------------------------------------------------
void test_enumerate() {

  //
  // standard use case with zipping included
  //

  // prepare the data set
  constexpr std::size_t N = 7U;

  std::array<int, N> twice;
  std::vector<double> thrice(N + 1);

  for (std::size_t i = 0; i < N; ++i) {
    twice[i] = 2 * i;
    thrice[i] = 3.0 * i;
  } // for
  thrice[N] = 3.0 * N;

  //
  // iteration using the first element as lead
  //
  unsigned int iLoop = 0;
  for (auto&& [i, a, b]: util::enumerate(twice, thrice)) {

    BOOST_TEST(i == iLoop);

    BOOST_TEST(a == twice[iLoop]);
    BOOST_TEST(&a == &(twice[iLoop]));

    BOOST_TEST(b == thrice[iLoop]);
    BOOST_TEST(&b == &thrice[iLoop]);

    ++iLoop;
  } // for
  BOOST_TEST(iLoop == twice.size());

  //
  // iteration using the second element as lead
  //

  // (make the second object shorter, so that we can check easily it leads)
  thrice.pop_back();
  thrice.pop_back();
  assert(thrice.size() == N - 1);

  iLoop = 0;
  for (auto&& [i, a, b]: util::enumerate<1>(twice, thrice)) {

    BOOST_TEST(i == iLoop);

    BOOST_TEST(a == twice[iLoop]);
    BOOST_TEST(&a == &(twice[iLoop]));

    BOOST_TEST(b == thrice[iLoop]);
    BOOST_TEST(&b == &thrice[iLoop]);

    ++iLoop;
  } // for
  BOOST_TEST(iLoop == thrice.size());


} // test_enumerate()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(enumerate_testcase) {

  test_enumerate();

} // BOOST_AUTO_TEST_CASE(enumerate_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
