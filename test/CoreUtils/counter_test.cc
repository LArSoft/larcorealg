/**
 * @file   counter_test.cc
 * @brief  Test of `util::counter()` and support utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 14, 2019
 *
 */

// testing library
#include "larcorealg/CoreUtils/counter.h"

// Boost libraries
#define BOOST_TEST_MODULE ( counter_test )
#include <boost/test/unit_test.hpp>

// C/C++ libraries
#include <vector>
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
void test_count_iterator_documentation() {

  // the promise:
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * std::vector<int> data;
   * for (util::count_iterator it;; ++it) {
   *   if (*it >= 10) break;
   *   data.push_back(*it); // implicit conversion `std::size_t` to `int`
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * this infinite loop will print `0` on the first iteration, `1` on the second
   * and so forth, until at the eleventh iteration, just before it can print
   * `10`, the loop is forcibly broken.
   *
   */
  std::vector<int> data;
  // BUG Clang 5.0 can't correctly apply explicit (or implicit) deduction rules; Clang 7.0 fixes it
  // for (util::count_iterator it;; ++it) {
  for (util::count_iterator<> it;; ++it) {
    if (*it >= 10) break;
    data.push_back(*it); // implicit conversion `std::size_t` to `int`
  }

  // the test:
  BOOST_TEST(data.size() == 10U);
  for (std::size_t i = 0; i < data.size(); ++i) {

    BOOST_TEST(data[i] == (int) i);

  } // for

} // test_count_iterator_documentation()


//------------------------------------------------------------------------------
void test_counter_documentation() {

  // the promise:
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<std::size_t> data;
   * for (auto i: util::counter(4, 8)) {
   *   data.push_back(i);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will insert in `data` the numbers from `4` to `7`, just like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * for (std::size_t i = 4; i < 8; ++i) {
   *   data.push_back(i);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * would.
   */
  std::vector<std::size_t> data;
  for (auto i: util::counter(4, 8)) {
    data.push_back(i);
  }

  // the test:
  std::vector<std::size_t> control_data;
  for (std::size_t i = 4; i < 8; ++i) {
    control_data.push_back(i);
  }

  BOOST_TEST(data.size() == control_data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    BOOST_TEST(data[i] == control_data[i]);
  }

} // test_counter_documentation()


// -----------------------------------------------------------------------------
void test_infinite_counter_documentation() {

  // the promise:
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<unsigned char> data;
   * for (auto ch: util::infinite_counter<unsigned char>()) {
   *   if (data.size() >= 512U) break;
   *   data.push_back(ch);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  std::vector<unsigned char> data;
  for (auto ch: util::infinite_counter<unsigned char>()) {
    if (data.size() >= 512U) break;
    data.push_back(ch);
  }

  // the test:
  constexpr std::size_t N = 1U << 8 * sizeof(unsigned char);
  static_assert(N == 256U); // just in case...

  BOOST_TEST(data.size() == N * 2);
  for (std::size_t i = 0; i < data.size(); ++i) {

    BOOST_TEST(data[i]  == static_cast<unsigned char>(i % N));

  } // for

} // test_infinite_counter_documentation()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(counter_testcase) {

  test_count_iterator_documentation();
  test_counter_documentation();
  test_infinite_counter_documentation();

} // BOOST_AUTO_TEST_CASE(counter_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
