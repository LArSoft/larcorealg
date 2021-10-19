/**
 * @file    NumericUtils_test.cc
 * @brief   Unit test for NumericUtils functions.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 30, 2018
 * @see     NumericUtils.h
 */

// Boost libraries
#define BOOST_TEST_MODULE ( NumericUtils_test )
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/MetaUtils.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

// C/C++ standard libraries
#include <limits> // std::numeric_limits<>
#include <type_traits> // std::declval(), std::is_same<>, ...


//------------------------------------------------------------------------------
template <typename A, typename B, typename D = A>
void test_absDiff() {

  static_assert(std::is_same<
    decltype(std::declval<A>() - std::declval<B>()),
    decltype(std::declval<B>() - std::declval<A>())
    >(), "Difference between types is asymmetric."
    );

  A const a = 5;
  B const b = 6;

  auto absDeltaAB = util::absDiff(a, b);
  auto absDeltaBA = util::absDiff(b, a);

  static_assert(std::is_same<decltype(absDeltaAB), decltype(absDeltaBA)>(),
    "Results of |a-b| and |b-a| have different type!");

  static_assert(std::is_same<decltype(absDeltaAB), D>(),
    "Results of |a-b| has unexpected type!");

  BOOST_TEST(absDeltaAB == D(1));
  BOOST_TEST(absDeltaBA == D(1));

  // test that a very large B compares kind-of-correctly with a small A
  // (the difference needs to be large)
  B const m = std::numeric_limits<B>::max() - b;

  auto absDeltaAM = util::absDiff(a, m);
  auto absDeltaMA = util::absDiff(m, a);

  BOOST_TEST(absDeltaAM > D(2*(a + b)));
  BOOST_TEST(absDeltaMA > D(2*(a + b)));

} // test_absDiff()


template <typename A, typename B, typename D = std::add_const_t<A>>
void test_constexpr_absDiff() {

  static_assert(std::is_same<
    decltype(std::declval<A>() - std::declval<B>()),
    decltype(std::declval<B>() - std::declval<A>())
    >(), "Difference between types is asymmetric."
    );

  constexpr A a = 5;
  constexpr B b = 6;

  constexpr auto absDeltaAB = util::absDiff(a, b);
  constexpr auto absDeltaBA = util::absDiff(b, a);

  static_assert(std::is_same<decltype(absDeltaAB), decltype(absDeltaBA)>(),
    "Results of |a-b| and |b-a| have different type!");

  static_assert(std::is_same<decltype(absDeltaAB), D>(),
    "Results of |a-b| has unexpected type!");

  static_assert(absDeltaAB == 1, "|5-6| != 1?!?");
  static_assert(absDeltaBA == 1, "|6-5| != 1?!?");

} // test_constexpr_absDiff()


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(absDiffTestCase) {

  BOOST_TEST_INFO("Testing <int,int>");
  test_absDiff<int, int>();
  test_constexpr_absDiff<int, int>();

  BOOST_TEST_INFO("Testing <unsigned int,unsigned int>");
  test_absDiff<unsigned int, unsigned int>();
  test_constexpr_absDiff<unsigned int, unsigned int>();

} // BOOST_AUTO_TEST_CASE(absDiffTestCase)

//------------------------------------------------------------------------------

// fun with C!!
static_assert(5U - 6U >= 0,                         "ERROR 1");
static_assert(5U - 6  >= 0,                         "ERROR 2");
static_assert(5U < 6U,                              "ERROR 3");
static_assert(5U < 6,                               "ERROR 4");
// static_assert(5U > -6,                              "ERROR 5"); // !!!
static_assert(std::is_unsigned<decltype(5U - 6)>(), "ERROR 6");

//------------------------------------------------------------------------------
