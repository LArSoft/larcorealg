/**
 * @file    DebugUtils_test.cc
 * @brief   Unit test for some of the utilities in `DebugUtils.h`.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 22, 2019
 * @see     `larcorealg/CoreUtils/DebugUtils.h`
 *
 * Test of backtrace print functions is on its own: `printBacktrace_test.cc`.
 */

// Boost libraries
#define BOOST_TEST_MODULE ( DebugUtils_test )
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/DebugUtils.h"

// C/C++ standard libraries
#include <type_traits>
#include <array>


//------------------------------------------------------------------------------
//--- static tests (would fail at compile time)
//------------------------------------------------------------------------------

/*
 * This is not a real test, but rather an example.
 * If `DoTest` is `true`, compilation will fail.
 *
 * The purpose of the example is to see what is the exact type `element_type`
 * of the collection type passed to `OurClass`, but only when the collection
 * type is not constant.
 *
 * In addition, `DoTest` needs to be changed to `true` for this to actually
 * happen.
 */

template <typename Coll>
struct OurClass {

  using Collection_t = Coll;

  using value_type = typename Collection_t::element_type;

  lar::debug::static_assert_on
    <value_type, std::is_const_v<std::remove_reference_t<Coll>>>
    debugVar;

}; // struct OurClass


void static_assert_on_test() {

  // this should never trigger a static assertion failure:
  (void) OurClass<std::unique_ptr<double>>();

  // this triggers a static assertion failure (#if 1, which we disabled)
#if 0
  (void) OurClass<std::unique_ptr<int[10]> const>();
#endif

} // static_assert_on_test()


//------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ReferencesTestCase) {

  static_assert_on_test();

} // BOOST_AUTO_TEST_CASE(ReferencesTestCase)

//------------------------------------------------------------------------------
