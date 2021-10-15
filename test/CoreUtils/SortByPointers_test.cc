/**
 * @file    SortByPointers_test.cc
 * @brief   Unit test for `SortByPointers.h` utilities
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    September 28th, 2017
 * @see     SortByPointers.h
 */

// Boost libraries
#define BOOST_TEST_MODULE ( SortByPointers_test )
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/SortByPointers.h"

// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::sort()
#include <type_traits> // std::decay_t, std::is_same
#include <cstdlib> // std::abs()


//------------------------------------------------------------------------------
template <typename Data>
struct AbsSorter {

  bool operator()(Data a, Data b) const { return std::abs(a) < std::abs(b); }
  bool operator()(Data const* a, Data const* b) const
    { return this->operator()(*a, *b); }

}; // struct AbsSorter


void test_SortByPointers() {

  std::vector<int> data = { 8, -7, 5, 9, -2 };

  AbsSorter<int> absSorter;

  auto sortedData = data;
  std::sort(sortedData.begin(), sortedData.end(), absSorter);

  util::SortByPointers(data,
    [absSorter](auto& coll){ std::sort(coll.begin(), coll.end(), absSorter); }
    );

  BOOST_CHECK_EQUAL_COLLECTIONS
    (data.cbegin(), data.cend(), sortedData.cbegin(), sortedData.cend());

} // test_SortByPointers()


//------------------------------------------------------------------------------
void test_makePointerVector() {

  std::vector<int> data = { 8, -7, 5 };

  std::vector<int*> expectedDataPtr = { &(data[0]), &(data[1]), &(data[2]) };

  auto dataPtr = util::makePointerVector(data);

  static_assert(
    std::is_same<std::decay_t<decltype(dataPtr)>, std::vector<int*>>(),
      "Unexpected data type from makePointerVector()");

  BOOST_CHECK_EQUAL_COLLECTIONS(
    dataPtr.cbegin(), dataPtr.cend(),
    expectedDataPtr.cbegin(), expectedDataPtr.cend()
    );

} // test_makePointerVector()


//------------------------------------------------------------------------------
void test_MoveFromPointers() {

  std::vector<int> data = { 1, 2, 3, 4 };

  std::vector<int*> const dataPtr
    = { &(data[2]), &(data[3]), &(data[0]), &(data[1]) };
  std::vector<int> expectedMovedData = { 3, 4, 1, 2 };

  std::vector<int> movedData;
  util::MoveFromPointers(movedData, dataPtr);

  BOOST_CHECK_EQUAL_COLLECTIONS(
    movedData.cbegin(), movedData.cend(),
    expectedMovedData.cbegin(), expectedMovedData.cend()
    );

} // test_MoveFromPointers()


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(SortByPointers_testcase) {

  test_SortByPointers();

} // BOOST_AUTO_TEST_CASE(SortByPointers_testcase)


BOOST_AUTO_TEST_CASE(makePointerVector_testcase) {

  test_makePointerVector();

} // BOOST_AUTO_TEST_CASE(makePointerVector_testcase)


BOOST_AUTO_TEST_CASE(MoveFromPointers_testcase) {

  test_MoveFromPointers();

} // BOOST_AUTO_TEST_CASE(MoveFromPointers_testcase)

