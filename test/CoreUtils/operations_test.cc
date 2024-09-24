/**
 * @file   larcorealg/test/CoreUtils/operations_test.cc
 * @brief  Unit test for `larcorealg/CoreUtils/operations.h`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 5, 2019
 * @see    `larcorealg/CoreUtils/operations.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE (operations_test)
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"
#include "larcorealg/CoreUtils/operations.h"

// External libraries
#include "range/v3/view/zip.hpp"

// C/C++ standard libraries
#include <algorithm>
#include <memory>  // std::unique_ptr<>
#include <numeric> // std::iota()
#include <vector>

//------------------------------------------------------------------------------
void test_AddressTaker_documentation()
{
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int*> ptrs(data.size());
   * std::transform
   *   (data.begin(), data.end(), ptrs.begin(), util::AddressTaker{});
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  std::vector<int*> ptrs(data.size());
  std::transform(data.begin(), data.end(), ptrs.begin(), util::AddressTaker{});

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the test
  BOOST_TEST(ptrs.size() == data.size());
  for (auto&& [value, ptr] : ranges::views::zip(data, ptrs)) {
    BOOST_TEST(ptr);
    BOOST_TEST(*ptr == value);
    BOOST_TEST(ptr == &value);
  }
}

//------------------------------------------------------------------------------
void test_takeAddress_documentation()
{
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);

  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int*> ptrs(data.size());
   * std::transform
   *   (data.begin(), data.end(), ptrs.begin(), util::takeAddress());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  std::vector<int*> ptrs(data.size());
  std::transform(data.begin(), data.end(), ptrs.begin(), util::takeAddress());

  // the test
  BOOST_TEST(ptrs.size() == data.size());
  for (auto&& [value, ptr] : ranges::views::zip(data, ptrs)) {
    BOOST_TEST(ptr);
    BOOST_TEST(*ptr == value);
    BOOST_TEST(ptr == &value);
  }
}

//------------------------------------------------------------------------------
void test_takeAddress()
{
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);

  std::vector<int const*> dataPtr;
  std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtr), util::takeAddress());

  for (auto&& [value, ptr] : ranges::views::zip(data, dataPtr)) {
    BOOST_TEST(ptr);
    BOOST_TEST(*ptr == value);
    BOOST_TEST(ptr == &value);
  }
}

//------------------------------------------------------------------------------
void test_takeAddress_whyBother()
{
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);

  std::vector<int const*> dataPtr;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // `std::addressof()` approach
  using addressof_t = int const* (*)(int const&);

  std::transform(
    data.cbegin(), data.cend(), std::back_inserter(dataPtr), ((addressof_t)&std::addressof));

  for (auto&& [value, ptr] : ranges::views::zip(data, dataPtr)) {
    BOOST_TEST(ptr);
    BOOST_TEST(*ptr == value);
    BOOST_TEST(ptr == &value);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // lambda approach
  dataPtr.clear();

  auto takeAddress = [](auto&& ref) { return std::addressof(ref); };

  std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtr), takeAddress);

  for (auto&& [value, ptr] : ranges::views::zip(data, dataPtr)) {
    BOOST_TEST(ptr);
    BOOST_TEST(*ptr == value);
    BOOST_TEST(ptr == &value);
  }
}

//------------------------------------------------------------------------------
void test_Dereferencer_documentation()
{
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);

  std::vector<int*> ptrs(data.size());
  std::transform(data.begin(), data.end(), ptrs.begin(), util::takeAddress());

  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> values(ptrs.size());
   * std::transform
   *   (ptrs.cbegin(), ptrs.cend(), values.begin(), util::Dereferencer{});
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  std::vector<int> values(ptrs.size());
  std::transform(ptrs.cbegin(), ptrs.cend(), values.begin(), util::Dereferencer{});

  // the test
  BOOST_TEST(values.size() == data.size());
  for (auto&& [value, orig] : ranges::views::zip(data, values)) {
    BOOST_TEST(value == orig);
  }
}

//------------------------------------------------------------------------------
void test_dereference_documentation()
{
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);

  std::vector<int*> ptrs(data.size());
  std::transform(data.begin(), data.end(), ptrs.begin(), util::takeAddress());

  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> values(ptrs.size());
   * std::transform
   *   (ptrs.cbegin(), ptrs.cend(), values.begin(), util::dereference());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  std::vector<int> values(ptrs.size());
  std::transform(ptrs.cbegin(), ptrs.cend(), values.begin(), util::dereference());

  // the test
  BOOST_TEST(values.size() == data.size());
  for (auto&& [value, orig] : ranges::views::zip(data, values)) {
    BOOST_TEST(value == orig);
  }
}

//------------------------------------------------------------------------------
void test_dereference_C_ptr()
{
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);

  std::vector<int const*> dataPtrs;
  std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtrs), util::takeAddress());

  std::vector<int> dataAgain;
  std::transform(
    dataPtrs.cbegin(), dataPtrs.cend(), std::back_inserter(dataAgain), util::dereference());

  BOOST_TEST(dataAgain.size() == data.size());
  for (auto&& [value, valueAgain] : ranges::views::zip(data, dataAgain)) {
    BOOST_TEST(valueAgain == value);
    BOOST_TEST(&valueAgain != &value);
  }
}

//------------------------------------------------------------------------------
void test_dereference_unique_ptr()
{
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);

  std::vector<std::unique_ptr<int>> dataPtrs;
  dataPtrs.reserve(data.size());
  for (auto&& value : data)
    dataPtrs.push_back(std::make_unique<int>(value));

  std::vector<int> dataAgain;
  std::transform(
    dataPtrs.cbegin(), dataPtrs.cend(), std::back_inserter(dataAgain), util::dereference());

  BOOST_TEST(dataAgain.size() == data.size());
  for (auto&& [value, valueAgain] : ranges::views::zip(data, dataAgain)) {
    BOOST_TEST(valueAgain == value);
    BOOST_TEST(&valueAgain != &value);
  }
}

//------------------------------------------------------------------------------
void test_dereference_uncopiable()
{
  struct ToughInt : private lar::UncopiableAndUnmovableClass {
    int value = 0;
  };

  ToughInt value;
  ToughInt const* pValue = &value;

  BOOST_TEST(pValue == &(util::dereference()(pValue)));
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(takeAddress_testcase)
{
  test_takeAddress();
  test_takeAddress_documentation();
  test_AddressTaker_documentation();
  test_takeAddress_whyBother();
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(dereference_testcase)
{
  test_dereference_C_ptr();
  test_dereference_unique_ptr();
  test_dereference_uncopiable();
  test_dereference_documentation();
  test_Dereferencer_documentation();
}
