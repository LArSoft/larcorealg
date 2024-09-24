/**
 * @file   get_elements_test.cc
 * @brief  Test of `util::get_elements()`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 13, 2019
 */

// testing library
#include "larcorealg/CoreUtils/get_elements.h"

// Boost libraries
#define BOOST_TEST_MODULE (get_elements_test)
#include <boost/test/unit_test.hpp>

// External libraries
#include "range/v3/view/zip.hpp"

// C/C++ libraries
#include <cstddef> // std::size_t
#include <tuple>
#include <utility> // std::as_const(), std::is_same_v
#include <vector>

// -----------------------------------------------------------------------------
void test_get_elements()
{
  struct A {};

  std::vector<std::tuple<int, void const*, A>> data{{0, nullptr, A{}}, {1, nullptr, A{}}};

  //
  // read-only access
  //
  {
    auto iItem = data.cbegin();
    for (auto&& [index, a] : util::get_elements<0U, 2U>(std::as_const(data))) {
      static_assert(std::is_same_v<decltype(index), int const&>);
      static_assert(std::is_same_v<decltype(a), A const&>);
      BOOST_TEST(index == std::get<0U>(*iItem));
      BOOST_TEST(&index == &std::get<0U>(*iItem));
      BOOST_TEST(&a == &std::get<2U>(*iItem));
      ++iItem;
    }
    BOOST_TEST((iItem == data.cend()));
  }

  //
  // read/write access
  //
  {
    auto iItem = data.cbegin();
    for (auto&& [index, a] : util::get_elements<0U, 2U>(data)) {
      static_assert(std::is_same_v<decltype(index), int&>);
      static_assert(std::is_same_v<decltype(a), A&>);
      BOOST_TEST(index == std::get<0U>(*iItem));
      BOOST_TEST(&index == &std::get<0U>(*iItem));
      BOOST_TEST(&a == &std::get<2U>(*iItem));
      ++iItem;
    }
    BOOST_TEST((iItem == data.cend()));
  }
}

void test_get_const_elements()
{
  struct A {};

  std::vector<std::tuple<int, void const*, A>> const data{{0, nullptr, A{}}, {1, nullptr, A{}}};

  //
  // read-only access
  //
  {
    auto iItem = data.cbegin();
    for (auto&& [index, a] : util::get_elements<0U, 2U>(std::as_const(data))) {
      static_assert(std::is_same_v<decltype(index), int const&>);
      static_assert(std::is_same_v<decltype(a), A const&>);
      BOOST_TEST(index == std::get<0U>(*iItem));
      BOOST_TEST(&index == &std::get<0U>(*iItem));
      BOOST_TEST(&a == &std::get<2U>(*iItem));
      ++iItem;
    }
    BOOST_TEST((iItem == data.cend()));
  }

  //
  // read/write access
  //
  {
    auto iItem = data.cbegin();
    for (auto&& [index, a] : util::get_elements<0U, 2U>(data)) {
      static_assert(std::is_same_v<decltype(index), int const&>);
      static_assert(std::is_same_v<decltype(a), A const&>);
      BOOST_TEST(index == std::get<0U>(*iItem));
      BOOST_TEST(&index == &std::get<0U>(*iItem));
      BOOST_TEST(&a == &std::get<2U>(*iItem));
      ++iItem;
    }
    BOOST_TEST((iItem == data.cend()));
  }
}

// -----------------------------------------------------------------------------
void test_get_elements_single()
{
  struct A {};

  std::vector<std::tuple<int, void const*, A>> data{{0, nullptr, A{}}, {1, nullptr, A{}}};

  //
  // read-only access
  //
  {
    auto iItem = data.cbegin();
    for (auto&& index : util::get_elements<0U>(std::as_const(data))) {
      static_assert(std::is_same_v<decltype(index), int const&>);
      BOOST_TEST(index == std::get<0U>(*iItem));
      BOOST_TEST(&index == &std::get<0U>(*iItem));
      ++iItem;
    }
    BOOST_TEST((iItem == data.cend()));
  }

  //
  // read/write access
  //
  {
    auto iItem = data.cbegin();
    for (auto&& index : util::get_elements<0U>(data)) {
      static_assert(std::is_same_v<decltype(index), int&>);
      BOOST_TEST(index == std::get<0U>(*iItem));
      BOOST_TEST(&index == &std::get<0U>(*iItem));
      ++iItem;
    }
    BOOST_TEST((iItem == data.cend()));
  }
}

void test_get_const_elements_single()
{
  struct A {};

  std::vector<std::tuple<int, void const*, A>> const data{{0, nullptr, A{}}, {1, nullptr, A{}}};

  //
  // read-only access
  //
  {
    auto iItem = data.cbegin();
    for (auto&& index : util::get_elements<0U>(std::as_const(data))) {
      static_assert(std::is_same_v<decltype(index), int const&>);
      BOOST_TEST(index == std::get<0U>(*iItem));
      BOOST_TEST(&index == &std::get<0U>(*iItem));
      ++iItem;
    }
    BOOST_TEST((iItem == data.cend()));
  }

  //
  // read/write access
  //
  {
    auto iItem = data.cbegin();
    for (auto&& index : util::get_elements<0U>(data)) {
      static_assert(std::is_same_v<decltype(index), int const&>);
      BOOST_TEST(index == std::get<0U>(*iItem));
      BOOST_TEST(&index == &std::get<0U>(*iItem));
      ++iItem;
    }
    BOOST_TEST((iItem == data.cend()));
  }
}

// -----------------------------------------------------------------------------
void get_elements_documentation_test()
{
  /*
   * The promise:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * std::vector<std::tuple<char, int, float>> data
   *   { { 'z', 0, 1.0F }, { 'o', 1, 3.0F }, { 't', 2, 9.0F } };
   * std::map<char, double> factors;
   *
   * for (auto const& [ letter, factor ]: util::get_elements<0U, 2U>(data)) {
   *   factors.emplace(letter, factor);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * and
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * std::vector<int> exponents;
   *
   * for (int exponent: util::get_elements<1U>(data)) {
   *   exponents.push_back(exponent);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  std::vector<std::tuple<char, int, float>> data{{'z', 0, 1.0F}, {'o', 1, 3.0F}, {'t', 2, 9.0F}};
  std::map<char, double> factors;

  for (auto const& [letter, factor] : util::get_elements<0U, 2U>(data)) {
    factors.emplace(letter, factor);
  }

  std::vector<int> exponents;

  for (int exponent : util::get_elements<1U>(data)) {
    exponents.push_back(exponent);
  }

  //
  // check
  //
  std::map<char, double> const expected_factors{{'z', 1.0F}, {'o', 3.0F}, {'t', 9.0F}};
  std::vector<int> const expected_exponents{0, 1, 2};

  BOOST_TEST(factors.size() == expected_factors.size());
  for (auto const& [f, expected_f] : ranges::views::zip(factors, expected_factors)) {
    BOOST_TEST(f.first == expected_f.first);
    BOOST_TEST(f.second == expected_f.second);
  }

  BOOST_TEST(exponents == expected_exponents);
}

// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(get_elements_testcase)
{
  test_get_elements();
  test_get_const_elements();
  test_get_elements_single();
  test_get_const_elements_single();

  get_elements_documentation_test();
}

// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
