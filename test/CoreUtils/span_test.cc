/**
 * @file   span_test.cc
 * @brief  Unit test for `util::span` class.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/CoreUtils/span.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE (span_test)
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/operations.h"
#include "larcorealg/CoreUtils/span.h"
#include "larcorealg/CoreUtils/zip.h"

#include "range/v3/view/enumerate.hpp"
#include "range/v3/view/zip.hpp"

// C/C++ standard libraries
#include <iostream>    // std::cout
#include <memory>      // std::unique_ptr<>
#include <numeric>     // std::iota()
#include <type_traits> // std::is_same<>
#include <vector>

//------------------------------------------------------------------------------
template <typename Iter, typename Cont>
void test_span(Cont& v)
{
  auto r = util::make_span(v);

  // check the type of the object
  static_assert(std::is_same<decltype(r), util::span<Iter>>());
  static_assert(std::is_same<decltype(r), util::span<Iter, Iter>>());
  static_assert(std::is_same<typename decltype(r)::value_type, typename Iter::value_type>());
  static_assert(std::is_same<typename decltype(r)::reference, typename Iter::reference>());

  BOOST_TEST((r.begin() == v.begin()));
  BOOST_TEST((r.end() == v.end()));

  BOOST_TEST(r.empty() == v.empty());
  BOOST_TEST(r.size() == v.size());

  auto iV = v.cbegin();
  unsigned int n = 0;
  for (auto i : r) {
    BOOST_TEST(i == *iV);
    ++n;
    ++iV;
  }
  BOOST_TEST(n == v.size());
}

//------------------------------------------------------------------------------
template <typename Iter, typename Cont>
void test_const_span(Cont& v)
{
  auto r = util::make_const_span(v);

  // check the type of the object
  static_assert(std::is_same<decltype(r), util::span<Iter>>());
  static_assert(std::is_same<decltype(r), util::span<Iter, Iter>>());
  static_assert(std::is_same<typename decltype(r)::value_type, typename Iter::value_type>());
  static_assert(std::is_same<typename decltype(r)::reference, typename Iter::reference>());

  BOOST_TEST((r.begin() == v.cbegin()));
  BOOST_TEST((r.end() == v.cend()));

  BOOST_TEST(r.empty() == v.empty());
  BOOST_TEST(r.size() == v.size());

  auto iV = v.cbegin();
  unsigned int n = 0;
  for (auto i : r) {
    BOOST_TEST(i == *iV);
    ++n;
    ++iV;
  }
  BOOST_TEST(n == v.size());
}

//------------------------------------------------------------------------------
template <typename Iter, typename Cont>
void test_adapted_span(Cont& v)
{
  /*
   * In this test we create a container of pointers to `v` elements and iterate
   * `v` values through that one. The iteration is supposed to be transparent,
   * with the iteration code not expressing that we are passing via pointers.
   */

  std::vector<typename Iter::pointer> ptrs;
  std::transform(v.begin(), v.end(), std::back_inserter(ptrs), [](auto& v) { return &v; });

  auto r =
    util::make_adapted_span(ptrs, [](auto&& iter) { return boost::make_indirect_iterator(iter); });

  // check the type of the object
  static_assert(std::is_same<typename decltype(r)::value_type, typename Iter::value_type>());

  BOOST_TEST(r.empty() == v.empty());
  BOOST_TEST(r.size() == v.size());

  unsigned int n = 0;
  for (auto&& [rangeValue, value] : ranges::views::zip(r, v)) {
    BOOST_TEST(&rangeValue == &value);
    BOOST_TEST(rangeValue == value);
    ++n;
  }
  BOOST_TEST(n == v.size());
}

//------------------------------------------------------------------------------
template <typename Iter, typename Cont>
void test_transformed_span(Cont& v)
{
  /*
   * In this test we create a container of pointers to `v` elements and iterate
   * `v` values through that one. The iteration is supposed to be transparent,
   * with the iteration code not expressing that we are passing via pointers,
   * because of our transformation dereferencing the pointers.
   */

  using pointer_t = typename std::iterator_traits<Iter>::pointer;
  using reference_t = typename std::iterator_traits<Iter>::reference;

  std::vector<pointer_t> ptrs;
  std::transform(v.begin(), v.end(), std::back_inserter(ptrs), [](auto& v) { return &v; });

  auto r = util::make_transformed_span(ptrs, [](auto* ptr) -> reference_t { return *ptr; });

  // check the type of the object
  static_assert(std::is_same<typename decltype(r)::value_type, typename Iter::value_type>());

  BOOST_TEST(r.empty() == v.empty());
  BOOST_TEST(r.size() == v.size());

  unsigned int n = 0;
  for (auto&& [rangeValue, value] : util::zip(r, v)) {
    BOOST_TEST(&rangeValue == &value);
    BOOST_TEST(rangeValue == value);
    ++n;
  }
  BOOST_TEST(n == v.size());
}

//------------------------------------------------------------------------------
template <typename Iter, typename Cont>
void test_transformed_span_with_unmoveable_values(Cont& v)
{
  /*
   * In this test we create a container of unique_ptr to `v` elements and
   * iterate `v` values through it. The iteration is supposed to be transparent,
   * with the iteration code not expressing that we are passing via pointers,
   * because of our transformation dereferencing the pointers.
   */

  using value_t = typename std::iterator_traits<Iter>::value_type;
  using pointer_t = std::unique_ptr<value_t>;

  std::vector<pointer_t> ptrs;
  std::transform(v.begin(), v.end(), std::back_inserter(ptrs), [](auto& v) {
    return std::make_unique<value_t>(v);
  });

  auto r = util::make_transformed_span(ptrs, util::dereference());

  // check the type of the object
  static_assert(std::is_same<typename decltype(r)::value_type, value_t>());

  BOOST_TEST(r.empty() == v.empty());
  BOOST_TEST(r.size() == v.size());

  unsigned int n = 0;
  for (auto&& [rangeValue, value] : ranges::views::zip(r, v)) {
    BOOST_TEST(&rangeValue == &value);
    BOOST_TEST(rangeValue == value);
    ++n;
  }
  BOOST_TEST(n == v.size());
}

//------------------------------------------------------------------------------

struct SpanDocumentation1TestClass {
  using Data_t = std::vector<int>;

  int sum(util::span<Data_t::iterator> r)
  {
    decltype(r)::value_type s = 0;
    for (auto v : r)
      s += v;
    return s;
  }

  void analyse()
  {
    Data_t v(10U);
    std::iota(v.begin(), v.end(), 1); // { 1, 2, 3, 4 ,5 ,6, 7, 8, 9, 10 }

    auto span5 = util::span(v.begin(), v.begin() + 5);
    std::cout << "Sum of 5 numbers: " << sum(span5) << std::endl;

    auto span8 = util::span(v.begin(), v.begin() + 8);
    std::cout << "Sum of 8 numbers: " << sum(span8) << std::endl;

    auto full_span = util::make_span(v);
    std::cout << "Sum of all numbers: " << sum(full_span) << std::endl;
  }

}; // struct SpanDocumentation1TestClass

//------------------------------------------------------------------------------

struct makeAdaptedSpanDocumentation1TestClass {

  float accumulate(std::vector<std::unique_ptr<float>> const& v)
  {
    using src_iterator = std::vector<std::unique_ptr<float>>::const_iterator;

    float sum = 0.0F;
    for (float v : util::make_adapted_span(v, boost::make_indirect_iterator<src_iterator>))
      sum += v;

    return sum;
  }

  void analyse()
  {
    float const salt = 3.0;
    constexpr std::size_t N = 10;

    std::vector<std::unique_ptr<float>> data;
    for (auto i : util::counter(N))
      data.push_back(std::make_unique<float>(salt * i));

    float sum = accumulate(data);
    float const expectedSum = (N * (N - 1) / 2) * salt;

    BOOST_TEST(sum == expectedSum);
  }

}; // struct makeAdaptedSpanDocumentation1TestClass

//------------------------------------------------------------------------------

struct makeTransformedSpanDocumentation1TestClass {

  float accumulate(std::vector<std::unique_ptr<float>> const& v)
  {
    float sum = 0.0F;
    for (float v : util::make_transformed_span(v, [](auto& ptr) { return *ptr; }))
      sum += v;

    return sum;
  }

  void scale(std::vector<std::unique_ptr<float>>& v, float factor)
  {
    for (float& v : util::make_transformed_span(v, [](auto& ptr) -> float& { return *ptr; }))
      v *= factor;
  }

  void analyse_accumulate()
  {
    float const salt = 3.0;
    constexpr std::size_t N = 10;

    std::vector<std::unique_ptr<float>> data;
    for (auto i : util::counter(N))
      data.push_back(std::make_unique<float>(salt * i));

    float sum = accumulate(data);
    float const expectedSum = (N * (N - 1) / 2) * salt;

    BOOST_TEST(sum == expectedSum);
  }

  void analyse_scale()
  {
    constexpr std::size_t N = 10;
    float const factor = 3.0F;

    std::vector<std::unique_ptr<float>> data;
    for (auto i : util::counter(N))
      data.push_back(std::make_unique<float>(i));

    scale(data, factor);

    for (auto&& [i, ptr] : data | ranges::views::enumerate)
      BOOST_TEST(*ptr == factor * i);
  }

  void analyse()
  {
    analyse_accumulate();
    analyse_scale();
  }

}; // struct makeTransformedSpanDocumentation1TestClass

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(span_testcase)
{
  using TestVector_t = std::vector<int>;

  TestVector_t ev;
  test_span<TestVector_t::iterator>(ev);
  test_const_span<TestVector_t::const_iterator>(ev);

  TestVector_t v3{1, 2, 3};
  test_span<TestVector_t::iterator>(v3);
  test_const_span<TestVector_t::const_iterator>(v3);

  TestVector_t const cv4{1, 2, 3, 4};
  test_span<TestVector_t::const_iterator>(cv4);
  test_const_span<TestVector_t::const_iterator>(cv4);
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(adapted_span_testcase)
{
  using TestVector_t = std::vector<int>;

  TestVector_t ev;
  test_adapted_span<TestVector_t::iterator>(ev);

  TestVector_t v3{1, 2, 3};
  test_adapted_span<TestVector_t::iterator>(v3);

  TestVector_t const cv4{1, 2, 3, 4};
  test_adapted_span<TestVector_t::const_iterator>(cv4);
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(transformed_span_testcase)
{
  using TestVector_t = std::vector<int>;

  TestVector_t ev;
  test_transformed_span<TestVector_t::iterator>(ev);
  test_transformed_span_with_unmoveable_values<TestVector_t::iterator>(ev);

  TestVector_t v3{1, 2, 3};
  test_transformed_span<TestVector_t::iterator>(v3);

  TestVector_t const cv4{1, 2, 3, 4};
  test_transformed_span<TestVector_t::const_iterator>(cv4);
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(span_documentation_testcase)
{
  // this only checks that the example compiles and runs
  SpanDocumentation1TestClass doc1;
  doc1.analyse();
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(adapted_span_documentation_testcase)
{
  makeAdaptedSpanDocumentation1TestClass doc1;
  doc1.analyse();
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(transformed_span_documentation_testcase)
{
  makeTransformedSpanDocumentation1TestClass doc1;
  doc1.analyse();
}

//------------------------------------------------------------------------------
