/**
 * @file   geometrydatacontainers_test.cxx
 * @brief  Unit test for GeometryDataContainers.h library.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 2nd, 2018
 */

// Boost libraries
#define BOOST_TEST_MODULE (geometry data containers test)
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::static_assert_on()
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/Geometry/GeometryDataContainers.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

//------------------------------------------------------------------------------
template <typename T>
struct Summer {

  T sum = T{0};

  void operator()(T v) { sum += v; }

  T get() const { return sum; }
  void reset() { sum = T{0}; }

}; // struct Summer

//------------------------------------------------------------------------------
void TPCDataContainerTest(geo::TPCDataContainer<int> data, // copy here is intentional
                          std::size_t const NCryostats,
                          std::size_t const NTPCs)
{
  std::size_t const N = NCryostats * NTPCs;

  static_assert(data.dimensions() == 2U);
  BOOST_TEST(data.dimSize<0U>() == NCryostats);
  BOOST_TEST(data.dimSize<1U>() == NTPCs);
  BOOST_TEST(data.dimSize<2U>() == 0U);
  BOOST_TEST(data.dimSize<3U>() == 0U);

  BOOST_TEST(!data.empty());
  BOOST_TEST(data.size() == N);
  BOOST_TEST(data.capacity() >= N);

  for (auto c : util::counter<unsigned int>(NCryostats))
    for (auto t : util::counter<unsigned int>(NTPCs))
      BOOST_TEST((data[{c, t}]) == 0);

  BOOST_TEST(data.firstID() == geo::TPCID(0, 0));
  BOOST_TEST(data.lastID() == geo::TPCID(1, 2));

  std::size_t expected_index = 0U;

  // simple R/W iteration test
  for (auto& value : data) {
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    geo::TPCID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);

  // ID/data pair R/W iteration test
  expected_index = 0U;
  for (auto&& [ID, value] : data.items()) {
    static_assert(std::is_same_v<decltype(ID), geo::TPCID>);
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    geo::TPCID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);

  BOOST_TEST(data.hasTPC({0, 0}));
  BOOST_TEST(data.hasTPC({0, 1}));
  BOOST_TEST(data.hasTPC({0, 2}));
  BOOST_TEST(!data.hasTPC({0, 3}));
  BOOST_TEST(!data.hasTPC({0, 4}));
  BOOST_TEST(data.hasTPC({1, 0}));
  BOOST_TEST(data.hasTPC({1, 1}));
  BOOST_TEST(data.hasTPC({1, 2}));
  BOOST_TEST(!data.hasTPC({1, 3}));
  BOOST_TEST(!data.hasTPC({1, 4}));
  BOOST_TEST(!data.hasTPC({2, 0}));
  BOOST_TEST(!data.hasTPC({2, 1}));
  BOOST_TEST(!data.hasTPC({2, 2}));
  BOOST_TEST(!data.hasTPC({2, 3}));
  BOOST_TEST(!data.hasTPC({2, 4}));

  BOOST_TEST(data.hasCryostat(geo::TPCID{0, 0}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{0, 1}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{0, 2}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{0, 3}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{0, 4}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{1, 0}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{1, 1}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{1, 2}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{1, 3}));
  BOOST_TEST(data.hasCryostat(geo::TPCID{1, 4}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{2, 0}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{2, 1}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{2, 2}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{2, 3}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{2, 4}));

  data[{0, 0}] = 4;
  BOOST_TEST((data[{0, 0}]) == 4);
  BOOST_TEST(data.at({0, 0}) == 4);
  data[{0, 0}] = 5;
  BOOST_TEST((data[{0, 0}]) == 5);
  BOOST_TEST(data.at({0, 0}) == 5);

  data[{0, 1}] = 6;
  BOOST_TEST((data[{0, 1}]) == 6);
  BOOST_TEST(data.at({0, 1}) == 6);

  BOOST_TEST((data[{0, 0}]) == 5);

  data[{0, 2}] = 7;
  BOOST_TEST((data[{0, 2}]) == 7);
  BOOST_TEST(data.at({0, 2}) == 7);

  BOOST_TEST((data[{0, 0}]) == 5);
  BOOST_TEST((data[{0, 1}]) == 6);

  data[{1, 0}] = 15;
  BOOST_TEST((data[{1, 0}]) == 15);
  BOOST_TEST(data.at({1, 0}) == 15);

  BOOST_TEST((data[{0, 0}]) == 5);
  BOOST_TEST((data[{0, 1}]) == 6);
  BOOST_TEST((data[{0, 2}]) == 7);

  data[{1, 1}] = 16;
  BOOST_TEST((data[{1, 1}]) == 16);
  BOOST_TEST(data.at({1, 1}) == 16);

  BOOST_TEST((data[{0, 0}]) == 5);
  BOOST_TEST((data[{0, 1}]) == 6);
  BOOST_TEST((data[{0, 2}]) == 7);
  BOOST_TEST((data[{1, 0}]) == 15);

  data[{1, 2}] = 17;
  BOOST_TEST((data[{1, 2}]) == 17);
  BOOST_TEST(data.at({1, 2}) == 17);

  BOOST_TEST((data[{0, 0}]) == 5);
  BOOST_TEST((data[{0, 1}]) == 6);
  BOOST_TEST((data[{0, 2}]) == 7);
  BOOST_TEST((data[{1, 0}]) == 15);
  BOOST_TEST((data[{1, 1}]) == 16);

  BOOST_CHECK_THROW(data.at({0, 3}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 4}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 3}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 4}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 3}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 4}), std::out_of_range);

  BOOST_TEST(data.first() == 5);
  data.first() = -5;
  BOOST_TEST((data[{0, 0}]) == -5);
  BOOST_TEST(data.first() == -5);
  data.first() = 5;

  BOOST_TEST(data.last() == 17);
  data.last() = -17;
  BOOST_TEST((data[{1U, 2U}]) == -17);
  BOOST_TEST(data.last() == -17);
  data.last() = 17;

  auto const& constData = data;

  BOOST_TEST(constData.size() == N);

  static_assert(std::decay_t<decltype(constData)>::dimensions() == 2U);
  BOOST_TEST(constData.dimSize<0U>() == NCryostats);
  BOOST_TEST(constData.dimSize<1U>() == NTPCs);
  BOOST_TEST(constData.dimSize<2U>() == 0U);
  BOOST_TEST(constData.dimSize<3U>() == 0U);

  BOOST_TEST(std::addressof(constData.first()) == std::addressof(data.first()));
  BOOST_TEST(std::addressof(constData.last()) == std::addressof(data.last()));

  BOOST_TEST((constData[{0, 0}]) == (data[{0, 0}]));
  BOOST_TEST((constData[{0, 1}]) == (data[{0, 1}]));
  BOOST_TEST((constData[{0, 2}]) == (data[{0, 2}]));
  BOOST_TEST((constData[{1, 0}]) == (data[{1, 0}]));
  BOOST_TEST((constData[{1, 1}]) == (data[{1, 1}]));
  BOOST_TEST((constData[{1, 2}]) == (data[{1, 2}]));
  BOOST_TEST(constData.at({0, 0}) == data.at({0, 0}));
  BOOST_TEST(constData.at({0, 1}) == data.at({0, 1}));
  BOOST_TEST(constData.at({0, 2}) == data.at({0, 2}));
  BOOST_TEST(constData.at({1, 0}) == data.at({1, 0}));
  BOOST_TEST(constData.at({1, 1}) == data.at({1, 1}));
  BOOST_TEST(constData.at({1, 2}) == data.at({1, 2}));

  BOOST_CHECK_THROW(constData.at({0, 3}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 4}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 3}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 4}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 3}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 4}), std::out_of_range);

  auto const cb = constData.begin();
  auto const ce = constData.end();
  BOOST_TEST(static_cast<size_t>(ce - cb) == N);

  // simple read-only iteration test
  expected_index = 0U;
  for (auto& value : constData) {
    static_assert(
      std::is_same_v<decltype(value), std::decay_t<decltype(constData)>::const_reference>);

    geo::TPCID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);

  // ID/data pair read-only iteration test
  expected_index = 0U;
  for (auto&& [ID, value] : constData.items()) {
    static_assert(std::is_same_v<decltype(ID), geo::TPCID>);
    static_assert(
      std::is_same_v<decltype(value), std::decay_t<decltype(constData)>::const_reference>);

    geo::TPCID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);

  data.fill(14);
  for (auto c : util::counter<unsigned int>(NCryostats))
    for (auto t : util::counter<unsigned int>(NTPCs))
      BOOST_TEST((data[{c, t}]) == 14);

  data.reset();
  for (auto c : util::counter<unsigned int>(NCryostats))
    for (auto t : util::counter<unsigned int>(NTPCs))
      BOOST_TEST((data[{c, t}]) == 0);
} // TPCDataContainerTest()

BOOST_AUTO_TEST_SUITE(geometrydatacontainers_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(TPCDataContainerTestCase)
{
  constexpr std::size_t NCryostats = 2U;
  constexpr std::size_t NTPCs = 3U;

  geo::TPCDataContainer<int> data1(NCryostats, NTPCs);
  TPCDataContainerTest(data1, NCryostats, NTPCs);
}

//------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
