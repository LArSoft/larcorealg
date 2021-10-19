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
#include "larcorealg/Geometry/GeometryDataContainers.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::static_assert_on()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"


//------------------------------------------------------------------------------
template <typename T>
struct Summer {

  T sum = T { 0 };

  void operator() (T v) { sum += v; }

  T get() const { return sum; }
  void reset() { sum = T{0}; }

}; // struct Summer

//------------------------------------------------------------------------------
void TPCDataContainerTest(
  geo::TPCDataContainer<int> data, // copy here is intentional
  std::size_t const NCryostats, std::size_t const NTPCs
) {

  std::size_t const N = NCryostats * NTPCs;

  static_assert(data.dimensions() == 2U);
  BOOST_TEST(data.dimSize<0U>() == NCryostats);
  BOOST_TEST(data.dimSize<1U>() == NTPCs);
  BOOST_TEST(data.dimSize<2U>() == 0U);
  BOOST_TEST(data.dimSize<3U>() == 0U);

  BOOST_TEST(!data.empty());
  BOOST_TEST(data.size() == N);
  BOOST_TEST(data.capacity() >= N);

  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto t: util::counter<unsigned int>(NTPCs))
      BOOST_TEST((data[{ c, t }]) == 0);

  BOOST_TEST(data.firstID() == geo::TPCID(0, 0));
  BOOST_TEST(data.lastID() == geo::TPCID(1, 2));


  std::size_t expected_index = 0U;

  // simple R/W iteration test
  for (auto& value: data) {
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    geo::TPCID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);

  // ID/data pair R/W iteration test
  expected_index = 0U;
  for (auto&& [ ID, value ]: data.items()) {
    static_assert(std::is_same_v<decltype(ID), geo::TPCID>);
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    geo::TPCID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);


  BOOST_TEST( data.hasTPC({ 0,  0}));
  BOOST_TEST( data.hasTPC({ 0,  1}));
  BOOST_TEST( data.hasTPC({ 0,  2}));
  BOOST_TEST(!data.hasTPC({ 0,  3}));
  BOOST_TEST(!data.hasTPC({ 0,  4}));
  BOOST_TEST( data.hasTPC({ 1,  0}));
  BOOST_TEST( data.hasTPC({ 1,  1}));
  BOOST_TEST( data.hasTPC({ 1,  2}));
  BOOST_TEST(!data.hasTPC({ 1,  3}));
  BOOST_TEST(!data.hasTPC({ 1,  4}));
  BOOST_TEST(!data.hasTPC({ 2,  0}));
  BOOST_TEST(!data.hasTPC({ 2,  1}));
  BOOST_TEST(!data.hasTPC({ 2,  2}));
  BOOST_TEST(!data.hasTPC({ 2,  3}));
  BOOST_TEST(!data.hasTPC({ 2,  4}));

  BOOST_TEST( data.hasCryostat(geo::TPCID{ 0,  0}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 0,  1}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 0,  2}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 0,  3}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 0,  4}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 1,  0}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 1,  1}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 1,  2}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 1,  3}));
  BOOST_TEST( data.hasCryostat(geo::TPCID{ 1,  4}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{ 2,  0}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{ 2,  1}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{ 2,  2}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{ 2,  3}));
  BOOST_TEST(!data.hasCryostat(geo::TPCID{ 2,  4}));

  data[{0, 0}] = 4;
  BOOST_TEST((data[{0, 0}]) == 4);
  BOOST_TEST(data.at({0, 0}) == 4);
  data[{0, 0}] = 5;
  BOOST_TEST((data[{0, 0}]) == 5);
  BOOST_TEST(data.at({0, 0}) == 5);

  data[{0, 1}] = 6;
  BOOST_TEST((data[{0, 1}]) == 6);
  BOOST_TEST(data.at({0, 1}) == 6);

  BOOST_TEST((data[{0, 0}]) ==  5);

  data[{0, 2}] = 7;
  BOOST_TEST((data[{0, 2}]) == 7);
  BOOST_TEST(data.at({0, 2}) == 7);

  BOOST_TEST((data[{0, 0}]) ==  5);
  BOOST_TEST((data[{0, 1}]) ==  6);

  data[{1, 0}] = 15;
  BOOST_TEST((data[{1, 0}]) == 15);
  BOOST_TEST(data.at({1, 0}) == 15);

  BOOST_TEST((data[{0, 0}]) ==  5);
  BOOST_TEST((data[{0, 1}]) ==  6);
  BOOST_TEST((data[{0, 2}]) ==  7);

  data[{1, 1}] = 16;
  BOOST_TEST((data[{1, 1}]) == 16);
  BOOST_TEST(data.at({1, 1}) == 16);

  BOOST_TEST((data[{0, 0}]) ==  5);
  BOOST_TEST((data[{0, 1}]) ==  6);
  BOOST_TEST((data[{0, 2}]) ==  7);
  BOOST_TEST((data[{1, 0}]) == 15);

  data[{1, 2}] = 17;
  BOOST_TEST((data[{1, 2}]) == 17);
  BOOST_TEST(data.at({1, 2}) == 17);

  BOOST_TEST((data[{0, 0}]) ==  5);
  BOOST_TEST((data[{0, 1}]) ==  6);
  BOOST_TEST((data[{0, 2}]) ==  7);
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
  data.first() =  5;

  BOOST_TEST(data.last() == 17);
  data.last() = -17;
  BOOST_TEST((data[{1U, 2U}]) == -17);
  BOOST_TEST(data.last() == -17);
  data.last() =  17;

  auto const& constData = data;

  BOOST_TEST(constData.size() == N);

  static_assert(std::decay_t<decltype(constData)>::dimensions() == 2U);
  BOOST_TEST(constData.dimSize<0U>() == NCryostats);
  BOOST_TEST(constData.dimSize<1U>() == NTPCs);
  BOOST_TEST(constData.dimSize<2U>() == 0U);
  BOOST_TEST(constData.dimSize<3U>() == 0U);

  BOOST_TEST
    (std::addressof(constData.first()) == std::addressof(data.first()));
  BOOST_TEST
    (std::addressof(constData.last()) == std::addressof(data.last()));

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
  for (auto& value: constData) {
    static_assert(std::is_same_v
      <decltype(value), std::decay_t<decltype(constData)>::const_reference>
      );

    geo::TPCID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);

  // ID/data pair read-only iteration test
  expected_index = 0U;
  for (auto&& [ ID, value ]: constData.items()) {
    static_assert(std::is_same_v<decltype(ID), geo::TPCID>);
    static_assert(std::is_same_v
      <decltype(value), std::decay_t<decltype(constData)>::const_reference>
      );

    geo::TPCID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);


  data.fill(14);
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto t: util::counter<unsigned int>(NTPCs))
      BOOST_TEST((data[{ c, t }]) == 14);

  data.apply([](int& v){ v *= 2; });
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto t: util::counter<unsigned int>(NTPCs))
      BOOST_TEST((data[{ c, t }]) == 28);

  Summer<int> summer;
  static_assert(std::is_same_v<decltype(data.apply(summer)), Summer<int>&>);
  data.apply(summer);
  BOOST_TEST(summer.get() == N * 28);

  summer.reset();
  static_assert
    (std::is_same_v<decltype(constData.apply(summer)), Summer<int>&>);
  constData.apply(summer);
  BOOST_TEST(summer.get() == N * 28);

  auto summer1 = data.apply(Summer<int>{});
  BOOST_TEST(summer1.get() == N * 28);

  auto summer2 = constData.apply(Summer<int>{});
  BOOST_TEST(summer2.get() == N * 28);

  data.reset();
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto t: util::counter<unsigned int>(NTPCs))
      BOOST_TEST((data[{ c, t }]) == 0);

  data.clear();
  BOOST_TEST(data.empty());

} // TPCDataContainerTest()


//------------------------------------------------------------------------------
void PlaneDataContainerTest(
  geo::PlaneDataContainer<int> data, // copy here is intentional
  std::size_t const NCryostats,
  std::size_t const NTPCs,
  std::size_t const NPlanes
) {

  std::size_t const N = NCryostats * NTPCs * NPlanes;

  static_assert(data.dimensions() == 3U);
  BOOST_TEST(data.dimSize<0U>() == NCryostats);
  BOOST_TEST(data.dimSize<1U>() == NTPCs);
  BOOST_TEST(data.dimSize<2U>() == NPlanes);
  BOOST_TEST(data.dimSize<3U>() == 0U);

  BOOST_TEST(!data.empty());
  BOOST_TEST(data.size() == N);
  BOOST_TEST(data.capacity() >= N);

  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto t: util::counter<unsigned int>(NTPCs))
      for (auto p: util::counter<unsigned int>(NPlanes))
        BOOST_TEST((data[{ c, t, p }]) == 0);

  BOOST_TEST(data.firstID() == geo::PlaneID(0, 0, 0));
  BOOST_TEST(data.lastID() == geo::PlaneID(1, 2, 1));


  std::size_t expected_index = 0U;

  // simple R/W iteration test
  for (auto& value: data) {
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    geo::PlaneID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);

  // ID/data pair R/W iteration test
  expected_index = 0U;
  for (auto&& [ ID, value ]: data.items()) {
    static_assert(std::is_same_v<decltype(ID), geo::PlaneID>);
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    geo::PlaneID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);

  BOOST_TEST( data.hasPlane({ 0, 0, 0}));
  BOOST_TEST( data.hasPlane({ 0, 0, 1}));
  BOOST_TEST(!data.hasPlane({ 0, 0, 2}));
  BOOST_TEST( data.hasPlane({ 0, 1, 0}));
  BOOST_TEST( data.hasPlane({ 0, 1, 1}));
  BOOST_TEST(!data.hasPlane({ 0, 1, 2}));
  BOOST_TEST( data.hasPlane({ 0, 2, 0}));
  BOOST_TEST( data.hasPlane({ 0, 2, 1}));
  BOOST_TEST(!data.hasPlane({ 0, 2, 2}));
  BOOST_TEST(!data.hasPlane({ 0, 3, 0}));
  BOOST_TEST(!data.hasPlane({ 0, 3, 1}));
  BOOST_TEST(!data.hasPlane({ 0, 3, 2}));
  BOOST_TEST(!data.hasPlane({ 0, 4, 0}));
  BOOST_TEST(!data.hasPlane({ 0, 4, 1}));
  BOOST_TEST(!data.hasPlane({ 0, 4, 2}));
  BOOST_TEST( data.hasPlane({ 1, 0, 0}));
  BOOST_TEST( data.hasPlane({ 1, 0, 1}));
  BOOST_TEST(!data.hasPlane({ 1, 0, 2}));
  BOOST_TEST( data.hasPlane({ 1, 1, 0}));
  BOOST_TEST( data.hasPlane({ 1, 1, 1}));
  BOOST_TEST(!data.hasPlane({ 1, 1, 2}));
  BOOST_TEST( data.hasPlane({ 1, 2, 0}));
  BOOST_TEST( data.hasPlane({ 1, 2, 1}));
  BOOST_TEST(!data.hasPlane({ 1, 2, 2}));
  BOOST_TEST(!data.hasPlane({ 1, 3, 0}));
  BOOST_TEST(!data.hasPlane({ 1, 3, 1}));
  BOOST_TEST(!data.hasPlane({ 1, 3, 2}));
  BOOST_TEST(!data.hasPlane({ 1, 4, 0}));
  BOOST_TEST(!data.hasPlane({ 1, 4, 1}));
  BOOST_TEST(!data.hasPlane({ 1, 4, 2}));
  BOOST_TEST(!data.hasPlane({ 2, 0, 0}));
  BOOST_TEST(!data.hasPlane({ 2, 0, 1}));
  BOOST_TEST(!data.hasPlane({ 2, 0, 2}));
  BOOST_TEST(!data.hasPlane({ 2, 1, 0}));
  BOOST_TEST(!data.hasPlane({ 2, 1, 1}));
  BOOST_TEST(!data.hasPlane({ 2, 1, 2}));
  BOOST_TEST(!data.hasPlane({ 2, 2, 0}));
  BOOST_TEST(!data.hasPlane({ 2, 2, 1}));
  BOOST_TEST(!data.hasPlane({ 2, 2, 2}));
  BOOST_TEST(!data.hasPlane({ 2, 3, 0}));
  BOOST_TEST(!data.hasPlane({ 2, 3, 1}));
  BOOST_TEST(!data.hasPlane({ 2, 3, 2}));
  BOOST_TEST(!data.hasPlane({ 2, 4, 0}));
  BOOST_TEST(!data.hasPlane({ 2, 4, 1}));
  BOOST_TEST(!data.hasPlane({ 2, 4, 2}));

  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 0, 0}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 0, 1}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 0, 2}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 1, 0}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 1, 1}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 1, 2}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 2, 0}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 2, 1}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 0, 2, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 0, 3, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 0, 3, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 0, 3, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 0, 4, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 0, 4, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 0, 4, 2}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 0, 0}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 0, 1}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 0, 2}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 1, 0}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 1, 1}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 1, 2}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 2, 0}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 2, 1}));
  BOOST_TEST( data.hasTPC(geo::PlaneID{ 1, 2, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 1, 3, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 1, 3, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 1, 3, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 1, 4, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 1, 4, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 1, 4, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 0, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 0, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 0, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 1, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 1, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 1, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 2, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 2, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 2, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 3, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 3, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 3, 2}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 4, 0}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 4, 1}));
  BOOST_TEST(!data.hasTPC(geo::PlaneID{ 2, 4, 2}));

  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 0, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 0, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 0, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 1, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 1, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 1, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 2, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 2, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 2, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 3, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 3, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 3, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 4, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 4, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 0, 4, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 0, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 0, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 0, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 1, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 1, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 1, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 2, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 2, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 2, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 3, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 3, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 3, 2}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 4, 0}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 4, 1}));
  BOOST_TEST( data.hasCryostat(geo::PlaneID{ 1, 4, 2}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 0, 0}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 0, 1}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 0, 2}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 1, 0}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 1, 1}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 1, 2}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 2, 0}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 2, 1}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 2, 2}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 3, 0}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 3, 1}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 3, 2}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 4, 0}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 4, 1}));
  BOOST_TEST(!data.hasCryostat(geo::PlaneID{ 2, 4, 2}));


  data[{0, 0, 0}] = 4;
  BOOST_TEST(  (data[{0, 0, 0}]) ==   4);
  BOOST_TEST(data.at({0, 0, 0}) ==    4);
  data[{0, 0, 0}] = 5;
  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(data.at({0, 0, 0}) ==    5);

  data[{0, 0, 1}] = 6;
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(data.at({0, 0, 1}) ==    6);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);

  data[{0, 1, 0}] = 15;
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(data.at({0, 1, 0}) ==   15);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);

  data[{0, 1, 1}] = 16;
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);
  BOOST_TEST(data.at({0, 1, 1}) ==   16);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);

  data[{0, 2, 0}] = 25;
  BOOST_TEST(  (data[{0, 2, 0}]) ==  25);
  BOOST_TEST(data.at({0, 2, 0}) ==   25);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);

  data[{0, 2, 1}] = 26;
  BOOST_TEST(  (data[{0, 2, 1}]) ==  26);
  BOOST_TEST(data.at({0, 2, 1}) ==   26);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);
  BOOST_TEST(  (data[{0, 2, 0}]) ==  25);

  data[{1, 0, 0}] = 105;
  BOOST_TEST(  (data[{1, 0, 0}]) == 105);
  BOOST_TEST(data.at({1, 0, 0}) ==  105);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);
  BOOST_TEST(  (data[{0, 2, 0}]) ==  25);
  BOOST_TEST(  (data[{0, 2, 1}]) ==  26);

  data[{1, 0, 1}] = 106;
  BOOST_TEST(  (data[{1, 0, 1}]) == 106);
  BOOST_TEST(data.at({1, 0, 1}) ==  106);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);
  BOOST_TEST(  (data[{0, 2, 0}]) ==  25);
  BOOST_TEST(  (data[{0, 2, 1}]) ==  26);
  BOOST_TEST(  (data[{1, 0, 0}]) == 105);

  data[{1, 1, 0}] = 115;
  BOOST_TEST(  (data[{1, 1, 0}]) == 115);
  BOOST_TEST(data.at({1, 1, 0}) ==  115);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);
  BOOST_TEST(  (data[{0, 2, 0}]) ==  25);
  BOOST_TEST(  (data[{0, 2, 1}]) ==  26);
  BOOST_TEST(  (data[{1, 0, 0}]) == 105);
  BOOST_TEST(  (data[{1, 0, 1}]) == 106);

  data[{1, 1, 1}] = 116;
  BOOST_TEST(  (data[{1, 1, 1}]) == 116);
  BOOST_TEST(data.at({1, 1, 1}) ==  116);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);
  BOOST_TEST(  (data[{0, 2, 0}]) ==  25);
  BOOST_TEST(  (data[{0, 2, 1}]) ==  26);
  BOOST_TEST(  (data[{1, 0, 0}]) == 105);
  BOOST_TEST(  (data[{1, 0, 1}]) == 106);
  BOOST_TEST(  (data[{1, 1, 0}]) == 115);

  data[{1, 2, 0}] = 125;
  BOOST_TEST(  (data[{1, 2, 0}]) == 125);
  BOOST_TEST(data.at({1, 2, 0}) ==  125);

  BOOST_TEST((data[{0, 0, 0}]) ==   5);
  BOOST_TEST((data[{0, 0, 1}]) ==   6);
  BOOST_TEST((data[{0, 1, 0}]) ==  15);
  BOOST_TEST((data[{0, 1, 1}]) ==  16);
  BOOST_TEST((data[{0, 2, 0}]) ==  25);
  BOOST_TEST((data[{0, 2, 1}]) ==  26);
  BOOST_TEST((data[{1, 0, 0}]) == 105);
  BOOST_TEST((data[{1, 0, 1}]) == 106);
  BOOST_TEST((data[{1, 1, 0}]) == 115);
  BOOST_TEST((data[{1, 1, 1}]) == 116);

  data[{1, 2, 1}] = 126;
  BOOST_TEST(  (data[{1, 2, 1}]) == 126);
  BOOST_TEST(data.at({1, 2, 1}) ==  126);

  BOOST_TEST((data[{0, 0, 0}]) ==   5);
  BOOST_TEST((data[{0, 0, 1}]) ==   6);
  BOOST_TEST((data[{0, 1, 0}]) ==  15);
  BOOST_TEST((data[{0, 1, 1}]) ==  16);
  BOOST_TEST((data[{0, 2, 0}]) ==  25);
  BOOST_TEST((data[{0, 2, 1}]) ==  26);
  BOOST_TEST((data[{1, 0, 0}]) == 105);
  BOOST_TEST((data[{1, 0, 1}]) == 106);
  BOOST_TEST((data[{1, 1, 0}]) == 115);
  BOOST_TEST((data[{1, 1, 1}]) == 116);
  BOOST_TEST((data[{1, 2, 0}]) == 125);


  BOOST_CHECK_THROW(data.at({0, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 0, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 1, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 0, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 1, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 2, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 3, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 3, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 3, 2}), std::out_of_range);

  BOOST_TEST(data.first() == 5);
  data.first() = -5;
  BOOST_TEST((data[{0, 0, 0}]) == -5);
  BOOST_TEST(data.first() == -5);
  data.first() =  5;

  BOOST_TEST(data.last() == 126);
  data.last() = -126;
  BOOST_TEST((data[{1U, 2U, 1U}]) == -126);
  BOOST_TEST(data.last() == -126);
  data.last() =  126;

  auto const& constData = data;

  static_assert(data.dimensions() == 3U);
  BOOST_TEST(data.dimSize<0U>() == NCryostats);
  BOOST_TEST(data.dimSize<1U>() == NTPCs);
  BOOST_TEST(data.dimSize<2U>() == NPlanes);
  BOOST_TEST(data.dimSize<3U>() == 0U);

  BOOST_TEST
    (std::addressof(constData.first()) == std::addressof(data.first()));
  BOOST_TEST
    (std::addressof(constData.last()) == std::addressof(data.last()));

  BOOST_TEST((constData[{0, 0, 0}]) == (data[{0, 0, 0}]));
  BOOST_TEST((constData[{0, 0, 1}]) == (data[{0, 0, 1}]));
  BOOST_TEST((constData[{0, 1, 0}]) == (data[{0, 1, 0}]));
  BOOST_TEST((constData[{0, 1, 1}]) == (data[{0, 1, 1}]));
  BOOST_TEST((constData[{0, 2, 0}]) == (data[{0, 2, 0}]));
  BOOST_TEST((constData[{0, 2, 1}]) == (data[{0, 2, 1}]));
  BOOST_TEST((constData[{1, 0, 0}]) == (data[{1, 0, 0}]));
  BOOST_TEST((constData[{1, 0, 1}]) == (data[{1, 0, 1}]));
  BOOST_TEST((constData[{1, 1, 0}]) == (data[{1, 1, 0}]));
  BOOST_TEST((constData[{1, 1, 1}]) == (data[{1, 1, 1}]));
  BOOST_TEST((constData[{1, 2, 0}]) == (data[{1, 2, 0}]));
  BOOST_TEST((constData[{1, 2, 1}]) == (data[{1, 2, 1}]));
  BOOST_TEST(constData.at({0, 0, 0}) == data.at({0, 0, 0}));
  BOOST_TEST(constData.at({0, 0, 1}) == data.at({0, 0, 1}));
  BOOST_TEST(constData.at({0, 1, 0}) == data.at({0, 1, 0}));
  BOOST_TEST(constData.at({0, 1, 1}) == data.at({0, 1, 1}));
  BOOST_TEST(constData.at({0, 2, 0}) == data.at({0, 2, 0}));
  BOOST_TEST(constData.at({0, 2, 1}) == data.at({0, 2, 1}));
  BOOST_TEST(constData.at({1, 0, 0}) == data.at({1, 0, 0}));
  BOOST_TEST(constData.at({1, 0, 1}) == data.at({1, 0, 1}));
  BOOST_TEST(constData.at({1, 1, 0}) == data.at({1, 1, 0}));
  BOOST_TEST(constData.at({1, 1, 1}) == data.at({1, 1, 1}));
  BOOST_TEST(constData.at({1, 2, 0}) == data.at({1, 2, 0}));
  BOOST_TEST(constData.at({1, 2, 1}) == data.at({1, 2, 1}));

  BOOST_CHECK_THROW(constData.at({0, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 0, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 1, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 0, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 1, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 2, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 3, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 3, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 3, 2}), std::out_of_range);


  auto const cb = constData.begin();
  auto const ce = constData.end();
  BOOST_TEST(static_cast<size_t>(ce - cb) == N);

  // simple read-only iteration test
  expected_index = 0U;
  for (auto& value: constData) {
    static_assert(std::is_same_v
      <decltype(value), std::decay_t<decltype(constData)>::const_reference>
      );

    geo::PlaneID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);

  // ID/data pair read-only iteration test
  expected_index = 0U;
  for (auto&& [ ID, value ]: constData.items()) {
    static_assert(std::is_same_v<decltype(ID), geo::PlaneID>);
    static_assert(std::is_same_v
      <decltype(value), std::decay_t<decltype(constData)>::const_reference>
      );

    geo::PlaneID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);


  data.fill(14);
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto t: util::counter<unsigned int>(NTPCs))
      for (auto p: util::counter<unsigned int>(NPlanes))
        BOOST_TEST((data[{ c, t, p }]) == 14);

  data.apply([](int& v){ v *= 2; });
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto t: util::counter<unsigned int>(NTPCs))
      for (auto p: util::counter<unsigned int>(NPlanes))
        BOOST_TEST((data[{ c, t, p }]) == 28);

  Summer<int> summer;
  static_assert(std::is_same_v<decltype(data.apply(summer)), Summer<int>&>);
  data.apply(summer);
  BOOST_TEST(summer.get() == N * 28);

  summer.reset();
  static_assert
    (std::is_same_v<decltype(constData.apply(summer)), Summer<int>&>);
  constData.apply(summer);
  BOOST_TEST(summer.get() == N * 28);

  auto summer1 = data.apply(Summer<int>{});
  BOOST_TEST(summer1.get() == N * 28);

  auto summer2 = constData.apply(Summer<int>{});
  BOOST_TEST(summer2.get() == N * 28);

  data.reset();
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto t: util::counter<unsigned int>(NTPCs))
      for (auto p: util::counter<unsigned int>(NPlanes))
        BOOST_TEST((data[{ c, t, p }]) == 0);

  data.clear();
  BOOST_TEST(data.empty());

} // PlaneDataContainerTest()


BOOST_AUTO_TEST_SUITE(geometrydatacontainers_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(TPCDataContainerTestCase) {

  constexpr std::size_t NCryostats = 2U;
  constexpr std::size_t NTPCs      = 3U;

  //
  // size constructor
  //
  geo::TPCDataContainer<int> data1(NCryostats, NTPCs);
  TPCDataContainerTest(data1, NCryostats, NTPCs);

  //
  // default constructor + resize
  //
  geo::TPCDataContainer<int> data2;
  BOOST_TEST(data2.empty());

  data2.resizeAs(data1);
  TPCDataContainerTest(data2, NCryostats, NTPCs);

} // TPCDataContainerTestCase

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(PlaneDataContainerTestCase) {

  constexpr std::size_t NCryostats = 2U;
  constexpr std::size_t NTPCs      = 3U;
  constexpr std::size_t NPlanes    = 2U;

  //
  // size constructor
  //
  geo::PlaneDataContainer<int> data1(NCryostats, NTPCs, NPlanes);
  PlaneDataContainerTest(data1, NCryostats, NTPCs, NPlanes);

  //
  // default constructor + resize
  //
  geo::PlaneDataContainer<int> data2;
  BOOST_TEST(data2.empty());

  data2.resizeAs(data1);
  PlaneDataContainerTest(data2, NCryostats, NTPCs, NPlanes);

} // PlaneDataContainerTestCase

//------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
