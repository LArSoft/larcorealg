/**
 * @file   readoutdatacontainers_test.cxx
 * @brief  Unit test for `ReadoutDataContainers.h` library.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 7, 2019
 */

// Boost libraries
#define BOOST_TEST_MODULE (readout data containers test)
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/Geometry/ReadoutDataContainers.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"


//------------------------------------------------------------------------------
template <typename T>
struct Summer {

  T sum = T { 0 };

  void operator() (T v) { sum += v; }

  T get() const { return sum; }
  void reset() { sum = T{0}; }

}; // struct Summer

//------------------------------------------------------------------------------
void TPCsetDataContainerTest(
  readout::TPCsetDataContainer<int> data, // copy is intentional
  std::size_t const NCryostats,
  std::size_t const NTPCsets
) {

  std::size_t const N = NCryostats * NTPCsets;

  BOOST_TEST(!data.empty());
  BOOST_TEST(data.size() == N);
  BOOST_TEST(data.capacity() >= N);

  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto s: util::counter<unsigned short int>(NTPCsets))
      BOOST_TEST((data[{ c, s }]) == 0);

  BOOST_TEST(data.firstID() == readout::TPCsetID(0, 0));
  BOOST_TEST(data.lastID() == readout::TPCsetID(1, 2));


  std::size_t expected_index = 0U;

  // simple R/W iteration test
  for (auto& value: data) {
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    readout::TPCsetID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);

  // ID/data pair R/W iteration test
  expected_index = 0U;
  for (auto&& [ ID, value ]: data.items()) {
    static_assert(std::is_same_v<decltype(ID), readout::TPCsetID>);
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    readout::TPCsetID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);

  BOOST_TEST( data.hasTPCset({ 0,  0}));
  BOOST_TEST( data.hasTPCset({ 0,  1}));
  BOOST_TEST( data.hasTPCset({ 0,  2}));
  BOOST_TEST(!data.hasTPCset({ 0,  3}));
  BOOST_TEST(!data.hasTPCset({ 0,  4}));
  BOOST_TEST( data.hasTPCset({ 1,  0}));
  BOOST_TEST( data.hasTPCset({ 1,  1}));
  BOOST_TEST( data.hasTPCset({ 1,  2}));
  BOOST_TEST(!data.hasTPCset({ 1,  3}));
  BOOST_TEST(!data.hasTPCset({ 1,  4}));
  BOOST_TEST(!data.hasTPCset({ 2,  0}));
  BOOST_TEST(!data.hasTPCset({ 2,  1}));
  BOOST_TEST(!data.hasTPCset({ 2,  2}));
  BOOST_TEST(!data.hasTPCset({ 2,  3}));
  BOOST_TEST(!data.hasTPCset({ 2,  4}));

  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 0,  0}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 0,  1}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 0,  2}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 0,  3}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 0,  4}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 1,  0}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 1,  1}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 1,  2}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 1,  3}));
  BOOST_TEST( data.hasCryostat(readout::TPCsetID{ 1,  4}));
  BOOST_TEST(!data.hasCryostat(readout::TPCsetID{ 2,  0}));
  BOOST_TEST(!data.hasCryostat(readout::TPCsetID{ 2,  1}));
  BOOST_TEST(!data.hasCryostat(readout::TPCsetID{ 2,  2}));
  BOOST_TEST(!data.hasCryostat(readout::TPCsetID{ 2,  3}));
  BOOST_TEST(!data.hasCryostat(readout::TPCsetID{ 2,  4}));

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

    readout::TPCsetID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);

  // ID/data pair read-only iteration test
  expected_index = 0U;
  for (auto&& [ ID, value ]: constData.items()) {
    static_assert(std::is_same_v<decltype(ID), readout::TPCsetID>);
    static_assert(std::is_same_v
      <decltype(value), std::decay_t<decltype(constData)>::const_reference>
      );

    readout::TPCsetID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);


  data.fill(14);
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto s: util::counter<unsigned short int>(NTPCsets))
      BOOST_TEST((data[{ c, s }]) == 14);

  data.apply([](int& v){ v *= 2; });
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto s: util::counter<unsigned short int>(NTPCsets))
      BOOST_TEST((data[{ c, s }]) == 28);

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
    for (auto s: util::counter<unsigned short int>(NTPCsets))
      BOOST_TEST((data[{ c, s }]) == 0);

  data.clear();
  BOOST_TEST(data.empty());

} // TPCsetDataContainerTest()


//------------------------------------------------------------------------------
void ROPDataContainerTest(
  readout::ROPDataContainer<int> data, // copy is intentional
  std::size_t const NCryostats,
  std::size_t const NTPCsets,
  std::size_t const NROPs
) {

  std::size_t const N = NCryostats * NTPCsets * NROPs;

  BOOST_TEST(!data.empty());
  BOOST_TEST(data.size() == N);
  BOOST_TEST(data.capacity() >= N);

  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto s: util::counter<unsigned short int>(NTPCsets))
      for (auto r: util::counter<unsigned int>(NROPs))
        BOOST_TEST((data[{ c, s, r }]) == 0);

  BOOST_TEST(data.firstID() == readout::ROPID(0, 0, 0));
  BOOST_TEST(data.lastID() == readout::ROPID(1, 2, 1));


  std::size_t expected_index = 0U;

  // simple R/W iteration test
  for (auto& value: data) {
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    readout::ROPID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);

  // ID/data pair R/W iteration test
  expected_index = 0U;
  for (auto&& [ ID, value ]: data.items()) {
    static_assert(std::is_same_v<decltype(ID), readout::ROPID>);
    static_assert(std::is_same_v<decltype(value), decltype(data)::reference>);

    readout::ROPID const expected_ID = data.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == data[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(data.size() == expected_index);


  BOOST_TEST( data.hasROP({ 0, 0, 0}));
  BOOST_TEST( data.hasROP({ 0, 0, 1}));
  BOOST_TEST(!data.hasROP({ 0, 0, 2}));
  BOOST_TEST( data.hasROP({ 0, 1, 0}));
  BOOST_TEST( data.hasROP({ 0, 1, 1}));
  BOOST_TEST(!data.hasROP({ 0, 1, 2}));
  BOOST_TEST( data.hasROP({ 0, 2, 0}));
  BOOST_TEST( data.hasROP({ 0, 2, 1}));
  BOOST_TEST(!data.hasROP({ 0, 2, 2}));
  BOOST_TEST(!data.hasROP({ 0, 3, 0}));
  BOOST_TEST(!data.hasROP({ 0, 3, 1}));
  BOOST_TEST(!data.hasROP({ 0, 3, 2}));
  BOOST_TEST(!data.hasROP({ 0, 4, 0}));
  BOOST_TEST(!data.hasROP({ 0, 4, 1}));
  BOOST_TEST(!data.hasROP({ 0, 4, 2}));
  BOOST_TEST( data.hasROP({ 1, 0, 0}));
  BOOST_TEST( data.hasROP({ 1, 0, 1}));
  BOOST_TEST(!data.hasROP({ 1, 0, 2}));
  BOOST_TEST( data.hasROP({ 1, 1, 0}));
  BOOST_TEST( data.hasROP({ 1, 1, 1}));
  BOOST_TEST(!data.hasROP({ 1, 1, 2}));
  BOOST_TEST( data.hasROP({ 1, 2, 0}));
  BOOST_TEST( data.hasROP({ 1, 2, 1}));
  BOOST_TEST(!data.hasROP({ 1, 2, 2}));
  BOOST_TEST(!data.hasROP({ 1, 3, 0}));
  BOOST_TEST(!data.hasROP({ 1, 3, 1}));
  BOOST_TEST(!data.hasROP({ 1, 3, 2}));
  BOOST_TEST(!data.hasROP({ 1, 4, 0}));
  BOOST_TEST(!data.hasROP({ 1, 4, 1}));
  BOOST_TEST(!data.hasROP({ 1, 4, 2}));
  BOOST_TEST(!data.hasROP({ 2, 0, 0}));
  BOOST_TEST(!data.hasROP({ 2, 0, 1}));
  BOOST_TEST(!data.hasROP({ 2, 0, 2}));
  BOOST_TEST(!data.hasROP({ 2, 1, 0}));
  BOOST_TEST(!data.hasROP({ 2, 1, 1}));
  BOOST_TEST(!data.hasROP({ 2, 1, 2}));
  BOOST_TEST(!data.hasROP({ 2, 2, 0}));
  BOOST_TEST(!data.hasROP({ 2, 2, 1}));
  BOOST_TEST(!data.hasROP({ 2, 2, 2}));
  BOOST_TEST(!data.hasROP({ 2, 3, 0}));
  BOOST_TEST(!data.hasROP({ 2, 3, 1}));
  BOOST_TEST(!data.hasROP({ 2, 3, 2}));
  BOOST_TEST(!data.hasROP({ 2, 4, 0}));
  BOOST_TEST(!data.hasROP({ 2, 4, 1}));
  BOOST_TEST(!data.hasROP({ 2, 4, 2}));

  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 0, 0}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 0, 1}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 0, 2}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 1, 0}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 1, 1}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 1, 2}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 2, 0}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 2, 1}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 0, 2, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 0, 3, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 0, 3, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 0, 3, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 0, 4, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 0, 4, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 0, 4, 2}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 0, 0}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 0, 1}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 0, 2}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 1, 0}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 1, 1}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 1, 2}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 2, 0}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 2, 1}));
  BOOST_TEST( data.hasTPCset(readout::ROPID{ 1, 2, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 1, 3, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 1, 3, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 1, 3, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 1, 4, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 1, 4, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 1, 4, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 0, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 0, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 0, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 1, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 1, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 1, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 2, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 2, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 2, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 3, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 3, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 3, 2}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 4, 0}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 4, 1}));
  BOOST_TEST(!data.hasTPCset(readout::ROPID{ 2, 4, 2}));

  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 0, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 0, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 0, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 1, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 1, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 1, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 2, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 2, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 2, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 3, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 3, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 3, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 4, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 4, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 0, 4, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 0, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 0, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 0, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 1, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 1, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 1, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 2, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 2, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 2, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 3, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 3, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 3, 2}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 4, 0}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 4, 1}));
  BOOST_TEST( data.hasCryostat(readout::ROPID{ 1, 4, 2}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 0, 0}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 0, 1}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 0, 2}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 1, 0}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 1, 1}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 1, 2}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 2, 0}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 2, 1}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 2, 2}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 3, 0}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 3, 1}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 3, 2}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 4, 0}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 4, 1}));
  BOOST_TEST(!data.hasCryostat(readout::ROPID{ 2, 4, 2}));


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

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);
  BOOST_TEST(  (data[{0, 2, 0}]) ==  25);
  BOOST_TEST(  (data[{0, 2, 1}]) ==  26);
  BOOST_TEST(  (data[{1, 0, 0}]) == 105);
  BOOST_TEST(  (data[{1, 0, 1}]) == 106);
  BOOST_TEST(  (data[{1, 1, 0}]) == 115);
  BOOST_TEST(  (data[{1, 1, 1}]) == 116);

  data[{1, 2, 1}] = 126;
  BOOST_TEST(  (data[{1, 2, 1}]) == 126);
  BOOST_TEST(data.at({1, 2, 1}) ==  126);

  BOOST_TEST(  (data[{0, 0, 0}]) ==   5);
  BOOST_TEST(  (data[{0, 0, 1}]) ==   6);
  BOOST_TEST(  (data[{0, 1, 0}]) ==  15);
  BOOST_TEST(  (data[{0, 1, 1}]) ==  16);
  BOOST_TEST(  (data[{0, 2, 0}]) ==  25);
  BOOST_TEST(  (data[{0, 2, 1}]) ==  26);
  BOOST_TEST(  (data[{1, 0, 0}]) == 105);
  BOOST_TEST(  (data[{1, 0, 1}]) == 106);
  BOOST_TEST(  (data[{1, 1, 0}]) == 115);
  BOOST_TEST(  (data[{1, 1, 1}]) == 116);
  BOOST_TEST(  (data[{1, 2, 0}]) == 125);


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

    readout::ROPID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);

  // ID/data pair read-only iteration test
  expected_index = 0U;
  for (auto&& [ ID, value ]: constData.items()) {
    static_assert(std::is_same_v<decltype(ID), readout::ROPID>);
    static_assert(std::is_same_v
      <decltype(value), std::decay_t<decltype(constData)>::const_reference>
      );

    readout::ROPID const expected_ID = constData.mapper().ID(expected_index);
    BOOST_TEST(ID == expected_ID);
    BOOST_TEST(value == constData[expected_ID]);

    ++expected_index;
  } // for
  BOOST_TEST(constData.size() == expected_index);


  data.fill(14);
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto s: util::counter<unsigned short int>(NTPCsets))
      for (auto r: util::counter<unsigned int>(NROPs))
        BOOST_TEST((data[{ c, s, r }]) == 14);

  data.apply([](int& v){ v *= 2; });
  for (auto c: util::counter<unsigned int>(NCryostats))
    for (auto s: util::counter<unsigned short int>(NTPCsets))
      for (auto r: util::counter<unsigned int>(NROPs))
        BOOST_TEST((data[{ c, s, r }]) == 28);

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
    for (auto s: util::counter<unsigned short int>(NTPCsets))
      for (auto r: util::counter<unsigned int>(NROPs))
        BOOST_TEST((data[{ c, s, r }]) == 0);

  data.clear();
  BOOST_TEST(data.empty());


} // ROPDataContainerTest()


BOOST_AUTO_TEST_SUITE(readoutdatacontainers_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(TPCsetDataContainerTestCase) {

  constexpr std::size_t NCryostats = 2U;
  constexpr std::size_t NTPCsets   = 3U;

  //
  // size constructor
  //
  readout::TPCsetDataContainer<int> data1(NCryostats, NTPCsets);
  TPCsetDataContainerTest(data1, NCryostats, NTPCsets);

  //
  // default constructor + resize
  //
  readout::TPCsetDataContainer<int> data2;
  BOOST_TEST(data2.empty());

  data2.resizeAs(data1);
  TPCsetDataContainerTest(data2, NCryostats, NTPCsets);

} // TPCsetDataContainerTestCase


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(ROPDataContainerTestCase) {

  constexpr std::size_t NCryostats = 2U;
  constexpr std::size_t NTPCsets   = 3U;
  constexpr std::size_t NROPs      = 2U;

  //
  // size constructor
  //
  readout::ROPDataContainer<int> data1(NCryostats, NTPCsets, NROPs);
  ROPDataContainerTest(data1, NCryostats, NTPCsets, NROPs);

  //
  // default constructor + resize
  //
  readout::ROPDataContainer<int> data2;
  BOOST_TEST(data2.empty());

  data2.resizeAs(data1);
  ROPDataContainerTest(data2, NCryostats, NTPCsets, NROPs);

} // ROPDataContainerTestCase


//------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
