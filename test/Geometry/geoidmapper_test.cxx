/**
 * @file   geoidmapper_test.cxx
 * @brief  Unit test for `geo::GeoIDmapper`.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   October 26, 2019
 * @see    `larcorealg/Geometry/GeometryDataContainers.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE (geo ID mapper test)
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/Geometry/GeometryIDmapper.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"


//------------------------------------------------------------------------------
void TPCIDmappingTest(
  geo::TPCIDmapper<> mapper, // copy here is intentional
  std::size_t const NCryostats, std::size_t const NTPCs
) {

  using Mapper_t = geo::TPCIDmapper<>;

  std::size_t const N = NCryostats * NTPCs;

  static_assert(mapper.dimensions() == 2U);
  BOOST_CHECK_EQUAL(mapper.dimSize<0U>(), NCryostats);
  BOOST_CHECK_EQUAL(mapper.dimSize<1U>(), NTPCs);
  BOOST_CHECK_EQUAL(mapper.dimSize<2U>(), 0U);
  BOOST_CHECK_EQUAL(mapper.dimSize<3U>(), 0U);

  BOOST_CHECK(!mapper.empty());
  BOOST_CHECK_EQUAL(mapper.size(), N);

  auto expected_index = Mapper_t::index_type{ 0 };
  for (auto c: util::counter<unsigned int>(NCryostats)) {
    for (auto t: util::counter<unsigned int>(NTPCs)) {
      geo::TPCID const expected_ID { c, t };

      auto const& ID = mapper.ID(expected_index);

      BOOST_CHECK_EQUAL(mapper.index(expected_ID), expected_index);
      BOOST_CHECK_EQUAL(ID, expected_ID);
      BOOST_CHECK(ID.isValid);

      ++expected_index;
    } // for TPCs
  } // for cryostats
  BOOST_CHECK_EQUAL(mapper.size(), expected_index);
  BOOST_CHECK(!mapper.ID(expected_index).isValid);

  BOOST_CHECK_EQUAL(mapper.firstID(), geo::TPCID(0, 0));
  BOOST_CHECK_EQUAL(mapper.lastID(), geo::TPCID(1, 2));

  BOOST_CHECK( mapper.hasElement({ 0, 0 }));
  BOOST_CHECK( mapper.hasElement({ 0, 1 }));
  BOOST_CHECK( mapper.hasElement({ 0, 2 }));
  BOOST_CHECK(!mapper.hasElement({ 0, 3 }));
  BOOST_CHECK(!mapper.hasElement({ 0, 4 }));
  BOOST_CHECK( mapper.hasElement({ 1, 0 }));
  BOOST_CHECK( mapper.hasElement({ 1, 1 }));
  BOOST_CHECK( mapper.hasElement({ 1, 2 }));
  BOOST_CHECK(!mapper.hasElement({ 1, 3 }));
  BOOST_CHECK(!mapper.hasElement({ 1, 4 }));
  BOOST_CHECK(!mapper.hasElement({ 2, 0 }));
  BOOST_CHECK(!mapper.hasElement({ 2, 1 }));
  BOOST_CHECK(!mapper.hasElement({ 2, 2 }));
  BOOST_CHECK(!mapper.hasElement({ 2, 3 }));
  BOOST_CHECK(!mapper.hasElement({ 2, 4 }));

  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 0, 0 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 0, 1 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 0, 2 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 0, 3 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 0, 4 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 1, 0 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 1, 1 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 1, 2 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 1, 3 }));
  BOOST_CHECK( mapper.hasElement<geo::CryostatID>(geo::TPCID{ 1, 4 }));
  BOOST_CHECK(!mapper.hasElement<geo::CryostatID>(geo::TPCID{ 2, 0 }));
  BOOST_CHECK(!mapper.hasElement<geo::CryostatID>(geo::TPCID{ 2, 1 }));
  BOOST_CHECK(!mapper.hasElement<geo::CryostatID>(geo::TPCID{ 2, 2 }));
  BOOST_CHECK(!mapper.hasElement<geo::CryostatID>(geo::TPCID{ 2, 3 }));
  BOOST_CHECK(!mapper.hasElement<geo::CryostatID>(geo::TPCID{ 2, 4 }));


  auto const& constMapper = mapper;

  BOOST_CHECK_EQUAL(constMapper.dimSize<0U>(), NCryostats);
  BOOST_CHECK_EQUAL(constMapper.dimSize<1U>(), NTPCs);
  BOOST_CHECK_EQUAL(constMapper.dimSize<2U>(), 0U);
  BOOST_CHECK_EQUAL(constMapper.dimSize<3U>(), 0U);


  mapper.clear();
  BOOST_CHECK(mapper.empty());

} // TPCIDmappingTest()


//------------------------------------------------------------------------------
void PlaneIDmappingTest(
  geo::PlaneIDmapper<> mapper, // copy here is intentional
  std::size_t const NCryostats,
  std::size_t const NTPCs,
  std::size_t const NPlanes
) {

  using Mapper_t = geo::PlaneIDmapper<>;

  std::size_t const N = NCryostats * NTPCs * NPlanes;

  static_assert(mapper.dimensions() == 3U);
  BOOST_CHECK_EQUAL(mapper.dimSize<0U>(), NCryostats);
  BOOST_CHECK_EQUAL(mapper.dimSize<1U>(), NTPCs);
  BOOST_CHECK_EQUAL(mapper.dimSize<2U>(), NPlanes);
  BOOST_CHECK_EQUAL(mapper.dimSize<3U>(), 0U);

  BOOST_CHECK(!mapper.empty());
  BOOST_CHECK_EQUAL(mapper.size(), N);

  auto expected_index = Mapper_t::index_type{ 0 };
  for (auto c: util::counter<unsigned int>(NCryostats)) {
    for (auto t: util::counter<unsigned int>(NTPCs)) {
      for (auto p: util::counter<unsigned int>(NPlanes)) {
        geo::PlaneID const expected_ID { c, t, p };

        BOOST_CHECK_EQUAL(mapper.index(expected_ID), expected_index);
        BOOST_CHECK_EQUAL(mapper.ID(expected_index), expected_ID);

        ++expected_index;
      } // for planes
    } // for TPCs
  } // for cryostats
  BOOST_CHECK_EQUAL(mapper.size(), expected_index);

  BOOST_CHECK_EQUAL(mapper.firstID(), geo::PlaneID(0, 0, 0));
  BOOST_CHECK_EQUAL(mapper.lastID(), geo::PlaneID(1, 2, 1));

  BOOST_CHECK( mapper.hasPlane({ 0, 0, 0}));
  BOOST_CHECK( mapper.hasPlane({ 0, 0, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 0, 2}));
  BOOST_CHECK( mapper.hasPlane({ 0, 1, 0}));
  BOOST_CHECK( mapper.hasPlane({ 0, 1, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 1, 2}));
  BOOST_CHECK( mapper.hasPlane({ 0, 2, 0}));
  BOOST_CHECK( mapper.hasPlane({ 0, 2, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 2, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 3, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 3, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 3, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 4, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 4, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 0, 4, 2}));
  BOOST_CHECK( mapper.hasPlane({ 1, 0, 0}));
  BOOST_CHECK( mapper.hasPlane({ 1, 0, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 0, 2}));
  BOOST_CHECK( mapper.hasPlane({ 1, 1, 0}));
  BOOST_CHECK( mapper.hasPlane({ 1, 1, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 1, 2}));
  BOOST_CHECK( mapper.hasPlane({ 1, 2, 0}));
  BOOST_CHECK( mapper.hasPlane({ 1, 2, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 2, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 3, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 3, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 3, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 4, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 4, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 1, 4, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 0, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 0, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 0, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 1, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 1, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 1, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 2, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 2, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 2, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 3, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 3, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 3, 2}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 4, 0}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 4, 1}));
  BOOST_CHECK(!mapper.hasPlane({ 2, 4, 2}));

  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 0, 0}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 0, 1}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 0, 2}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 1, 0}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 1, 1}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 1, 2}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 2, 0}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 2, 1}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 0, 2, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 0, 3, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 0, 3, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 0, 3, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 0, 4, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 0, 4, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 0, 4, 2}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 0, 0}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 0, 1}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 0, 2}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 1, 0}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 1, 1}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 1, 2}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 2, 0}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 2, 1}));
  BOOST_CHECK( mapper.hasTPC(geo::PlaneID{ 1, 2, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 1, 3, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 1, 3, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 1, 3, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 1, 4, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 1, 4, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 1, 4, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 0, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 0, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 0, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 1, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 1, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 1, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 2, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 2, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 2, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 3, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 3, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 3, 2}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 4, 0}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 4, 1}));
  BOOST_CHECK(!mapper.hasTPC(geo::PlaneID{ 2, 4, 2}));

  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 0, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 0, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 0, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 1, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 1, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 1, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 2, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 2, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 2, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 3, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 3, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 3, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 4, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 4, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 0, 4, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 0, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 0, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 0, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 1, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 1, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 1, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 2, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 2, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 2, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 3, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 3, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 3, 2}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 4, 0}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 4, 1}));
  BOOST_CHECK( mapper.hasCryostat(geo::PlaneID{ 1, 4, 2}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 0, 0}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 0, 1}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 0, 2}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 1, 0}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 1, 1}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 1, 2}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 2, 0}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 2, 1}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 2, 2}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 3, 0}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 3, 1}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 3, 2}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 4, 0}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 4, 1}));
  BOOST_CHECK(!mapper.hasCryostat(geo::PlaneID{ 2, 4, 2}));


  auto const& constMapper = mapper;

  BOOST_CHECK_EQUAL(constMapper.dimSize<0U>(), NCryostats);
  BOOST_CHECK_EQUAL(constMapper.dimSize<1U>(), NTPCs);
  BOOST_CHECK_EQUAL(constMapper.dimSize<2U>(), NPlanes);
  BOOST_CHECK_EQUAL(constMapper.dimSize<3U>(), 0U);

  mapper.clear();
  BOOST_CHECK(mapper.empty());

} // PlaneIDmappingTest()


BOOST_AUTO_TEST_SUITE(geoidmapper_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(TPCIDmappingTestCase) {

  constexpr std::size_t NCryostats = 2U;
  constexpr std::size_t NTPCs      = 3U;

  //
  // size constructor
  //
  geo::TPCIDmapper<> mapper1(NCryostats, NTPCs);
  TPCIDmappingTest(mapper1, NCryostats, NTPCs);

  //
  // default constructor + resize
  //
  geo::TPCIDmapper<> mapper2;
  BOOST_CHECK(mapper2.empty());

  mapper2.resizeAs(mapper1);
  TPCIDmappingTest(mapper2, NCryostats, NTPCs);

} // TPCIDmappingTestCase


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(PlaneIDmappingTestCase) {

  constexpr std::size_t NCryostats = 2U;
  constexpr std::size_t NTPCs      = 3U;
  constexpr std::size_t NPlanes    = 2U;

  //
  // size constructor
  //
  geo::PlaneIDmapper<> mapper1(NCryostats, NTPCs, NPlanes);
  PlaneIDmappingTest(mapper1, NCryostats, NTPCs, NPlanes);

  //
  // default constructor + resize
  //
  geo::PlaneIDmapper<> mapper2;
  BOOST_CHECK(mapper2.empty());

  mapper2.resizeAs(mapper1);
  PlaneIDmappingTest(mapper2, NCryostats, NTPCs, NPlanes);

} // PlaneIDmappingTestCase


//------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
