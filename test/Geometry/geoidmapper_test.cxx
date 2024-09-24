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
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/Geometry/GeometryIDmapper.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

//------------------------------------------------------------------------------
void TPCIDmappingTest(geo::TPCIDmapper<> mapper, // copy here is intentional
                      std::size_t const NCryostats,
                      std::size_t const NTPCs)
{
  using Mapper_t = geo::TPCIDmapper<>;

  std::size_t const N = NCryostats * NTPCs;

  static_assert(mapper.dimensions() == 2U);
  BOOST_TEST(mapper.dimSize<0U>() == NCryostats);
  BOOST_TEST(mapper.dimSize<1U>() == NTPCs);
  BOOST_TEST(mapper.dimSize<2U>() == 0U);
  BOOST_TEST(mapper.dimSize<3U>() == 0U);

  BOOST_TEST(!mapper.empty());
  BOOST_TEST(mapper.size() == N);

  auto expected_index = Mapper_t::index_type{0};
  for (auto c : util::counter<unsigned int>(NCryostats)) {
    for (auto t : util::counter<unsigned int>(NTPCs)) {
      geo::TPCID const expected_ID{c, t};

      auto const& ID = mapper.ID(expected_index);

      BOOST_TEST(mapper.index(expected_ID) == expected_index);
      BOOST_TEST(ID == expected_ID);
      BOOST_TEST(ID.isValid);

      ++expected_index;
    } // for TPCs
  }   // for cryostats
  BOOST_TEST(mapper.size() == expected_index);
  BOOST_TEST(!mapper.ID(expected_index).isValid);

  BOOST_TEST(mapper.firstID() == geo::TPCID(0, 0));
  BOOST_TEST(mapper.lastID() == geo::TPCID(1, 2));

  BOOST_TEST(mapper.hasElement({0, 0}));
  BOOST_TEST(mapper.hasElement({0, 1}));
  BOOST_TEST(mapper.hasElement({0, 2}));
  BOOST_TEST(!mapper.hasElement({0, 3}));
  BOOST_TEST(!mapper.hasElement({0, 4}));
  BOOST_TEST(mapper.hasElement({1, 0}));
  BOOST_TEST(mapper.hasElement({1, 1}));
  BOOST_TEST(mapper.hasElement({1, 2}));
  BOOST_TEST(!mapper.hasElement({1, 3}));
  BOOST_TEST(!mapper.hasElement({1, 4}));
  BOOST_TEST(!mapper.hasElement({2, 0}));
  BOOST_TEST(!mapper.hasElement({2, 1}));
  BOOST_TEST(!mapper.hasElement({2, 2}));
  BOOST_TEST(!mapper.hasElement({2, 3}));
  BOOST_TEST(!mapper.hasElement({2, 4}));

  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{0, 0}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{0, 1}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{0, 2}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{0, 3}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{0, 4}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{1, 0}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{1, 1}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{1, 2}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{1, 3}));
  BOOST_TEST(mapper.hasElement<geo::CryostatID>(geo::TPCID{1, 4}));
  BOOST_TEST(!mapper.hasElement<geo::CryostatID>(geo::TPCID{2, 0}));
  BOOST_TEST(!mapper.hasElement<geo::CryostatID>(geo::TPCID{2, 1}));
  BOOST_TEST(!mapper.hasElement<geo::CryostatID>(geo::TPCID{2, 2}));
  BOOST_TEST(!mapper.hasElement<geo::CryostatID>(geo::TPCID{2, 3}));
  BOOST_TEST(!mapper.hasElement<geo::CryostatID>(geo::TPCID{2, 4}));

  auto const& constMapper = mapper;

  BOOST_TEST(constMapper.dimSize<0U>() == NCryostats);
  BOOST_TEST(constMapper.dimSize<1U>() == NTPCs);
  BOOST_TEST(constMapper.dimSize<2U>() == 0U);
  BOOST_TEST(constMapper.dimSize<3U>() == 0U);
}

BOOST_AUTO_TEST_SUITE(geoidmapper_test)

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(TPCIDmappingTestCase)
{
  constexpr std::size_t NCryostats = 2U;
  constexpr std::size_t NTPCs = 3U;

  geo::TPCIDmapper<> mapper1(NCryostats, NTPCs);
  TPCIDmappingTest(mapper1, NCryostats, NTPCs);
}

//------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
