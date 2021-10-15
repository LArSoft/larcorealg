/**
 * @file   GeometryGeoIDTestAlg.cxx
 * @brief  Unit test for geometry iterators
 * @date   October 31, 2016
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 * The methods require a Boost test environment.
 */

// LArSoft libraries
#include "test/Geometry/GeometryGeoIDTestAlg.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/GeometryCore.h"

// Boost libraries
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <type_traits>



//-----------------------------------------------------------------------------
unsigned int geo::GeometryGeoIDTestAlg::Run() const {
  // All the tests

  CryostatGeoIDTest();
  TPCGeoIDTest();
  PlaneGeoIDTest();
  WireGeoIDTest();

  return 0;
} // GeometryGeoIDTestAlg::Run()



//-----------------------------------------------------------------------------
void geo::GeometryGeoIDTestAlg::CryostatGeoIDTest() const {

  auto iCryo = geom->begin_cryostat_id();
  for (geo::CryostatGeo const& cryostat: geom->IterateCryostats()) {

    geo::CryostatID const& ID = cryostat.ID();

    // the ID of this CryostatGeo is the expected one in a sequential scheme:
    BOOST_CHECK_EQUAL(ID, *iCryo);

    // the ID of this CryostatGeo is associated to the CryostatGeo itself
    auto const& cryostatFromID = geom->Cryostat(ID);
    BOOST_CHECK_EQUAL(&cryostat, &cryostatFromID);


    ++iCryo;
  } // for cryostat

} // GeometryGeoIDTestAlg::CryostatGeoIDTest()



//-----------------------------------------------------------------------------
void geo::GeometryGeoIDTestAlg::TPCGeoIDTest() const {

  auto iTPC = geom->begin_TPC_id();
  for (geo::TPCGeo const& tpc: geom->IterateTPCs()) {

    geo::TPCID const& ID = tpc.ID();

    // the ID of this TPCGeo is the expected one in a sequential scheme:
    BOOST_CHECK_EQUAL(ID, *iTPC);

    // the ID of this TPCGeo is associated to the TPCGeo itself
    auto const& TPCFromID = geom->TPC(ID);
    BOOST_CHECK_EQUAL(&tpc, &TPCFromID);


    ++iTPC;
  } // for TPC

} // GeometryGeoIDTestAlg::TPCGeoIDTest()



//-----------------------------------------------------------------------------
void geo::GeometryGeoIDTestAlg::PlaneGeoIDTest() const {

  auto iPlane = geom->begin_plane_id();
  for (geo::PlaneGeo const& plane: geom->IteratePlanes()) {

    geo::PlaneID const& ID = plane.ID();

    // the ID of this PlaneGeo is the expected one in a sequential scheme:
    BOOST_CHECK_EQUAL(ID, *iPlane);

    // the ID of this PlaneGeo is associated to the PlaneGeo itself
    auto const& planeFromID = geom->Plane(ID);
    BOOST_CHECK_EQUAL(&plane, &planeFromID);


    ++iPlane;
  } // for plane

} // GeometryGeoIDTestAlg::PlaneGeoIDTest()



//-----------------------------------------------------------------------------
void geo::GeometryGeoIDTestAlg::WireGeoIDTest() const {

  auto iWire = geom->begin_wire_id();
  for (geo::WireGeo const& wire [[gnu::unused]]: geom->IterateWires()) {

    // this test is disabled since geo::WireID does not support ID()
    // (as of LArSoft 6.13)
    /*
    geo::WireID const& ID = wire.ID();

    // the ID of this WireGeo is the expected one in a sequential scheme:
    BOOST_CHECK_EQUAL(ID, *iWire);

    // the ID of this WireGeo is associated to the WireGeo itself
    auto const& wireFromID = geom->Wire(ID);
    BOOST_CHECK_EQUAL(&wire, &wireFromID);
    */

    ++iWire;
  } // for wire

} // GeometryGeoIDTestAlg::WireGeoIDTest()


//-----------------------------------------------------------------------------
