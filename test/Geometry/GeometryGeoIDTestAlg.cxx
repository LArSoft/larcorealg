/**
 * @file   GeometryGeoIDTestAlg.cxx
 * @brief  Unit test for geometry iterators
 * @date   October 31, 2016
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 * The methods require a Boost test environment.
 */

// LArSoft libraries
#include "GeometryGeoIDTestAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Boost libraries
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <type_traits>

//-----------------------------------------------------------------------------
unsigned int geo::GeometryGeoIDTestAlg::Run() const
{
  // All the tests

  CryostatGeoIDTest();
  TPCGeoIDTest();
  PlaneGeoIDTest();
  WireGeoIDTest();

  return 0;
}

//-----------------------------------------------------------------------------
void geo::GeometryGeoIDTestAlg::CryostatGeoIDTest() const
{
  auto iCryo = geom->begin<CryostatID>();
  for (CryostatGeo const& cryostat : geom->Iterate<CryostatGeo>()) {

    CryostatID const& ID = cryostat.ID();

    // the ID of this CryostatGeo is the expected one in a sequential scheme:
    BOOST_TEST(ID == *iCryo);

    // the ID of this CryostatGeo is associated to the CryostatGeo itself
    auto const& cryostatFromID = geom->Cryostat(ID);
    BOOST_TEST(&cryostat == &cryostatFromID);

    ++iCryo;
  } // for cryostat
}

//-----------------------------------------------------------------------------
void geo::GeometryGeoIDTestAlg::TPCGeoIDTest() const
{
  auto iTPC = geom->begin<TPCID>();
  for (TPCGeo const& tpc : geom->Iterate<TPCGeo>()) {

    TPCID const& ID = tpc.ID();

    // the ID of this TPCGeo is the expected one in a sequential scheme:
    BOOST_TEST(ID == *iTPC);

    // the ID of this TPCGeo is associated to the TPCGeo itself
    auto const& TPCFromID = geom->TPC(ID);
    BOOST_TEST(&tpc == &TPCFromID);

    ++iTPC;
  } // for TPC
}

//-----------------------------------------------------------------------------
void geo::GeometryGeoIDTestAlg::PlaneGeoIDTest() const
{
  auto iPlane = wireReadoutGeom->begin<PlaneID>();
  for (PlaneGeo const& plane : wireReadoutGeom->Iterate<PlaneGeo>()) {

    PlaneID const& ID = plane.ID();

    // the ID of this PlaneGeo is the expected one in a sequential scheme:
    BOOST_TEST(ID == *iPlane);

    // the ID of this PlaneGeo is associated to the PlaneGeo itself
    auto const& planeFromID = wireReadoutGeom->Plane(ID);
    BOOST_TEST(&plane == &planeFromID);

    ++iPlane;
  } // for plane
}

//-----------------------------------------------------------------------------
void geo::GeometryGeoIDTestAlg::WireGeoIDTest() const
{
  auto iWire = wireReadoutGeom->begin<WireID>();
  for (WireGeo const& wire [[maybe_unused]] : wireReadoutGeom->Iterate<WireGeo>()) {
    ++iWire;
  }
}

//-----------------------------------------------------------------------------
