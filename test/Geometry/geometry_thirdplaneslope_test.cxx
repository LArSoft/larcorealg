/**
 * @file   geometry_thirdplaneslope_test.cxx
 * @brief  Simple unit test on a standard detector
 * @date   June 2nd, 2015
 * @author petrillo@fnal.gov
 *
 * Usage: just run the executable.
 * Boost unit testing environment keeps the arguments secret anyway.
 */

// Boost test libraries; we want to define this stuff as soon as possible
#define BOOST_TEST_MODULE GeometryThirdPlaneSlopeTest

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_base.h"
#include "larcorealg/TestUtils/boost_unit_test_base.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"

// utility libraries

// C/C++ standard libraries

using boost::test_tools::tolerance;

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// in the specific, the type of the channel mapping and a proper test name,
// used for output only; BasicGeometryEnvironmentConfiguration can read the
// configuration file name from command line, and
// BoostCommandLineConfiguration<> makes it initialize in time for Boost
// to catch it when instanciating the fixture.
struct StandardGeometryConfiguration:
  public testing::BoostCommandLineConfiguration<
    testing::BasicGeometryEnvironmentConfiguration<geo::ChannelMapStandardAlg>
    >
{
  /// Constructor: overrides the application name
  StandardGeometryConfiguration()
    { SetApplicationName("GeometryThirdPlaneSlopeTest"); }
}; // class StandardGeometryConfiguration

/*
 * Our fixture is based on GeometryTesterFixture, configured with the object
 * above.
 */
using SimpleGeometryTestFixture
  = testing::GeometryTesterEnvironment<StandardGeometryConfiguration>;



//------------------------------------------------------------------------------
//---  The tests
//---
//
// Note on Boost fixture options:
// - BOOST_FIXTURE_TEST_SUITE will provide the tester as a always-present data
//   member in the environment, as "Tester()"; but a new fixture, with a new
//   geometry and a new tester, is initialized on each test case
// - BOOST_GLOBAL_FIXTURE does not provide tester access, so one has to get it
//   as GeometryIteratorTestFixture::GlobalTester(); on the other hand, the
//   fixture is initialized only when a new global one is explicitly created.
//

BOOST_FIXTURE_TEST_SUITE(GeometryIterators, SimpleGeometryTestFixture)
// BOOST_GLOBAL_FIXTURE(SimpleGeometryTestFixture)


BOOST_AUTO_TEST_CASE( AllTests )
{
  geo::GeometryCore const& geom = *(Provider<geo::GeometryCore>());

  const double angle_u = 1. / 3. * util::pi<double>();
  const double angle_v = 2. / 3. * util::pi<double>();
  const double angle_w = 1. / 2. * util::pi<double>();

  BOOST_TEST_MESSAGE(
    "Wire angles: u=" << angle_u << " v=" << angle_v << " => w=" << angle_w
    );

  const double slope_u = 1. / std::sqrt(3);
  const double slope_v = 1. / std::sqrt(3);

  const double expected_slope_w = 0.5;

  double slope_w = geom.ComputeThirdPlaneSlope
    (angle_u, slope_u, angle_v, slope_v, angle_w);

  BOOST_TEST_MESSAGE(
    "Slopes: s(u)=" << slope_u << " s(v)=" << slope_v << " => s(w)=" << slope_w
    );

  BOOST_TEST(slope_w == expected_slope_w, 0.01% tolerance());

} // BOOST_AUTO_TEST_CASE( AllTests )


BOOST_AUTO_TEST_SUITE_END()
