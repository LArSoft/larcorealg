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
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"
#include "larcorealg/TestUtils/boost_unit_test_base.h"
#include "larcorealg/TestUtils/geometry_unit_test_base.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// art libraries
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory>

using boost::test_tools::tolerance;

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed; in the specific, the type of the
// channel mapping and a proper test name, used for output only;
// BasicGeometryEnvironmentConfiguration can read the configuration file name from command
// line, and BoostCommandLineConfiguration<> makes it initialize in time for Boost to
// catch it when instanciating the fixture.
struct StandardGeometryConfiguration
  : testing::BoostCommandLineConfiguration<testing::BasicGeometryEnvironmentConfiguration> {
  StandardGeometryConfiguration() { SetApplicationName("GeometryThirdPlaneSlopeTest"); }
};

// Our fixture is based on GeometryTesterFixture, configured with the object above.
using SimpleGeometryTestFixture =
  testing::GeometryTesterEnvironment<StandardGeometryConfiguration, geo::GeoObjectSorterStandard>;

//------------------------------------------------------------------------------
//---  The tests
//---

BOOST_FIXTURE_TEST_SUITE(GeometryIterators, SimpleGeometryTestFixture)

BOOST_AUTO_TEST_CASE(AllTests)
{
  auto const wireReadoutGeom = lar::standalone::SetupReadout({}, Provider<geo::GeometryCore>());

  double const angle_u = 1. / 3. * util::pi<double>();
  double const angle_v = 2. / 3. * util::pi<double>();
  double const angle_w = 1. / 2. * util::pi<double>();

  BOOST_TEST_MESSAGE("Wire angles: u=" << angle_u << " v=" << angle_v << " => w=" << angle_w);

  double const slope_u = 1. / std::sqrt(3);
  double const slope_v = 1. / std::sqrt(3);

  double const expected_slope_w = 0.5;

  double const slope_w =
    wireReadoutGeom->ComputeThirdPlaneSlope(angle_u, slope_u, angle_v, slope_v, angle_w);

  BOOST_TEST_MESSAGE("Slopes: s(u)=" << slope_u << " s(v)=" << slope_v << " => s(w)=" << slope_w);

  BOOST_TEST(slope_w == expected_slope_w, 0.01 % tolerance());
}

BOOST_AUTO_TEST_SUITE_END()
