/**
 * @file   geometry_standardchannelmapping_test.cxx
 * @brief  Unit test of channel mapping on a standard detector
 * @date   June 26th, 2015
 * @author petrillo@fnal.gov
 *
 * Usage: just run the executable.
 * Or plug a FHiCL file in the command line.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate a main()
// function; Boost is pulled in by boost_unit_test_base.h
#define BOOST_TEST_MODULE GeometryStandardChannelMappingTest

// LArSoft libraries
#include "ChannelMapStandardTestAlg.h"
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"
#include "larcorealg/TestUtils/boost_unit_test_base.h"
#include "larcorealg/TestUtils/geometry_unit_test_base.h"

// art libraries
#include "fhiclcpp/ParameterSet.h"

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed; in the specific, the type of the
// channel mapping and a proper test name, used for output only;
// BasicGeometryEnvironmentConfiguration can read the configuration file name from command
// line, and BoostCommandLineConfiguration<> makes it initialize in time for Boost to
// catch it when instanciating the fixture.
struct StandardGeometryConfiguration
  : public testing::BoostCommandLineConfiguration<testing::BasicGeometryEnvironmentConfiguration> {
  StandardGeometryConfiguration() { SetApplicationName("GeometryStandardChannelMappingTest"); }
};

/**
 * Our fixture is based on GeometryTesterEnvironment, configured with the object above.
 * It provides to the testing environment:
 * - `Tester()`, returning a configured instance of the test algorithm;
 * - `GlobalTester()`, (static) returning a global configured instance of the test
 *   algorithm.
 *
 * The testing::TestSharedGlobalResource<> facility provides a singleton instance of the
 * tester algorithm, shared through the whole program and with this object too.
 *
 * This sharing allows the fixture to be used either as global or as per-suite.  In the
 * former case, the BOOST_AUTO_TEST_CASE's will access the global test algotithm instance
 * through the static call to `GeometryStandardChannelMappingTestFixture::GlobalTester()`;
 * in the latter case, it will access the local tester via the member function `Tester()`.
 * In this case, whether `GlobalTester()` and `Tester()` point to the same tester depends
 * on Boost unit test implementation.
 */
class GeometryStandardChannelMappingTestFixture
  : private testing::GeometryTesterEnvironment<StandardGeometryConfiguration,
                                               geo::GeoObjectSorterStandard> {
  using Tester_t = geo::ChannelMapStandardTestAlg;
  using TesterRegistry_t = testing::TestSharedGlobalResource<Tester_t>;

public:
  GeometryStandardChannelMappingTestFixture()
  {
    AcquireProvider(lar::standalone::SetupReadout({}, Geometry()));
    tester_ptr = std::make_shared<Tester_t>(Geometry(), Provider<geo::WireReadoutGeom>());
    TesterRegistry_t::ProvideDefaultSharedResource(tester_ptr);
  }

  Tester_t& Tester() { return *tester_ptr; }
  static Tester_t& GlobalTester() { return TesterRegistry_t::Resource(); }

private:
  std::shared_ptr<Tester_t> tester_ptr; ///< our tester (may be shared)
};

//------------------------------------------------------------------------------
//---  The tests
//---

BOOST_GLOBAL_FIXTURE(GeometryStandardChannelMappingTestFixture);

BOOST_AUTO_TEST_CASE(TPCsetMappingTestCase)
{
  GeometryStandardChannelMappingTestFixture::GlobalTester().TPCsetMappingTest();
}

BOOST_AUTO_TEST_CASE(ROPMappingTestCase)
{
  GeometryStandardChannelMappingTestFixture::GlobalTester().ROPMappingTest();
}

BOOST_AUTO_TEST_CASE(ChannelMappingTestCase)
{
  GeometryStandardChannelMappingTestFixture::GlobalTester().ChannelMappingTest();
}
