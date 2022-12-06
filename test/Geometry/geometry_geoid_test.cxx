/**
 * @file   geometry_geoid_test.cxx
 * @brief  Unit test for geometry ID consistency on a standard detector
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   October 31, 2016
 *
 * Usage: just run the executable.
 * Or plug a FHiCL file in the command line.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate
// a main() function; Boost is pulled in by boost_unit_test_base.h
#define BOOST_TEST_MODULE GeometryGeoIDTest

// LArSoft libraries
#include "GeometryGeoIDTestAlg.h"
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/TestUtils/boost_unit_test_base.h"
#include "larcorealg/TestUtils/geometry_unit_test_base.h"

// utility libraries

// C/C++ standard libraries

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// in the specific, the type of the channel mapping and a proper test name,
// used for output only; BasicGeometryEnvironmentConfiguration can read the
// configuration file name from command line, and
// BoostCommandLineConfiguration<> makes it initialize in time for Boost
// to catch it when instanciating the fixture.
struct StandardGeometryConfiguration
  : public testing::BoostCommandLineConfiguration<testing::BasicGeometryEnvironmentConfiguration> {
  /// Constructor: overrides the application name
  StandardGeometryConfiguration() { SetApplicationName("GeometryGeoIDUnitTest"); }
}; // class StandardGeometryConfiguration

/*
 * Our fixture is based on GeometryTesterEnvironment, configured with the object
 * above.
 * It provides to the testing environment:
 * - `Tester()`, returning a configured instance of the test algorithm;
 * - `GlobalTester()`, (static) returning a global configured instance of the
 *   test algorithm.
 *
 * The testing::TestSharedGlobalResource<> facility provides a singleton
 * instance of the tester algorithm, shared through the whole program and with
 * this object too.
 *
 * This sharing allows the fixture to be used either as global or as per-suite.
 * In the former case, the BOOST_AUTO_TEST_CASE's will access the global test
 * algorithm instance through the static call to
 * `GeometryGeoIDTestFixture::GlobalTester()`; in the latter case, it will
 * access the local tester via the member function `Tester()`.
 * In this case, whether `GlobalTester()` and `Tester()` point to the same
 * tester depends on Boost unit test implementation.
 */
class GeometryGeoIDTestFixture
  : private testing::GeometryTesterEnvironment<StandardGeometryConfiguration,
                                               geo::GeoObjectSorterStandard> {
  using Tester_t = geo::GeometryGeoIDTestAlg;

  using TesterRegistry_t = testing::TestSharedGlobalResource<Tester_t>;

public:
  /// Constructor: initialize the tester with the Geometry from base class
  GeometryGeoIDTestFixture()
  {
    // create a new tester
    tester_ptr = std::make_shared<Tester_t>(Geometry());
    // if no tester is default yet, share ours:
    TesterRegistry_t::ProvideDefaultSharedResource(tester_ptr);
  }

  /// Retrieves the local tester
  Tester_t& Tester() { return *(tester_ptr.get()); }

  /// Retrieves the global tester
  static Tester_t& GlobalTester() { return TesterRegistry_t::Resource(); }

private:
  std::shared_ptr<Tester_t> tester_ptr; ///< our tester (may be shared)
};                                      // class GeometryGeoIDTestFixture

//------------------------------------------------------------------------------
//---  The tests
//---
//
// BOOST_GLOBAL_FIXTURE does not provide tester access, so one has to get it
// as GeometryGeoIDTestFixture::GlobalTester(); on the other hand, the
// fixture is initialised only when a new global one is explicitly created.
//

BOOST_GLOBAL_FIXTURE(GeometryGeoIDTestFixture);

BOOST_AUTO_TEST_CASE(CryostatGeoIDTest)
{
  GeometryGeoIDTestFixture::GlobalTester().CryostatGeoIDTest();
} // BOOST_AUTO_TEST_CASE( CryostatGeoIDTest )

BOOST_AUTO_TEST_CASE(TPCGeoIDTest)
{
  GeometryGeoIDTestFixture::GlobalTester().TPCGeoIDTest();
} // BOOST_AUTO_TEST_CASE( TPCGeoIDTest )

BOOST_AUTO_TEST_CASE(PlaneGeoIDTest)
{
  GeometryGeoIDTestFixture::GlobalTester().PlaneGeoIDTest();
} // BOOST_AUTO_TEST_CASE( PlaneGeoIDTest )

BOOST_AUTO_TEST_CASE(WireGeoIDTest)
{
  GeometryGeoIDTestFixture::GlobalTester().WireGeoIDTest();
} // BOOST_AUTO_TEST_CASE( WireGeoIDTest )

// BOOST_AUTO_TEST_SUITE_END()
