/**
 * @file   geometry_iterator_test.cxx
 * @brief  Unit test for geometry iterators on a standard detector
 * @date   May 7th, 2015
 * @author petrillo@fnal.gov
 *
 * Usage: just run the executable.
 * Or plug a FHiCL file in the command line.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate a main()
// function; Boost is pulled in by boost_unit_test_base.h
#define BOOST_TEST_MODULE GeometryIteratorTest

// LArSoft libraries
#include "GeometryIteratorTestAlg.h"
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"
#include "larcorealg/TestUtils/boost_unit_test_base.h"
#include "larcorealg/TestUtils/geometry_unit_test_base.h"

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
  StandardGeometryConfiguration() { SetApplicationName("GeometryIteratorUnitTest"); }
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
 * through the static call to `GeometryIteratorTestFixture::GlobalTester()`; in the latter
 * case, it will access the local tester via the member function `Tester()`.  In this
 * case, whether `GlobalTester()` and `Tester()` point to the same tester depends on Boost
 * unit test implementation.
 */

class GeometryIteratorTestFixture
  : private testing::GeometryTesterEnvironment<StandardGeometryConfiguration,
                                               geo::GeoObjectSorterStandard> {
  using Tester_t = geo::GeometryIteratorTestAlg;
  using TesterRegistry_t = testing::TestSharedGlobalResource<Tester_t>;

public:
  /// Constructor: initialize the tester with the Geometry from base class
  GeometryIteratorTestFixture()
  {
    AcquireProvider(lar::standalone::SetupReadout({}, Geometry()));
    tester_ptr = std::make_shared<Tester_t>(Geometry(), Provider<geo::WireReadoutGeom>());
    TesterRegistry_t::ProvideDefaultSharedResource(tester_ptr);
  }

  /// Retrieves the local tester
  Tester_t& Tester() { return *tester_ptr; }

  /// Retrieves the global tester
  static Tester_t& GlobalTester() { return TesterRegistry_t::Resource(); }

private:
  std::shared_ptr<Tester_t> tester_ptr;
}; // class GeometryIteratorTestFixture

//------------------------------------------------------------------------------
//---  The tests
//---

BOOST_GLOBAL_FIXTURE(GeometryIteratorTestFixture);

BOOST_AUTO_TEST_CASE(CryostatIDIteratorsTest)
{
  GeometryIteratorTestFixture::GlobalTester().CryostatIDIteratorsTest();
}

BOOST_AUTO_TEST_CASE(CryostatIteratorsTest)
{
  GeometryIteratorTestFixture::GlobalTester().CryostatIteratorsTest();
}

BOOST_AUTO_TEST_CASE(TPCIDIteratorsTest)
{
  GeometryIteratorTestFixture::GlobalTester().TPCIDIteratorsTest();
}

BOOST_AUTO_TEST_CASE(TPCIteratorsTest)
{
  GeometryIteratorTestFixture::GlobalTester().TPCIteratorsTest();
}

BOOST_AUTO_TEST_CASE(PlaneIDIteratorsTest)
{
  GeometryIteratorTestFixture::GlobalTester().PlaneIDIteratorsTest();
}

BOOST_AUTO_TEST_CASE(PlaneIteratorsTest)
{
  GeometryIteratorTestFixture::GlobalTester().PlaneIteratorsTest();
}

BOOST_AUTO_TEST_CASE(WireIDIteratorsTest)
{
  GeometryIteratorTestFixture::GlobalTester().WireIDIteratorsTest();
}

BOOST_AUTO_TEST_CASE(WireIteratorsTest)
{
  GeometryIteratorTestFixture::GlobalTester().WireIteratorsTest();
}

// BOOST_AUTO_TEST_SUITE_END()
