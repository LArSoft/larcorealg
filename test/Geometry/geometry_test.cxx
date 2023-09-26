/**
 * @file   geometry_test.cxx
 * @brief  Unit test for geometry functionalities on a standard detector
 * @date   May 5th, 2015
 * @author petrillo@fnal.gov
 *
 * Usage:
 *   geometry_test  [ConfigurationFile [GeometryTestParameterSet]]
 *
 * By default, GeometryTestParameterSet is set to "physics.analysers.geotest".
 *
 * This unit test uses geometry_unit_test_base.h to build an environment with a geometry
 * set up.  For an example of use with Boost unit test module, see
 * geometry_iterator_test.cxx .
 */

// LArSoft libraries
#include "GeometryTestAlg.h"
#include "larcorealg/Geometry/AuxDetGeoObjectSorterStandard.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"
#include "larcorealg/TestUtils/geometry_unit_test_base.h"

// art libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed; we use an existing class provided
// for this purpose, since our test environment allows us to tailor it at run time.
using StandardGeometryConfiguration = testing::BasicGeometryEnvironmentConfiguration;

/*
 * GeometryTesterFixture, configured with the object above, is used in a
 * non-Boost-unit-test context.
 * It provides:
 * - `geo::GeometryCore const* Geometry()`
 * - `geo::GeometryCore const* GlobalGeometry()` (static member)
 */
using StandardGeometryTestEnvironment =
  testing::GeometryTesterEnvironment<StandardGeometryConfiguration, geo::GeoObjectSorterStandard>;

//------------------------------------------------------------------------------
//---  The tests
//---

/** ****************************************************************************
 * @brief Runs the test
 * @param argc number of arguments in argv
 * @param argv arguments to the function
 * @return number of detected errors (0 on success)
 * @throw cet::exception most of error situations throw
 *
 * The arguments in argv are:
 * 0. name of the executable ("Geometry_test")
 * 1. path to the FHiCL configuration file
 * 2. FHiCL path to the configuration of the geometry test
 *    (default: physics.analyzers.geotest)
 * 3. FHiCL path to the configuration of the geometry
 *    (default: services.Geometry)
 */
//------------------------------------------------------------------------------
int main(int argc, char const** argv)
{
  StandardGeometryConfiguration config("geometry_test");
  config.SetMainTesterParameterSetName("geotest");

  // parameter parsing
  int iParam = 0;

  // first argument: configuration file (mandatory)
  if (++iParam < argc) config.SetConfigurationPath(argv[iParam]);

  // second argument: path of the parameter set for geometry test configuration
  // (optional; default: "physics.analysers.geotest")
  if (++iParam < argc) config.SetMainTesterParameterSetPath(argv[iParam]);

  // third argument: path of the parameter set for geometry configuration
  // (optional; default: "services.Geometry" from the inherited object)
  if (++iParam < argc) config.SetGeometryParameterSetPath(argv[iParam]);

  // testing environment setup
  StandardGeometryTestEnvironment TestEnvironment(config);
  auto const wireReadoutGeom = lar::standalone::SetupReadout(
    TestEnvironment.ServiceParameters("WireReadout"), TestEnvironment.Geometry());
  auto const auxDetGeom =
    lar::standalone::SetupAuxDetGeometry(TestEnvironment.ServiceParameters("AuxDetGeometry"));

  // run the test algorithm

  // 1. we initialize it from the configuration in the environment,
  geo::GeometryTestAlg Tester{TestEnvironment.Geometry(),
                              wireReadoutGeom.get(),
                              auxDetGeom.get(),
                              TestEnvironment.TesterParameters()};

  // 2. then we run it!
  unsigned int nErrors = Tester.Run();

  // 3. And finally we cross fingers.
  if (nErrors > 0) { mf::LogError("geometry_test") << nErrors << " errors detected!"; }

  return nErrors;
} // main()
