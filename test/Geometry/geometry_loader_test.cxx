/**
 * @file   geometry_loader_test.cxx
 * @brief  Unit test for geometry functionalities on a standard detector.
 * @date   June 22, 2017
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 * Usage:
 *
 *     geometry_loader_test  ConfigurationFile [test configuration FHiCL path]
 *
 * The configuration file path must be complete, i.e. it must point directly to the
 * configuration file.  By default, the test configuration is expected in
 * `"physics.analysis.geotest"`.
 *
 * This unit test uses `StaticLoadGeometry()` to set up the geometry.
 */

// LArSoft libraries
#include "GeometryTestAlg.h"
#include "larcorealg/Geometry/StandaloneBasicSetup.h"    // SetupMessageFacility()...
#include "larcorealg/Geometry/StandaloneGeometrySetup.h" // SetupGeometry()

// utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <stdexcept> // std::runtime_error
#include <string>

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
 *    (default: physics.analysers.geotest)
 * 3. FHiCL path to the configuration of the geometry
 *    (default: services.Geometry)
 */

//------------------------------------------------------------------------------
int main(int argc, char const** argv)
{
  // parameter parsing
  std::string configPath;
  std::string geoTestConfigPath = "physics.analyzers.geotest";

  int iParam = 0;

  // first argument: configuration file (mandatory)
  if (++iParam < argc)
    configPath = argv[iParam];
  else
    throw std::runtime_error("No configuration file specified.");

  // second argument: test configuration path
  if (++iParam < argc) geoTestConfigPath = argv[iParam];

  // 1. testing environment setup
  using namespace lar::standalone;

  // parse a configuration file
  fhicl::ParameterSet pset = ParseConfiguration(configPath);

  // set up message facility
  SetupMessageFacility(pset, "geometry_loader_test");
  mf::SetContextIteration("setup");

  // set up geometry
  auto geom = SetupGeometry(pset.get<fhicl::ParameterSet>("services.Geometry"));

  auto wireReadoutGeom =
    SetupReadout(pset.get<fhicl::ParameterSet>("services.WireGeom", {}), geom.get());

  // update the context string for the messages
  mf::SetContextIteration("run");

  // 2. prepare the test algorithm
  geo::GeometryTestAlg Tester(
    geom.get(), wireReadoutGeom.get(), pset.get<fhicl::ParameterSet>(geoTestConfigPath));

  // 3. then we run it!
  unsigned int nErrors = Tester.Run();

  // 4. And finally we cross fingers.
  mf::SetContextIteration("end");
  if (nErrors > 0) { mf::LogError("geometry_test") << nErrors << " errors detected!"; }

  return nErrors;
} // main()
