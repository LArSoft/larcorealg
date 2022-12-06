/**
 * @file   driftvolumes_test.cxx
 * @brief  Test for neighbourhood discovery in simple geometries.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 13, 2017
 *
 * Usage:
 *
 *     driftvolumes_test configuration.fcl
 *
 * The configuration file must contain a geometry service provider configuration
 * under `services.Geometry`, That geometry configuration must use the standard
 * channel mapping.
 */

// LArSoft libraries
#include "larcorealg/CoreUtils/DumpUtils.h"      // lar::dump::vector3D()
#include "larcorealg/Geometry/DriftPartitions.h" // BuildDriftVolumes()
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/TestUtils/geometry_unit_test_base.h"

// utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// we use an existing class provided for this purpose, since our test
// environment allows us to tailor it at run time.
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
 * 2. FHiCL path to the configuration of the geometry
 *    (default: services.Geometry)
 *
 */
//------------------------------------------------------------------------------
int main(int argc, char const** argv)
{

  StandardGeometryConfiguration config("driftvolumes_test");

  //
  // parameter parsing
  //
  int iParam = 0;

  // first argument: configuration file (mandatory)
  if (++iParam < argc) config.SetConfigurationPath(argv[iParam]);

  // second argument: path of the parameter set for geometry configuration
  // (optional; default: "services.Geometry" from the inherited object)
  if (++iParam < argc) config.SetGeometryParameterSetPath(argv[iParam]);

  //
  // testing environment setup
  //
  StandardGeometryTestEnvironment TestEnvironment(config);
  auto const& geom = *TestEnvironment.Geometry();

  //
  // run the test algorithm
  //

  unsigned int nErrors = 0;
  for (auto const& cryo : geom.Iterate<geo::CryostatGeo>()) {
    auto partition = geo::buildDriftVolumes(cryo);
    mf::LogVerbatim("driftvolumes_test") << "Partition for cryostat " << cryo.ID() << ":";
    partition.print(mf::LogVerbatim("driftvolumes_test"));

    //
    // test that the partition topology is correct
    //
    for (geo::TPCGeo const& TPC : geom.Iterate<geo::TPCGeo>(cryo.ID())) {

      auto const& center = TPC.GetCenter();

      auto where = partition.TPCat(center);
      if (!where) {
        mf::LogProblem("driftvolumes_test")
          << "Center of TPC " << TPC.ID() << " " << lar::dump::vector3D(center)
          << " not assigned to any TPC!";
        ++nErrors;
      }
      else if (where->ID() != TPC.ID()) {
        mf::LogProblem log("driftvolumes_test");
        log << "Center of TPC " << TPC.ID() << " " << lar::dump::vector3D(center)
            << " assigned to TPC " << where->ID() << ":\n";
        where->PrintTPCInfo(log, "  ", /* verbosity */ 5);
        ++nErrors;
      }

      //
      // test ownership of a uniform distribution of points inside the TPC
      //
      constexpr int nXsteps = 5, nYsteps = 5, nZsteps = 5;
      const double xstep = TPC.HalfSizeX() / (nXsteps + 1);
      const double ystep = TPC.HalfSizeY() / (nYsteps + 1);
      const double zstep = TPC.HalfSizeZ() / (nZsteps + 1);
      for (int xs = -nXsteps; xs <= nXsteps; ++xs) {
        double const x = center.X() + xs * xstep;
        for (int ys = -nYsteps; ys <= nYsteps; ++ys) {
          double const y = center.Y() + ys * ystep;
          for (int zs = -nZsteps; zs <= nZsteps; ++zs) {
            double const z = center.Z() + zs * zstep;

            auto where = partition.TPCat({x, y, z});
            if (!where) {
              mf::LogProblem("driftvolumes_test")
                << "Point " << lar::dump::vector3D(center) << " within TPC " << TPC.ID() << " ("
                << xs << "," << ys << "," << zs << ") not assigned to any TPC!";
              ++nErrors;
            }
            else if (where->ID() != TPC.ID()) {
              mf::LogProblem log("driftvolumes_test");
              log << "Point " << lar::dump::vector3D(center) << " within TPC " << TPC.ID() << " ("
                  << xs << "," << ys << "," << zs << ") assigned to TPC " << where->ID() << ":\n";
              where->PrintTPCInfo(log, "  ", /* verbosity */ 5);
              ++nErrors;
            }

          } // for zs
        }   // for ys
      }     // for xs

    } // for TPCs

  } // for

  // 4. And finally we cross fingers.
  if (nErrors > 0) { mf::LogError("geometry_test") << nErrors << " errors detected!"; }

  return nErrors;
} // main()
