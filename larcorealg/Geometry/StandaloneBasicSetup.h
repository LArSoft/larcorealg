/**
 * @file   larcorealg/Geometry/StandaloneBasicSetup.h
 * @brief  Collection of functions for quick setup of basic facilities.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 22, 2017
 *
 * Currently the following functionality is provided:
 * * configuration parsing (`ParseConfiguration()`)
 * * message facility service setup (`SetupMessageFacility()`)
 *
 * A complete basic setup for `my test` may look like:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * using namespace lar::standalone;
 *
 * // parse the configuration file taken from the first command line argument:
 * fhicl::ParameterSet pset = ParseConfiguration(argv[1]);
 *
 * // set up the message facility service
 * SetupMessageFacility(pset, "my test");
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Additional setup utilities may be available in `lar::standalone` namespace.
 * Also, some service providers can be set up using `testing`setupProvider()`
 * functions (the function itself is defined in `ProviderTestHelpers.h`).
 *
 * Currently this is a header-only library.
 *
 */

#ifndef LARCOREALG_GEOMETRY_STANDALONEBASICSETUP_H
#define LARCOREALG_GEOMETRY_STANDALONEBASICSETUP_H

// utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h" // cet::filepath_lookup_after1

// C/C++ standard libraries
#include <string>


namespace lar {

  /// Utilities for use in an environment without _art_.
  namespace standalone {

    //--------------------------------------------------------------------------
    /**
     * @brief Parses a FHiCL configuration file.
     * @param configPath the path to the configuration file to be parsed
     * @param lookupPolicy class performing the proper lookup for FHiCL files
     * @return a parameter set containing the whole configuration from the file
     *
     * Some lookup policies are provided in cetlib (see
     * `cetlib/filepath_maker.h` file).
     */
    fhicl::ParameterSet ParseConfiguration
      (std::string configPath, cet::filepath_maker& lookupPolicy);


    //--------------------------------------------------------------------------
    /**
     * @brief Parses a FHiCL configuration file.
     * @param configPath the path to the configuration file to be parsed
     * @return a parameter set containing the whole configuration from the file
     *
     * The lookup policy for finding the FHiCL files is such that:
     * * `configPath` must be readily available: no special lookup is performed
     * * FHiCL files included (directly or indirectly) by `configPath` are
     *   searched for in the path list specified in the environment variable
     *   `FHICL_FILE_PATH`
     *
     */
    fhicl::ParameterSet ParseConfiguration(std::string configPath);


    //--------------------------------------------------------------------------
    /**
     * @brief Sets up the message facility service.
     * @param pset global configuration parameter set
     * @param applName (default: `"standalone"`) name of running the application
     *
     * The configuration is read from the path `services.message` (as for the
     * standard _art_ behaviour). Any configuration working in _art_ is expected
     * to work here as well.
     *
     * Technical details:
     * * "context singlet" is set to `main` (`mf::SetContextSinglet()`)
     * * "context iteration" is set to empty (`mf::SetContextIteration()`)
     */
    void SetupMessageFacility
      (fhicl::ParameterSet const& pset, std::string applName = "standalone");


    //--------------------------------------------------------------------------

  } // namespace standalone
} // namespace lar


//------------------------------------------------------------------------------
//---  inline implementation
//---
//------------------------------------------------------------------------------
inline fhicl::ParameterSet lar::standalone::ParseConfiguration
  (std::string configPath, cet::filepath_maker& lookupPolicy)
{
  fhicl::ParameterSet pset;
  pset = fhicl::ParameterSet::make(configPath, lookupPolicy);
  return pset;
} // ParseConfiguration(string, filepath_maker)


//------------------------------------------------------------------------------
inline fhicl::ParameterSet lar::standalone::ParseConfiguration
  (std::string configPath)
{
  cet::filepath_lookup_after1 policy("FHICL_FILE_PATH");
  return ParseConfiguration(configPath, policy);
} // ParseConfiguration(string)


//------------------------------------------------------------------------------
inline void lar::standalone::SetupMessageFacility
  (fhicl::ParameterSet const& pset, std::string applName /* = "standalone" */)
{
  mf::StartMessageFacility(pset.get<fhicl::ParameterSet>("services.message"));
  mf::SetApplicationName(applName);
  mf::SetContextSinglet("main");
  mf::SetContextIteration("");
} // lar::standalone::SetupMessageFacility()


//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_STANDALONEBASICSETUP_H
