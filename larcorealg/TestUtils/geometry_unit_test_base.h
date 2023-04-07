/**
 * @file   geometry_unit_test_base.h
 * @brief  Base class for objects initializing a geometry
 * @date   May 7th, 2015
 * @author petrillo@fnal.gov
 * @see    unit_test_base.h
 *
 * Provides an environment for easy set up of a Geometry-aware test.  Keep in mind that,
 * as much as I could push on flexibility, the channel mapping algorithm must be
 * hard-coded and, if using Boost unit test, the configuration file location must be hard
 * coded too (or you can use the provided configuration).
 *
 * For an example of usage, see larcore/test/Geometry/geometry_iterator_test.cxx
 *
 * The standard TesterEnvironment<> class can't handle geo::GeometryCore.  The reason is
 * twofold: for once, WireReadout service is not factorized, so we need to choose
 * explicitly the WireReadoutGeom implementation (here this is obtained by a template
 * argument). Another is that GeometryCore both is required and requires WireReadoutGeom
 * to have a complete initialisation. There are ways to overcome the issue at the cost of
 * added complication.
 *
 * Currently provides:
 * - BasicGeometryEnvironmentConfiguration: a test environment configuration
 * - GeometryTesterEnvironment: a prepacked geometry-aware test environment
 *
 */

#ifndef TEST_GEOMETRY_UNIT_TEST_BASE_H
#define TEST_GEOMETRY_UNIT_TEST_BASE_H

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderStandard.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/TestUtils/unit_test_base.h"

// utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// CET libraries
#include "cetlib/search_path.h"

// C/C++ standard libraries
#include <iostream> // for output before message facility is set up
#include <map>
#include <memory> // std::unique_ptr<>
#include <string>
#include <type_traits>

namespace testing {

  /** **************************************************************************
   * @brief Class holding a configuration for a test environment
   * @tparam CHANNELMAP the class used for channel mapping
   * @see GeometryTesterEnvironment
   *
   * This class needs to be fully constructed by the default constructor in order to be
   * useful as Boost unit test fixture.  It is supposed to be passed as a template
   * parameter to another class that can store an instance of it and extract configuration
   * information from it.
   */
  struct BasicGeometryEnvironmentConfiguration : public BasicEnvironmentConfiguration {

    /// Default constructor; this is what is used in Boost unit test
    BasicGeometryEnvironmentConfiguration() : BasicEnvironmentConfiguration() { LocalInit(); }

    /// Constructor: acquires parameters from the command line
    BasicGeometryEnvironmentConfiguration(int argc, char** argv)
      : BasicEnvironmentConfiguration(argc, argv)
    {
      LocalInit();
    }

    /// Constructor; accepts the name as parameter
    BasicGeometryEnvironmentConfiguration(std::string name) : BasicEnvironmentConfiguration(name)
    {
      LocalInit();
    }

    BasicGeometryEnvironmentConfiguration(int argc, char** argv, std::string name)
      : BasicEnvironmentConfiguration(argc, argv, name)
    {
      LocalInit();
    }

    /// @{
    /// @name Access to configuration
    /// FHiCL path for the geometry configuration
    std::string GeometryParameterSetPath() const
    {
      return ServiceParameterSetPath(GeometryServiceName());
    }

    /// A string describing the default parameter set to configure geometry
    std::string DefaultGeometryConfiguration() const
    {
      return DefaultServiceConfiguration(GeometryServiceName());
    }

    ///@}

    /// @{
    /// @name Set configuration

    /// Sets the FHiCL path for the geometry configuration
    void SetGeometryParameterSetPath(std::string path)
    {
      SetServiceParameterSetPath(GeometryServiceName(), path);
    }

    /// Sets a string describing the default parameter set to configure geometry
    void SetDefaultGeometryConfiguration(std::string cfg)
    {
      AddDefaultServiceConfiguration(GeometryServiceName(), cfg);
    }

    ///@}

    /// Returns the name of the service
    static std::string GeometryServiceName() { return "Geometry"; }

  protected:
    /// Initialize with some default values
    void LocalInit()
    {
      SetDefaultGeometryConfiguration(R"(
          SurfaceY:        200.  # in cm, vertical distance to the surface
          Name:            "lartpcdetector"
          GDML:            "LArTPCdetector.gdml"
          ROOT:            "LArTPCdetector.gdml"
          SortingParameters: { tool_type: GeoObjectSorterStandard }
          )");
    } // LocalInit()

  }; // class BasicGeometryEnvironmentConfiguration<>

  /** **************************************************************************
   * @brief Environment for a geometry test
   * @tparam ConfigurationClass a class providing compile-time configuration
   *
   * The test environment is set up on construction.
   *
   * The environment provides:
   * - Geometry() method to access geometry (as a constant pointer)
   * - Parameters() method returning the complete FHiCL configuration
   * - TesterParameters() method returning the configuration for the test
   *
   * This class or a derived one can be used as global fixture for unit tests that require
   * the presence of geometry (in the form of geo::GeometryCore instance).
   *
   * Unfortunately Boost does not give any control on the initialization of the object, so
   * everything must be ready to go as hard coded.  The ConfigurationClass class tries to
   * alleviate that.  That is another, small static class that GeometryTesterEnvironment
   * uses to get its parameters.
   *
   * The requirements for the ConfigurationClass are:
   * - `ChannelMapClass`: concrete type of channel mapping algorithm class
   * - `std::string ApplicationName()`: the application name
   * - `std::string ConfigurationPath()`: path to the configuration file
   * - `std::string GeometryParameterSetPath()`: FHiCL path to the configuration of the
   *   geometry; in art is `"services.Geometry"`
   * - `std::string TesterParameterSetPath()`: FHiCL path to the configuration of the
   *   geometry
   * - `std::string DefaultGeometryConfiguration()` returning a FHiCL string to be parsed
   *   to extract the default geometry configuration
   * - `std::string DefaultTesterConfiguration()` returning a FHiCL string to be parsed to
   *   extract the default test configuration
   *
   * Whether the configuration comes from a file or from the two provided defaults, it is
   * always expected within the parameter set paths: the default configuration must also
   * contain that path.
   *
   * Note that there is no room for polymorphism here since the setup happens on
   * construction.  Some methods are declared virtual in order to allow to tweak some
   * steps of the set up, but it's not trivial to create a derived class that works
   * correctly: the derived class must declare a new default constructor, and that default
   * constructor must call the protected constructor
   * (GeometryTesterEnvironment<ConfigurationClass>(no_setup))
   */
  template <typename ConfigurationClass, typename ObjectSorter>
  class GeometryTesterEnvironment : public TesterEnvironment<ConfigurationClass> {
    using TesterEnvironment_t = TesterEnvironment<ConfigurationClass>;

  public:
    /**
     * @brief Constructor: sets everything up and declares the test started
     *
     * The configuration is from a default-constructed ConfigurationClass.
     * This is suitable for use as Boost unit test fixture.
     */
    GeometryTesterEnvironment(bool bSetup = true) : TesterEnvironment_t(false)
    {
      if (bSetup) Setup();
    }

    //@{
    /**
     * @brief Setup from a configuration
     * @param configurer an instance of ConfigurationClass
     *
     * The configuration is from the specified configurer class.
     *
     * This constructor allows to use a non-default-constructed configuration.  This can't
     * be used (at best of my knowledge) when using this class as Boost unit test fixture.
     *
     * In the r-value-reference constructor, the configurer is moved.
     */
    GeometryTesterEnvironment(ConfigurationClass const& cfg_obj, bool bSetup = true)
      : TesterEnvironment_t(cfg_obj, false)
    {
      if (bSetup) Setup();
    }
    GeometryTesterEnvironment(ConfigurationClass&& cfg_obj, bool bSetup = true)
      : TesterEnvironment_t(std::move(cfg_obj), false)
    {
      if (bSetup) Setup();
    }
    //@}

    geo::GeometryCore const* Geometry() const noexcept
    {
      return this->template Provider<geo::GeometryCore>();
    }

    /// Destructor: closing remarks
    virtual ~GeometryTesterEnvironment();

  protected:
    /// The complete initialization, ran at construction by default
    void Setup();

    /// Creates a new geometry
    virtual std::unique_ptr<geo::GeometryCore> CreateNewGeometry() const;

  private:
    ConfigurationClass config; ///< instance of the configurer

  }; // class GeometryTesterEnvironment<>

  template <typename ObjectSorter>
  auto createSorter(fhicl::ParameterSet const& pset [[maybe_unused]])
  {
    if constexpr (std::is_constructible_v<ObjectSorter>) {
      return std::make_unique<ObjectSorter>();
    }
    else if constexpr (std::is_constructible_v<ObjectSorter, fhicl::ParameterSet>) {
      return std::make_unique<ObjectSorter>(pset);
    }
  }

  //****************************************************************************
  template <typename ConfigurationClass, typename ObjectSorter>
  GeometryTesterEnvironment<ConfigurationClass, ObjectSorter>::~GeometryTesterEnvironment()
  {
    mf::LogInfo("Test") << config.ApplicationName() << " completed.";
  }

  /** **************************************************************************
   * @brief Sets the geometry of the standard detector up
   *
   * This function sets up the geometry according to the provided information:
   * - the configuration must contain enough information to locate the geometry
   *   description file
   * - we trust that that geometry works well with the ChannelMapClass specified in
   *   ConfigurationClass
   *
   */
  template <typename ConfigurationClass, typename ObjectSorter>
  std::unique_ptr<geo::GeometryCore>
  GeometryTesterEnvironment<ConfigurationClass, ObjectSorter>::CreateNewGeometry() const
  {
    std::string ProviderParameterSetPath = this->Config().GeometryParameterSetPath();

    //
    // create the new geometry service provider
    //
    fhicl::ParameterSet ProviderConfig =
      this->Parameters().template get<fhicl::ParameterSet>(ProviderParameterSetPath);

    auto new_geom = std::make_unique<geo::GeometryCore>(
      ProviderConfig,
      std::make_unique<geo::GeometryBuilderStandard>(
        ProviderConfig.get<fhicl::ParameterSet>("Builder", {})),
      createSorter<ObjectSorter>(ProviderConfig.get<fhicl::ParameterSet>("SortingParameters", {})));

    std::string RelativePath = ProviderConfig.get<std::string>("RelativePath", "");

    std::string GDMLFileName = RelativePath + ProviderConfig.get<std::string>("GDML"),
                ROOTFileName = RelativePath + ProviderConfig.get<std::string>("ROOT");

    // Search all reasonable locations for the geometry file; we see if by any chance
    // art's FW_SEARCH_PATH directory is set and try there; if not, we do expect the path
    // to be complete enough for ROOT to cope.
    cet::search_path sp("FW_SEARCH_PATH");

    std::string ROOTfile;
    if (!sp.find_file(ROOTFileName, ROOTfile)) ROOTfile = ROOTFileName;

    // we really don't care of GDML file, since we are not going to run Geant4
    std::string GDMLfile;
    if (!sp.find_file(GDMLFileName, GDMLfile)) {
      mf::LogWarning("CreateNewGeometry") << "GDML file '" << GDMLfile << "' not found.";
    }

    // initialize the geometry with the files we have found
    new_geom->LoadGeometryFile(GDMLfile, ROOTfile);

    return new_geom;
  }

  template <typename ConfigurationClass, typename ObjectSorter>
  void GeometryTesterEnvironment<ConfigurationClass, ObjectSorter>::Setup()
  {
    TesterEnvironment_t::Setup();
    this->AcquireProvider(CreateNewGeometry());
    mf::LogInfo("Test") << config.ApplicationName() << " Geometry setup complete.";
  }

} // namespace testing

#endif // TEST_GEOMETRY_UNIT_TEST_BASE_H
