/**
 * @file   unit_test_base.h
 * @brief  Base class for unit tests using FHiCL configuration
 * @date   December 1st, 2015
 * @author petrillo@fnal.gov
 *
 * Provides an environment for easy set up of a message-facility-aware test.
 *
 * For an example of how to expand it to host services,
 * see geometry_unit_test_base.h
 *
 * Currently provides:
 * - BasicEnvironmentConfiguration: a test environment configuration
 * - TestSharedGlobalResource: mostly internal use
 * - TesterEnvironment: a prepacked test environment with some provider support
 *
 * This is a pure template header. It will require the following libraries:
 *
 * * `MF_MessageLogger`
 * * `fhiclcpp`
 * * `cetlib`
 *
 */

#ifndef TEST_UNIT_TEST_BASE_H
#define TEST_UNIT_TEST_BASE_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ProviderPack.h"
#include "larcorealg/TestUtils/ProviderList.h"
#include "larcorealg/TestUtils/ProviderTestHelpers.h"

// utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/parse.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// CET libraries
#include "cetlib/filepath_maker.h"
#include "cetlib/search_path.h"

// C/C++ standard libraries
#include <iostream> // for output before message facility is set up
#include <map>
#include <memory>    // std::unique_ptr<>
#include <stdexcept> // std::logic_error
#include <string>
#include <type_traits> // std::add_rvalue_reference()
#include <utility>     // std::move(), std::forward()

namespace testing {

  namespace details {

    /// Reads and makes available the command line parameters
    class CommandLineArguments {
    public:
      /// Constructor: automatically parses from Boost arguments
      CommandLineArguments() { Clear(); }

      /// Constructor: parses from specified arguments
      CommandLineArguments(int argc, char** argv) { ParseArguments(argc, argv); }

      /// Parses arguments
      void ParseArguments(int argc, char** argv);

      /// Returns the name of the executable as started
      std::string Executable() const { return exec_name; }

      /// Returns the list of non-Boost-test arguments on the command line
      std::vector<std::string> const& Arguments() const { return args; }

      /// Returns whether we have arguments up to the iArg-th (0-based)
      bool hasArgument(size_t iArg) const { return iArg < args.size(); }

      /// Returns the value of the iArg-th (0-based; no range check!)
      std::string const& Argument(size_t iArg) const { return args[iArg]; }

    private:
      std::string exec_name;         ///< name of the test executable (from argv[0])
      std::vector<std::string> args; ///< command line arguments (from argv[0])

      /// Erases the stored arguments
      void Clear()
      {
        exec_name.clear();
        args.clear();
      }

    }; // class CommandLineArguments

    inline void CommandLineArguments::ParseArguments(int argc, char** argv)
    {
      Clear();
      if (argc == 0) return;

      exec_name = argv[0];

      args.resize(argc - 1);
      std::copy(argv + 1, argv + argc, args.begin());

    } // CommandLineArguments:ParseArguments()

    // forward declaration
    template <typename TestEnv, typename Pack, typename... Provs>
    struct ProviderPackFiller;

  } // namespace details

  /** **************************************************************************
   * @brief Class holding a configuration for a test environment
   *
   * This class needs to be fully constructed by the default constructor
   * in order to be useful as Boost unit test fixture.
   * It is supposed to be passed as a template parameter to another class
   * that can store an instance of it and extract configuration information
   * from it.
   */
  struct BasicEnvironmentConfiguration {

    /// Default constructor; this is what is used in Boost unit test
    BasicEnvironmentConfiguration() { DefaultInit(); }

    /// Constructor: acquires parameters from the command line
    BasicEnvironmentConfiguration(int argc, char** argv) : BasicEnvironmentConfiguration()
    {
      ParseCommandLine(argc, argv);
    }

    /// Constructor; accepts the name as parameter
    BasicEnvironmentConfiguration(std::string name) : BasicEnvironmentConfiguration()
    {
      SetApplicationName(name);
    }

    BasicEnvironmentConfiguration(int argc, char** argv, std::string name)
      : BasicEnvironmentConfiguration(argc, argv)
    {
      SetApplicationName(name);
    }

    /// @{
    /// @name Access to configuration
    /// Name of the application
    std::string ApplicationName() const { return appl_name; }

    /// Path to the configuration file
    std::string ConfigurationPath() const { return config_path; }

    /// FHiCL path for the configuration of the test algorithm
    std::string TesterParameterSetPath(std::string name) const
    {
      auto iPath = test_paths.find(name);
      return (iPath == test_paths.end()) ? ("physics.analyzers." + name) : iPath->second;
    }

    /// Name of the test algorithm instance
    std::string MainTesterParameterSetName() const { return main_test_name; }

    /// FHiCL path for the configuration of the test algorithm
    std::string MainTesterParameterSetPath() const
    {
      return MainTesterParameterSetName().empty() ?
               "" :
               TesterParameterSetPath(MainTesterParameterSetName());
    }

    /// FHiCL path for the configuration of the service
    std::string ServiceParameterSetPath(std::string name) const
    {
      auto iPath = service_paths.find(name);
      return (iPath == service_paths.end()) ? ("services." + name) : iPath->second;
    } // ServiceParameterSetPath()

    /// A string describing default parameter set to configure specified test
    std::string DefaultTesterConfiguration(std::string tester_name) const
    {
      return analyzers_default_cfg.at(tester_name);
    }

    /// A string describing the default parameter set to configure the test
    std::string DefaultServiceConfiguration(std::string service_name) const
    {
      return services_default_cfg.at(service_name);
    }

    /// A string describing the full default parameter set
    std::string DefaultConfiguration() const { return BuildDefaultConfiguration(); }

    /// Returns the name of the executable as started
    std::string ExecutablePath() const { return arguments.Executable(); }

    /// Returns the list of non-Boost-test arguments on the command line
    std::vector<std::string> const& EexcutableArguments() const { return arguments.Arguments(); }

    ///@}

    /// @{
    /// @name Set configuration

    /// Sets the name of the application
    void SetApplicationName(std::string name) { appl_name = name; }

    /// Sets the path to the configuration file
    void SetConfigurationPath(std::string path) { config_path = path; }

    /// Sets the FHiCL name for the configuration of the test algorithm
    void SetMainTesterParameterSetName(std::string name) { main_test_name = name; }

    /// Sets the FHiCL path for the configuration of a test algorithm
    void SetTesterParameterSetPath(std::string test_name, std::string path)
    {
      test_paths[test_name] = path;
    }

    /// Sets the FHiCL path for the configuration of the main test algorithm
    void SetMainTesterParameterSetPath(std::string path)
    {
      if (MainTesterParameterSetName().empty()) {
        throw std::logic_error("Request setting configuration of non-existent main tester");
      }
      SetTesterParameterSetPath(MainTesterParameterSetName(), path);
    }

    /// Sets the FHiCL path for the configuration of a test algorithm
    void SetServiceParameterSetPath(std::string service_name, std::string path)
    {
      service_paths[service_name] = path;
    }

    /// Adds a default configuration for the specified service
    void AddDefaultServiceConfiguration(std::string service_name, std::string service_cfg)
    {
      services_default_cfg[service_name] = service_cfg;
    }

    /// Adds a default configuration for the specified tester
    void AddDefaultTesterConfiguration(std::string tester_name, std::string tester_cfg)
    {
      analyzers_default_cfg[tester_name] = tester_cfg;
    }

    /// Adds a default configuration for the main tester
    void AddDefaultTesterConfiguration(std::string tester_cfg)
    {
      if (MainTesterParameterSetName().empty()) {
        throw std::logic_error("Request adding configuration of non-existent main tester");
      }
      AddDefaultTesterConfiguration(MainTesterParameterSetName(), tester_cfg);
    }

    ///@}

  protected:
    using ConfigurationMap_t = std::map<std::string, std::string>;
    using PathMap_t = std::map<std::string, std::string>;

    std::string appl_name;      ///< name of the application
    std::string config_path;    ///< configuration file path
    std::string main_test_name; ///< name of main test algorithm
    std::string main_test_path; ///< path of main test algorithm configuration

    /// Returns the default test name
    static std::string DefaultApplicationName() { return "Test"; }

    /// Configuration of all the services
    ConfigurationMap_t services_default_cfg;
    /// Configuration of all the analyzer modules
    ConfigurationMap_t analyzers_default_cfg;

    /// Set of paths for tester configuration
    PathMap_t test_paths;
    /// Set of paths for service configuration
    PathMap_t service_paths;

    /// Extracts arguments from the command line, uses first one as config path
    void ParseCommandLine(int argc, char** argv)
    {
      arguments.ParseArguments(argc, argv);
      if (arguments.hasArgument(0)) SetConfigurationPath(arguments.Argument(0)); // first argument
    }

    /// Initialize with some default values
    void DefaultInit()
    {
      SetApplicationName(DefaultApplicationName());
      SetMainTesterParameterSetName("");
      // a destination which will react to all messages from DEBUG up:
      AddDefaultServiceConfiguration("message",
                                     R"(
           debugModules: [ '*' ]
           destinations : {
             stdout: {
               type:      cout
               threshold: DEBUG
               categories: {
                 default: {
                   limit: -1
                 }
               } // categories
             } // stdout
             statistics: {type:cout}
           } // destinations
         )");
    } // DefaultInit()

    /// A string describing the full default parameter set
    std::string BuildDefaultServiceConfiguration() const
    {
      return BuildServiceConfiguration(services_default_cfg);
    }

    /// A string describing the full default parameter set
    std::string BuildDefaultTestConfiguration() const
    {
      return BuildTestConfiguration(analyzers_default_cfg);
    }

    /// A string describing the full default parameter set
    std::string BuildDefaultConfiguration() const
    {
      return BuildConfiguration(services_default_cfg, analyzers_default_cfg);
    }

    /// A string with the service section from service parameter sets
    static std::string BuildServiceConfiguration(ConfigurationMap_t const& services)
    {
      std::string cfg;
      cfg += "\nservices: {";
      for (auto const& service_info : services) {
        cfg += "\n  " + service_info.first + ": {";
        cfg += "\n" + service_info.second;
        cfg += "\n  } # " + service_info.first;
      } // for services
      cfg += "\n} # services"
             "\n";
      return cfg;
    } // BuildServiceConfiguration()

    /// A string with the physics section from analyzer parameter sets
    static std::string BuildTestConfiguration(ConfigurationMap_t const& analyzers)
    {
      std::string cfg;
      cfg += "\nphysics: {"
             "\n  analyzers: {";
      for (auto const& module_info : analyzers) {
        cfg += "\n  " + module_info.first + ": {";
        cfg += "\n" + module_info.second;
        cfg += "\n  } # " + module_info.first;
      } // for analyzers
      cfg += "\n  } # analyzers"
             "\n} # physics";
      return cfg;
    } // BuildServiceConfiguration()

    /// A string describing the full default parameter set
    static std::string BuildConfiguration(ConfigurationMap_t const& services,
                                          ConfigurationMap_t const& modules)
    {
      std::string cfg;
      cfg += BuildServiceConfiguration(services);
      cfg += BuildTestConfiguration(modules);
      return cfg;
    } // BuildConfiguration()

  private:
    details::CommandLineArguments arguments; ///< command line arguments

  }; // class BasicEnvironmentConfiguration<>

  /** **************************************************************************
   * @brief Utility class providing singleton objects to the derived classes
   * @tparam RES the type of object (include constantness if needed)
   *
   * The object is expected to be shared.
   */
  template <typename RES>
  class TestSharedGlobalResource {
    using Resource_t = RES;

  public:
    using ResourcePtr_t = std::shared_ptr<Resource_t>;

    /// @name Add and share resources
    /// @{

    /// Adds a shared resource to the resource registry
    static void AddSharedResource(std::string res_name, ResourcePtr_t res_ptr)
    {
      Resources[res_name] = res_ptr;
    }

    /// Adds a shared resource to the resource registry (empty name)
    static void AddDefaultSharedResource(ResourcePtr_t res_ptr)
    {
      AddSharedResource(std::string(), res_ptr);
    }

    /// Registers a shared resource only if none exists yet
    template <typename... Args>
    static ResourcePtr_t ProvideSharedResource(std::string res_name, ResourcePtr_t res_ptr)
    {
      if (hasResource(res_name)) return ResourcePtr_t();
      AddSharedResource(res_name, res_ptr);
      return res_ptr;
    }

    /// Creates a shared resource as default only if none exists yet
    template <typename... Args>
    static ResourcePtr_t ProvideDefaultSharedResource(ResourcePtr_t res_ptr)
    {
      return ProvideSharedResource(std::string(), res_ptr);
    }

    //@{
    /// Adds a shared resource only if it is old_res_ptr
    static bool ReplaceSharedResource(std::string res_name,
                                      Resource_t const* old_res_ptr,
                                      ResourcePtr_t res_ptr)
    {
      ResourcePtr_t current_res_ptr = ShareResource();
      if (current_res_ptr.get() != old_res_ptr) return false;
      AddSharedResource(res_name, res_ptr);
      return true;
    }
    static bool ReplaceSharedResource(std::string res_name,
                                      ResourcePtr_t old_res_ptr,
                                      ResourcePtr_t res_ptr)
    {
      return ReplaceSharedResource(res_name, old_res_ptr.get(), res_ptr);
    }
    //@}

    //@{
    /// Adds a shared resource as default resource only if it is old_res_ptr
    static bool ReplaceDefaultSharedResource(Resource_t const* old_res_ptr, ResourcePtr_t res_ptr)
    {
      return ReplaceSharedResource(std::string(), old_res_ptr, res_ptr);
    }
    static bool ReplaceDefaultSharedResource(ResourcePtr_t old_res_ptr, ResourcePtr_t res_ptr)
    {
      return ReplaceSharedResource(std::string(), old_res_ptr, res_ptr);
    }
    //@}

    /// Constructs and registers a new resource with a specified name
    template <typename... Args>
    static ResourcePtr_t CreateResource(std::string res_name, Args&&... args)
    {
      ResourcePtr_t res_ptr(new Resource_t(std::forward<Args>(args)...));
      AddSharedResource(res_name, res_ptr);
      return res_ptr;
    }

    /// Constructs and registers a new resource with no name
    template <typename... Args>
    static void CreateDefaultResource(Args&&... args)
    {
      CreateResource(std::string(), std::forward<Args>(args)...);
    }

    /// Creates a shared resource only if none exists yet
    template <typename... Args>
    static ResourcePtr_t ProposeSharedResource(std::string res_name, Args&&... args)
    {
      return hasResource(res_name) ? ResourcePtr_t() :
                                     CreateResource(res_name, std::forward<Args>(args)...);
    }

    /// Creates a shared resource as default only if none exists yet
    template <typename... Args>
    static ResourcePtr_t ProposeDefaultSharedResource(Args&&... args)
    {
      return ProposeSharedResource(std::string(), std::forward<Args>(args)...);
    }

    /// @}

    /// @name Resource access
    /// @{

    /// Returns whether a resource exists
    /// @throws std::out_of_range if not available
    static bool hasResource(std::string name = "")
    {
      auto iRes = Resources.find(name);
      return (iRes != Resources.end()) && bool(iRes->second);
    }

    /// Retrieves the specified resource for sharing (nullptr if none)
    static ResourcePtr_t ShareResource(std::string name = "")
    {
      auto iRes = Resources.find(name);
      return (iRes == Resources.end()) ? ResourcePtr_t() : iRes->second;
    }

    /// Retrieves the specified resource, or throws if not available
    static Resource_t& Resource(std::string name = "") { return *(Resources.at(name).get()); }

    /// @}

    /// Destroys the specified resource (does nothing if no such resource)
    static Resource_t& DestroyResource(std::string name = "") { Resources.erase(name); }

  private:
    static std::map<std::string, ResourcePtr_t> Resources;

  }; // class TestSharedGlobalResource<>

  template <typename RES>
  std::map<std::string, typename TestSharedGlobalResource<RES>::ResourcePtr_t>
    TestSharedGlobalResource<RES>::Resources;

  /** **************************************************************************
   * @brief Environment for a test
   * @tparam ConfigurationClass a class providing compile-time configuration
   *
   * The test environment is set up on construction.
   *
   * The environment provides:
   * - Parameters() method returning the complete FHiCL configuration
   * - TesterParameters() method returning the configuration for the test
   *
   * This class or a derived one can be used as global fixture for unit tests.
   *
   * Unfortunately Boost does not give any control on the initialization of the
   * object, so everything must be ready to go as hard coded.
   * The ConfigurationClass class tries to alleviate that.
   * That is another, small static class that BasicTesterEnvironment uses to
   * get its parameters.
   *
   * The requirements for the ConfigurationClass are:
   * - `std::string ApplicationName()`: the application name
   * - `std::string ConfigurationPath()`: path to the configuration file
   * - `std::string MainTesterParameterSetName()`: name of the
   *   configuration of the main test (commodity)
   * - `std::string DefaultTesterConfiguration()` returning a FHiCL string
   *   to be parsed to extract the default test configuration
   *
   * Whether the configuration comes from a file or from the two provided
   * defaults, it is always expected within the parameter set paths:
   * the default configuration must also contain that path.
   *
   * Note that there is no room for polymorphism here since the setup happens
   * on construction.
   * Some methods are declared virtual in order to allow to tweak some steps
   * of the set up, but it's not trivial to create a derived class that works
   * correctly: the derived class must declare a new default constructor,
   * and that default constructor must call the protected constructor
   * (BasicTesterEnvironment<ConfigurationClass>(no_setup))
   */
  template <typename ConfigurationClass>
  class BasicTesterEnvironment {

  public:
    using Configuration_t = ConfigurationClass;

    /**
     * @brief Constructor: sets everything up and declares the test started
     * @param bSetup (_default: `true`_) call `Setup()` after construction
     *
     * The configuration is from a default-constructed ConfigurationClass.
     * This is suitable for use as Boost unit test fixture.
     */
    BasicTesterEnvironment(bool bSetup = true)
    {
      if (bSetup) Setup();
    }

    //@{
    /**
     * @brief Setup from a configuration
     * @param configurer an instance of `ConfigurationClass`
     * @param bSetup (_default: `true`_) call `Setup()` after construction
     *
     * The configuration is from the specified configurer class.
     *
     * This constructor allows to use a non-default-constructed configuration.
     * This can't be used (at best of my knowledge) when using this class as
     * Boost unit test fixture.
     *
     * In the r-value-reference constructor, the configurer is moved.
     */
    BasicTesterEnvironment(Configuration_t const& configurer, bool bSetup = true)
      : config(configurer)
    {
      if (bSetup) Setup();
    }
    BasicTesterEnvironment(Configuration_t&& configurer, bool bSetup = true) : config(configurer)
    {
      if (bSetup) Setup();
    }
    //@}

    /// Destructor: closing remarks
    virtual ~BasicTesterEnvironment();

    /// @{
    /// @name Configuration retrieval

    /// Returns the full configuration
    fhicl::ParameterSet const& Parameters() const { return params; }

    /// Returns the configuration of the specified service
    fhicl::ParameterSet ServiceParameters(std::string service_name) const
    {
      return params.get<fhicl::ParameterSet>(config.ServiceParameterSetPath(service_name));
    }

    /// Returns the configuration of the specified test
    fhicl::ParameterSet TesterParameters(std::string test_name) const
    {
      return params.get<fhicl::ParameterSet>(config.TesterParameterSetPath(test_name));
    }

    /// Returns the configuration of the main test (undefined if no main test)
    fhicl::ParameterSet TesterParameters() const
    {
      if (config.MainTesterParameterSetName().empty())
        return {};
      else
        return TesterParameters(config.MainTesterParameterSetName());
    }

    /// @}

    static fhicl::ParameterSet CompileParameterSet(std::string cfg);

  protected:
    /// Returns a read-only version of the configuration
    Configuration_t const& Config() const { return config; }

    /// The complete initialization, ran at construction by default
    virtual void Setup();

    /// Reads and translates the configuration
    virtual void Configure();

    /**
     * @brief Creates a full configuration for the test
     * @return a parameters set with the complete configuration
     */
    virtual fhicl::ParameterSet DefaultParameters() const
    {
      return CompileParameterSet(config.DefaultConfiguration());
    }

    //@{
    /// Sets up the message facility
    virtual void SetupMessageFacility(fhicl::ParameterSet const& pset,
                                      std::string appl_name = "") const;
    virtual void SetupMessageFacility() const
    {
      SetupMessageFacility(Parameters(), config.ApplicationName());
    }
    //@}

    /**
     * @brief Fills the test configuration from file or from default
     *
     * If a FHiCL configuration file is specified, the configuration of the test
     * is read from it according to the parameter set path of the test.
     * Otherwise, it is parsed from the default one provided by the configurer.
     */
    /// Parses from file and returns a FHiCL data structure
    static fhicl::ParameterSet ParseParameters(std::string config_path);

  private:
    /// Test environment options
    struct Options_t {
      bool MessageLevels = false; ///< print message levels on screen
    };                            // Options_t

    Configuration_t config; ///< instance of the configurer
    Options_t options;      ///< options for the test environment

    /// Parses the configuration, looking for the test environment options
    void ParseEnvironmentOptions();

    fhicl::ParameterSet params; ///< full configuration of the test

  }; // class BasicTesterEnvironment<>

  //****************************************************************************
  /**
   * @brief A test environment with some support for service providers
   * @tparam ConfigurationClass a class providing compile-time configuration
   *
   * This test environment extends BasicTesterEnvironment with some basic
   * support for service providers.
   *
   *
   * Service provider support
   * =========================
   *
   * This environment makes it available the method `Provider<Prov>()`, which
   * returns a pointer to the provider of type Prov.
   *
   * All providers must be set up _after_ the test environment is constructed.
   * The environment provides the following facilities:
   *
   * * SetupProvider() to set up a service provider with generic arguments
   * * SetupProviderFromService() to set up a service provider with a parameter
   *     set extracted from the configuration
   * * AcquireProvider() to register a service provider already available
   * * DropProvider() to destroy an existing provider
   *
   * The set up methods support a `For` variant (e.g. `SetupProviderFor()`) to
   * register the provider also under the type of its interface. For example,
   * if `LArPropertiesStandard` is an implementation of `LArProperties`,
   * the call:
   *
   *     env.SetupProviderFor<LArProperties, LArPropertiesStandard>(pset);
   *
   * will set up a `LArPropertiesStandard` provider just like
   *
   *     env.SetupProvider<LArPropertiesStandard>(pset);
   *
   * would, and it makes the provider available as `LArProperties` too, so that
   * both calls:
   *
   *     env.Provider<LArProperties>();
   *     env.Provider<LArPropertiesStandard>();
   *
   * are valid and return the same provider.
   *
   *
   * Use as test fixture
   * ====================
   *
   * The providers must be set up _after_ the test environment is constructed.
   * This also means an additional complication for fixtures that require to
   * be constructed in a final state, as it is the case for Boost unit test
   * suite fixtures.
   * In these cases, a class should publicly derive from TesterEnvironment, and
   * the necessary setup should be added into the constructor of this derived
   * class.
   *
   * Note that, as in the case of BasicTesterEnvironment, in this case there is
   * no room for polymorphism here since the setup need to happen on
   * construction.
   */
  template <typename ConfigurationClass>
  class TesterEnvironment : public BasicTesterEnvironment<ConfigurationClass> {
    using TesterEnvBase_t = BasicTesterEnvironment<ConfigurationClass>;
    using TesterEnv_t = TesterEnvironment<ConfigurationClass>;

  public:
    // inherit constructors
    using TesterEnvBase_t::TesterEnvBase_t;

    /**
     * @brief Sets a service provider up by calling its testing::setupProvider()
     * @tparam Prov type of provider
     * @tparam Args type of arguments for the setup function
     * @param args arguments for the setup function
     * @return a pointer to the provider set up
     * @throw runtime_error if the provider already exists
     * @see SetupProviderFor(), AcquireProvider()
     *
     * A provider of type Prov is created, set up and recorded.
     * Provider setup is delegated to `testing::setupProvider` function specific
     * to the provider itself (that is, `testing::setupProvider<Prov>(args...)`)
     * to which the setup arguments are forwarded.
     * If the provider already exists, an exception is thrown.
     */
    template <typename Prov, typename... Args>
    Prov* SetupProvider(Args&&... args)
    {
      if (!providers.setup<Prov>(std::forward<Args>(args)...))
        throw std::runtime_error("Provider already exists!");
      return providers.getPointer<Prov>();
    }

    /**
     * @brief Sets a service provider up by calling its testing::setupProvider()
     * @tparam Prov type of provider
     * @tparam Args type of arguments for the setup function
     * @param name the service name (for configuration retrieval)
     * @param args arguments for the setup function
     * @return a pointer to the provider set up
     * @see SetupProvider()
     * @throw runtime_error if the provider already exists
     *
     * A provider of type Prov is created, set up and recorded.
     * Provider setup is attempted by constructing the provider with a parameter
     * set from the registered configuration of service with specified `name`.
     */
    template <typename Prov, typename... Args>
    Prov* SetupProviderFromService(std::string name, Args&&... args)
    {
      return SetupProvider<Prov>(this->ServiceParameters(name, std::forward<Args>(args)...));
    }

    /**
     * @brief Acquires a service provider
     * @tparam Prov type of provider
     * @param prov the provider to be acquired
     * @return a pointer to the provider
     * @see SetupProvider()
     * @throw runtime_error if the provider already exists
     *
     * This method registers and takes ownership of the specified provider.
     * It is similar to SetupProvider() except that user is in charge of the
     * preliminary creation and setup of the provider.
     */
    template <typename Prov>
    Prov* AcquireProvider(std::unique_ptr<Prov>&& prov)
    {
      if (!providers.acquire(std::move(prov))) throw std::runtime_error("Provider already exists!");
      return providers.getPointer<Prov>();
    }

    /**
     * @brief Sets a provider up, recording it as implementation of Interface
     * @tparam Interface type of provider interface being implemented
     * @tparam Prov type of provider
     * @tparam Args type of arguments for the setup function
     * @param args arguments for the setup function
     * @return a pointer to the provider set up
     * @see SetupProvider()
     *
     * This method performs the same type of setup as SetupProvider().
     * In addition, it registers the provider as an implementation of Interface.
     * This means that the provider can be obtained not only with
     * `provider<Prov>()`, which returns a pointer to the actual class Prov,
     * but also as `provider<Interface>()`, which returns a pointer to the base
     * class Interface.
     */
    template <typename Interface, typename Prov, typename... Args>
    Prov* SetupProviderFor(Args&&... args)
    {
      auto prov = SetupProvider<Prov>(std::forward<Args>(args)...);
      providers.set_alias<Prov, Interface>();
      return prov;
    }

    /**
     * @brief Sets a provider up, recording it as implementation of Interface
     * @tparam Interface type of provider interface being implemented
     * @tparam Prov type of provider
     * @tparam Args type of arguments for the setup function
     * @param name the service name (for configuration retrieval)
     * @param args arguments for the setup function
     * @return a pointer to the provider set up
     * @see SetupProviderFromService(), SetupProviderFor()
     * @throw runtime_error if the provider already exists
     *
     * This method performs the same type of setup as
     * SetupProviderFromService().
     * In addition, it registers the provider as an implementation of Interface.
     * This means that the provider can be obtained not only with
     * `provider<Prov>()`, which returns a pointer to the actual class Prov,
     * but also as `provider<Interface>()`, which returns a pointer to the base
     * class Interface.
     */
    template <typename Interface, typename Prov, typename... Args>
    Prov* SetupProviderFromServiceFor(std::string name, Args&&... args)
    {
      auto* prov = SetupProviderFromService<Prov>(name, std::forward<Args>(args)...);
      providers.set_alias<Prov, Interface>();
      return prov;
    }

    /**
     * @brief Acquires a service provider implementing an interface
     * @tparam Prov type of provider
     * @tparam Interface type provider alias
     * @param prov the provider to be acquired
     * @return a pointer to the provider
     * @see SetupProviderFor(), AcquireProvider()
     *
     * This method registers and takes ownership of the specified provider,
     * like AcquireProvider() does. It also registers the provider as an
     * implementation of Interface class, as SetupProviderFor does.
     * It is similar to SetupProvider() except that user is in charge of the
     * preliminar creation and setup of the provider.
     */
    template <typename Interface, typename Prov>
    Prov* AcquireProviderFor(std::unique_ptr<Prov>&& prov)
    {
      auto prov_ptr = providers.acquire(prov);
      providers.set_alias<Prov, Interface>();
      return prov_ptr;
    }

    /**
     * @brief Oversimplified provider setup
     * @return a pointer to the provider
     * @tparam Prov provider type
     *
     * This is a one-step setup of the specified provider.
     *
     * It is available only if Prov provider comes with an implementation of
     * testing::SimpleEnvironmentSetupClass that explains how to set up an
     * environment.
     *
     */
    template <typename Prov>
    Prov* SimpleProviderSetup()
    {
      return simpleEnvironmentSetup<Prov>(*this);
    }

    /**
     * @brief Removes and destroys the specified provider
     * @tparam Prov type of provider to be destroyed
     * @throw runtime_error if the provider was not present
     */
    template <typename Prov>
    void DropProvider()
    {
      if (!providers.erase<Prov>()) throw std::runtime_error("Provider not present!");
    }

    /// Return the specified provider (throws if not available)
    template <typename Prov>
    Prov const* Provider() const
    {
      return providers.getPointer<Prov>();
    }

    /**
     * @brief Fills the specified provider pack with providers
     * @throw runtime_error and everything provider() method can throw
     * @see Provider()
     */
    template <typename... Provs>
    void FillProviderPack(lar::ProviderPack<Provs...>& pack) const
    {
      details::ProviderPackFiller<TesterEnv_t, lar::ProviderPack<Provs...>, Provs...>::fill(*this,
                                                                                            pack);
    } // FillProviderPack()

    /**
     * @brief Returns a provider pack for the specified provider
     * @tparam Prov type of the provider
     * @throw runtime_error and everything provider() method can throw
     * @see FillProviderPack()
     *
     * The provider is required to have a `providers_type` type defined as an
     * specialisation of lar::ProviderPack.
     */
    template <typename Prov>
    typename Prov::providers_type ProviderPackFor() const
    {
      typename Prov::providers_type pack;
      FillProviderPack(pack);
      return pack;
    } // ProviderPackFor()

  protected:
    ProviderList providers; ///< list of available providers
  };                        // class TesterEnvironment<>

  /**
   * @brief Constructs and returns a TesterEnvironment object
   * @tparam CONFIG type of configuration object (detected from arguments)
   * @tparam TESTENV type of the tester environment (default: TesterEnvironment)
   * @tparam ARGS pack of types of the remining constructor arguments (optional)
   * @param config the configuration object to be used
   * @param other_args the remaining arguments of the tester constructor
   *
   * This function creates and returns a tester environment.
   * By default, the tester environment class is TesterEnvironment<CONFIG>
   * and no additional constructor arguments are needed except for special
   * needs. The simplest way to use the function with an already available
   * configuration is:
   *
   *    auto TestEnv = testing::CreateTesterEnvironment(config);
   *
   * where TestEnv is assigned a specialization of TesterEnvironment.
   *
   * The template class TESTENV undergoes the following requirement:
   *
   *  - it must have a constructor using a CONFIG constant reference as first
   *    argument
   *
   * The CONFIG object is subject to no special requirements besides the ones
   * from TESTENV constructor.
   */
  // Note: this function is expected to be used with automatic type detection;
  // the rules on "universal references" dictate that if config is a (l-value)
  // reference, CONFIG itself is a l-value reference. We don't want to create
  // a TesterEnvironment<Config&>, so we explicitly remove the reference from
  // CONFIG (decay does that and aso removes the constantness, that we also
  // don't want to be embedded in CONFIG).
  template <typename CONFIG,
            typename TESTENV = TesterEnvironment<std::decay_t<CONFIG>>,
            typename... ARGS>
  TESTENV CreateTesterEnvironment(CONFIG&& config, ARGS... other_args)
  {
    return TESTENV(std::forward<CONFIG>(config), std::forward<ARGS>(other_args)...);
  }

  //****************************************************************************
  namespace details {
    /// Helper to fill a provider pack: main specialisation
    template <typename TestEnv, typename Pack, typename Prov, typename... Others>
    struct ProviderPackFiller<TestEnv, Pack, Prov, Others...> {
      static void fill(TestEnv const& env, Pack& pack)
      {
        pack.set(env.template Provider<Prov>());
        ProviderPackFiller<TestEnv, Pack, Others...>::fill(env, pack);
      } // fill()

    }; // ProviderPackFiller<TestEnv, Pack, Prov, Others...>

    // end-of-recursion specialisation
    template <typename TestEnv, typename Pack>
    struct ProviderPackFiller<TestEnv, Pack> {
      static void fill(TestEnv const&, Pack&) {}
    }; // ProviderPackFiller<>

  } // namespace details

  //****************************************************************************
  template <typename ConfigurationClass>
  BasicTesterEnvironment<ConfigurationClass>::~BasicTesterEnvironment()
  {
    mf::LogInfo("Test") << config.ApplicationName() << " completed.";
  }

  /** **************************************************************************
   * @brief Compiles a parameter set from a string
   * @return a parameters set with the complete configuration
   */
  template <typename ConfigurationClass>
  fhicl::ParameterSet BasicTesterEnvironment<ConfigurationClass>::CompileParameterSet(
    std::string cfg)
  {
    fhicl::ParameterSet global_pset;
    global_pset = fhicl::ParameterSet::make(cfg);
    return global_pset;
  }

  /** **************************************************************************
   * @brief Returns the configuration from a FHiCL file
   * @param config_path full path of the FHiCL configuration file
   * @return a parameters set with the complete configuration from the file
   */
  template <typename ConfigurationClass>
  fhicl::ParameterSet BasicTesterEnvironment<ConfigurationClass>::ParseParameters(
    std::string config_path)
  {
    // configuration file lookup policy
    char const* fhicl_env = getenv("FHICL_FILE_PATH");
    std::string search_path = fhicl_env ? std::string(fhicl_env) + ":" : ".:";
    cet::filepath_first_absolute_or_lookup_with_dot policy(search_path);

    // parse a configuration file; obtain intermediate form
    auto table = fhicl::parse_document(config_path, policy);
    return fhicl::ParameterSet::make(table);
  }

  /** **************************************************************************
   * @brief Fills the configuration
   *
   * The complete configuration (message facility and services) is read and
   * saved, hence accessible by Parameters() method.
   *
   * The configuration file path is taken by default from the first argument
   * of the test.
   * If that first argument is not present or empty, the default configuration
   * path is received from the configurer.
   * If the configuration path is still empty, a hard-coded configuration
   * is used; otherwise, the FHiCL file specified in that path is parsed and
   * used as full configuration.
   */
  template <typename ConfigurationClass>
  void BasicTesterEnvironment<ConfigurationClass>::Configure()
  {
    std::string config_path = config.ConfigurationPath();
    params = config_path.empty() ? DefaultParameters() : ParseParameters(config_path);
  }

  /** **************************************************************************
   * @brief Sets the message facility up
   *
   * Message facility configuration is expected in "services.message" parameter
   * set. If not there, the default configuration is used.
   */
  template <typename ConfigurationClass>
  void BasicTesterEnvironment<ConfigurationClass>::SetupMessageFacility(
    fhicl::ParameterSet const& pset,
    std::string appl_name /* = "" */) const
  {
    fhicl::ParameterSet mf_pset;

    if (!pset.get_if_present(config.ServiceParameterSetPath("message"), mf_pset)) {
      mf_pset = CompileParameterSet(config.DefaultServiceConfiguration("message"));
      std::cout << "Using default message facility configuration:\n"
                << mf_pset.to_indented_string(1) << std::endl;
    } // if no configuration is available

    mf::StartMessageFacility(mf_pset);
    if (!appl_name.empty()) mf::SetApplicationName(appl_name);
    mf::SetContextIteration("Initialization");
    if (options.MessageLevels) {
      std::cout << "Printing message levels in 'MessageFacility' category." << std::endl;

      mf::LogProblem("MessageFacility") << "Error messages are shown.";
      mf::LogPrint("MessageFacility") << "Warning messages are shown.";
      mf::LogVerbatim("MessageFacility") << "Info messages are shown.";
      mf::LogTrace("MessageFacility") << "Debug messages are shown.";
      MF_LOG_TRACE("MessageFacility")
        << "MF_LOG_TRACE/MF_LOG_DEBUG messages are not compiled away.";
    } // if print message levels
    mf::LogInfo("MessageFacility") << "MessageFacility started.";
    mf::SetContextSinglet("main");
  }

  template <typename ConfigurationClass>
  void BasicTesterEnvironment<ConfigurationClass>::Setup()
  {
    Configure();
    ParseEnvironmentOptions();
    SetupMessageFacility();

    // Optionally print the configuration
    {
      mf::LogInfo msg("Configuration");
      msg << "Complete configuration (";
      if (config.ConfigurationPath().empty())
        msg << "default";
      else
        msg << "'" << config.ConfigurationPath() << "'";
      msg << "):\n" << Parameters().to_indented_string(1);
    }
    mf::LogInfo("Test") << config.ApplicationName() << " base setup complete.";
  }

  template <typename ConfigurationClass>
  void BasicTesterEnvironment<ConfigurationClass>::ParseEnvironmentOptions()
  {

    struct OptionsFromConfig_t {
      fhicl::Atom<bool> messageLevels{
        fhicl::Name("messageLevels"),
        fhicl::Comment("prints a message per level (to verify the visible ones"),
        false // default: no
      };
      fhicl::Atom<bool> printOptions{
        fhicl::Name("printOptions"),
        fhicl::Comment("prints a the list of options (this is one of them!)"),
        false // default: no
      };
    }; // OptionsFromConfig_t

    struct ValidationHelper {
      static void printDummy(std::ostream& out)
      {
        fhicl::Table<OptionsFromConfig_t>(fhicl::Name("test"),
                                          fhicl::Comment("Test environment options"))
          .print_allowed_configuration(out);
      } // printDummy()

      static fhicl::Table<OptionsFromConfig_t> validate(fhicl::ParameterSet const& pset)
      {
        try {
          return fhicl::Table<OptionsFromConfig_t>(pset, {});
        }
        catch (...) {
          std::cerr << "Error parsing environment test options! Valid options:" << std::endl;
          ValidationHelper::printDummy(std::cerr);
          throw;
        }
      } // validate()
    };

    fhicl::ParameterSet pset = params.get<fhicl::ParameterSet>("test", {});

    // automatically performs validation
    fhicl::Table<OptionsFromConfig_t> configTable(ValidationHelper::validate(pset));

    if (configTable().printOptions()) {
      std::cout << "The following options can be passed to the test environment"
                << " by putting them in the \"test: {}\" table of the configuration file:"
                << std::endl;
      ValidationHelper::printDummy(std::cout);
    }

    options.MessageLevels = configTable().messageLevels();
  } // BasicTesterEnvironment<>::ParseEnvironmentOptions()

} // namespace testing

#endif // TEST_UNIT_TEST_BASE_H
