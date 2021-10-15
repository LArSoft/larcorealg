/**
 * @file   boost_unit_test_base.h
 * @brief  Base class for tests using Boost unit test library
 * @date   May 22th, 2015
 * @author petrillo@fnal.gov
 * @see    unit_test_base.h
 *
 * Provides an environment for easy set up of a Boost test.
 * This is based and wraps the objects in unit_test_base.h.
 * Since the wrapping is implemented by templates, derived classes than share
 * interface with the unit_test_base.h objects should also work with the
 * wrappers.
 *
 * For an example of usage, see larcore/test/Geometry/geometry_iterator_test.cxx
 *
 * This is a pure template header. It will require the same libraries as
 * unit_test_base.h .
 *
 */


#ifndef TEST_BOOST_UNIT_TEST_BASE_H
#define TEST_BOOST_UNIT_TEST_BASE_H

// LArSoft libraries
#include "larcorealg/TestUtils/unit_test_base.h"

// Boost libraries
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <string>


namespace testing {

  /** **************************************************************************
   * @brief Class holding a configuration for a Boost test fixture
   * @tparam CONFIGURATIONCLASS a base configuration class
   * @see BasicEnvironmentConfiguration, TesterEnvironment
   *
   * This class needs to be fully constructed by the default constructor
   * in order to be useful as Boost unit test fixture.
   * It is supposed to be passed as a template parameter to another class
   * that can store an instance of it and extract configuration information
   * from it.
   *
   * This template just adds to the standard construction of the wrapped class
   * a configuration that reads the parameters from the command line.
   * It also hides all the constructors except two.
   */
  template <typename CONFIGURATIONCLASS>
  struct BoostCommandLineConfiguration: public CONFIGURATIONCLASS {

    using Base_t = CONFIGURATIONCLASS;

    /// Default constructor; this is what is used in Boost unit test
    BoostCommandLineConfiguration(): Base_t()
      { ParseCommandLineFromBoost(); }

    /// Constructor; accepts the name as parameter
    BoostCommandLineConfiguration(std::string name): Base_t(name)
      { ParseCommandLineFromBoost(); }

      protected:

    /// Parses arguments as delivered by Boost
    void ParseCommandLineFromBoost()
      {
        Base_t::ParseCommandLine(
          boost::unit_test::framework::master_test_suite().argc,
          boost::unit_test::framework::master_test_suite().argv
          );
      }

  }; // class BoostCommandLineConfiguration<>


} // namespace testing

#endif // TEST_BOOST_UNIT_TEST_BASE_H
