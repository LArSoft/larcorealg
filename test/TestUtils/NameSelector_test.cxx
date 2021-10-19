/**
 * @file   NameSelector_test.cxx
 * @brief  Test of NameSelector class
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 24th, 2015
 */

// LArSoft libraries
#include "larcorealg/TestUtils/NameSelector.h"

// Boost libraries
#include "canvas/Utilities/Exception.h"

// Boost libraries
#define BOOST_TEST_MODULE ( NameSelector_test )
#include <boost/test/unit_test.hpp>

// C/C++ standard library
#include <string>
#include <iostream>

/**
 * In the following tests, the names are conventionally called "acceptXXX" if
 * the final result of the query is expected to be acceptance, "rejectXXX" if
 * the final result of the query is expected to be rejection, or "throwXXX" if
 * the final result of the query is expected to be the throw of an exception.
 *
 * If an error does happen, run the test with '--log_level=all' to get enough
 * information to track the error location:
 *
 * - for an unexpected exception (even without '--log_level=all'):
 *
 *   .../NameSelector_test.cxx(53): error in "DefaultThrowTest": check selector.Query(name) == testing::NameSelector::rsAccepted failed [2 != 1]
 *   unknown location(0): fatal error in "DefaultThrowTest": std::exception: ---- Configuration BEGIN
 *     NameSelector: name 'accept_unknown' not configured.
 *   ---- Configuration END
 *
 *   Boost error points to the test name ("DefaultThrowTest") and to the check
 *   that failed (line 53), while the message of the exception fills in with
 *   the name being checked
 *
 * - for a wrong result, running with '--log_level=all':
 *
 *   Testing 'accept_unknown'
 *   ../NameSelector_test.cxx(53): error in "DefaultRejectTest": check selector.Query(name) == testing::NameSelector::rsAccepted failed [0 != 1]
 *
 *   Boost points to the test and the failed check as before, while
 *   BOOST_TEST_MESSAGE informs about which name was being tested.
 *
 */
void CheckNames
  (testing::NameSelector const& selector, std::vector<std::string> const& names)
{

  for (std::string const& name: names) {
    BOOST_TEST_MESSAGE("Testing '" << name << "'");
    if (name.find("throw") == 0U) {
      BOOST_TEST(selector.Query(name) == testing::NameSelector::rsThrow);
      BOOST_CHECK_THROW(selector(name), art::Exception);
      BOOST_CHECK_THROW(selector.Accepted(name), art::Exception);
      BOOST_CHECK_THROW(selector.Rejected(name), art::Exception);
    }
    else if (name.find("reject") == 0U) {
      BOOST_TEST
        (selector.Query(name) == testing::NameSelector::rsRejected);
      BOOST_TEST(!selector(name));
      BOOST_TEST(!selector.Accepted(name));
      BOOST_TEST(selector.Rejected(name));
    }
    else { // accept
      BOOST_TEST
        (selector.Query(name) == testing::NameSelector::rsAccepted);
      BOOST_TEST(selector(name));
      BOOST_TEST(selector.Accepted(name));
      BOOST_TEST(!selector.Rejected(name));
    }
  } // for

} // CheckNames()

//------------------------------------------------------------------------------
//
// Simple use case
//
BOOST_AUTO_TEST_CASE(SimpleTest) {

  // default answer: yes
  testing::NameSelector selector;

  //
  // initialize with simple elements
  //
  std::cout << "SimpleTest" << std::endl;
  selector.ParseNames("accept_one", "+accept_two");
  selector.ParseName("-reject_three");

  //
  // visually verify the configuration
  //
  selector.PrintConfiguration(std::cout);
  std::cout << std::endl;

  //
  // check all the symbols, plus one
  //
  CheckNames(selector, {
    "accept_one", "accept_two", "reject_three", "accept_unknown"
    });

  //
  // verify that we did not miss anything
  //
  BOOST_TEST(selector.CheckQueryRegistry(std::cout));

} // BOOST_AUTO_TEST_CASE(SimpleTest)


//
// Test for check of missing queries
//
BOOST_AUTO_TEST_CASE(MissingQuery) {

  // default answer: yes
  testing::NameSelector selector;

  //
  // initialize with simple elements
  //
  std::cout << "MissingQuery" << std::endl;
  selector.ParseName("accept_one");

  //
  // verify that we did miss something
  //
  BOOST_TEST(!selector.CheckQueryRegistry(std::cout));

} // BOOST_AUTO_TEST_CASE(MissingQuery)


//
// Test for default throw
//
BOOST_AUTO_TEST_CASE(DefaultConstructorTest) {

  // default answer: yes
  testing::NameSelector selector;

  //
  // initialize with simple elements
  //
  std::cout << "DefaultConstructorTest" << std::endl;

  //
  // visually verify the configuration
  //
  selector.PrintConfiguration(std::cout);
  std::cout << std::endl;

  //
  // check all the symbols
  //
  CheckNames(selector, { "accept_one", "accept_two" });

} // BOOST_AUTO_TEST_CASE(DefaultConstructorTest)


//
// Test for default throw
//
BOOST_AUTO_TEST_CASE(DefaultThrowTest) {

  // default answer: yes
  testing::NameSelector selector(testing::NameSelector::rsThrow);

  //
  // initialize with simple elements
  //
  std::cout << "DefaultThrowTest" << std::endl;
  selector.ParseNames("accept_one", "-reject_three");

  //
  // check all the symbols
  //
  CheckNames(selector, { "accept_one", "reject_three", "throw_two" });

} // BOOST_AUTO_TEST_CASE(DefaultThrowTest)


//
// Test for default rejection
//
BOOST_AUTO_TEST_CASE(DefaultRejectTest) {

  // default answer: yes
  testing::NameSelector selector(testing::NameSelector::rsRejected);

  //
  // initialize with simple elements
  //
  std::cout << "DefaultRejectTest" << std::endl;
  selector.ParseNames("accept_one", "-reject_three");

  //
  // check all the symbols, plus one
  //
  CheckNames(selector, { "accept_one", "reject_two", "reject_unknown" });

} // BOOST_AUTO_TEST_CASE(DefaultRejectTest)


//
// Test for default acceptance
//
BOOST_AUTO_TEST_CASE(DefaultAcceptTest) {

  // default answer: yes
  testing::NameSelector selector(testing::NameSelector::rsAccepted);

  //
  // initialize with simple elements
  //
  std::cout << "DefaultAcceptTest" << std::endl;
  selector.ParseNames("accept_one", "-reject_three");

  //
  // check all the symbols, plus one
  //
  CheckNames(selector, { "accept_one", "accept_two", "accept_unknown" });

} // BOOST_AUTO_TEST_CASE(DefaultAcceptTest)


//
// Test for adding name set
//
BOOST_AUTO_TEST_CASE(SetDefinitionTest) {

  // default answer: yes
  testing::NameSelector selector;

  // set 2 does not follow the naming convention because it will be subtracted
  selector.AddToDefinition
    ("set1", "accept_set1_one", "+accept_set1_two", "-reject_set1_three");
  selector.AddToDefinition
    ("set2", "-accept_set2_one", "-accept_set2_two", "+reject_set2_three");
  selector.Define("set3", std::vector<std::string>
    { "accept_set3_one", "+accept_set3_two", "-reject_set3_three" }
    );

  //
  // initialize with simple elements
  //
  std::cout << "SetDefinitionTest" << std::endl;
  selector.ParseNames("accept_one", "set1", "-@set2", "+set3", "-reject_two");

  //
  // visually verify the configuration
  //
  selector.PrintConfiguration(std::cout);
  std::cout << std::endl;

  //
  // check all the symbols, plus one
  //
  CheckNames(selector, {
    "accept_set1_one", "accept_set1_two", "reject_set1_three",
    "accept_set2_one", "accept_set2_two", "reject_set2_three",
    "accept_one", "reject_two", "accept_unknown"
    });

} // BOOST_AUTO_TEST_CASE(SetDefinitionTest)


//
// Test for double specification
//
BOOST_AUTO_TEST_CASE(OverrideTest) {

  // default answer: yes
  testing::NameSelector selector;

  // set does not follow the naming convention because it will be subtracted
  selector.AddToDefinition("default",
    "accept_one", "-accept_two", "+reject_three"
    );

  //
  // initialize with simple elements
  //
  std::cout << "OverrideTest" << std::endl;
  selector.ParseNames("+reject_three", "-default", "accept_one", "accept_four");

  //
  // visually verify the configuration
  //
  selector.PrintConfiguration(std::cout);
  std::cout << std::endl;

  //
  // check all the symbols
  //
  CheckNames(selector, {
    "accept_one", "accept_two", "reject_three", "accept_four"
    });

} // BOOST_AUTO_TEST_CASE(OverrideTest)


//
// Global specifications
//
BOOST_AUTO_TEST_CASE(GlobalSpecTest) {

  // default answer: yes
  testing::NameSelector selector;

  // set 2 does not follow the naming convention because it will be subtracted
  selector.AddToDefinition("default",
    "accept_one", "-accept_two", "+reject_three"
    );

  //
  // initialize with simple elements
  //
  std::cout << "GlobalSpecTest" << std::endl;
  selector.ParseNames("accept_one", "-*", "-reject_three", "accept_four");

  //
  // visually verify the configuration
  //
  selector.PrintConfiguration(std::cout);
  std::cout << std::endl;

  //
  // check all the symbols, plus one
  //
  CheckNames(selector, {
    "accept_one", "reject_unknown", "reject_three", "accept_four"
    });

} // BOOST_AUTO_TEST_CASE(GlobalSpecTest)


//
// Clear-all directive
//
BOOST_AUTO_TEST_CASE(ClearAllTest) {

  // default answer: yes
  testing::NameSelector selector(testing::NameSelector::rsThrow);

  //
  // initialize
  //
  std::cout << "ClearAllTest" << std::endl;
  selector.ParseNames("throw_lost", "!", "-reject_three", "accept_four");

  //
  // visually verify the configuration
  //
  selector.PrintConfiguration(std::cout);
  std::cout << std::endl;

  //
  // check all the symbols, plus one
  //
  CheckNames(selector, {
    "throw_lost", "throw_unknown", "reject_three", "accept_four"
    });

} // BOOST_AUTO_TEST_CASE(ClearAllTest)
