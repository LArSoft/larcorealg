/**
 * @file   UncopiableAndUnmovableClass_test.cc
 * @brief  Tests the content of UncopiableAndUnmovableClass.h
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   April 28, 2016
 * @see    UncopiableAndUnmovableClass.h
 * 
 * This test takes no command line argument.
 * 
 */

/*
 * Boost Magic: define the name of the module;
 * and do that before the inclusion of Boost unit test headers
 * because it will change what they provide.
 * Among the those, there is a main() function and some wrapping catching
 * unhandled exceptions and considering them test failures, and probably more.
 */
#define BOOST_TEST_MODULE ( UncopiableAndUnmovableClass_test )

// LArSoft libraries
#include "larcore/CoreUtils/UncopiableAndUnmovableClass.h"

// Boost libraries
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// C/C++ standard libraries
#include <type_traits>


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(UncopiableAndUnmovableClassTest) {
   
   // check lar::UncopiableAndUnmovableClass class itself
   BOOST_CHECK
     (!std::is_copy_constructible<lar::UncopiableAndUnmovableClass>::value);
   BOOST_CHECK
     (!std::is_copy_assignable<lar::UncopiableAndUnmovableClass>::value);
   BOOST_CHECK
     (!std::is_move_constructible<lar::UncopiableAndUnmovableClass>::value);
   BOOST_CHECK
     (!std::is_move_assignable<lar::UncopiableAndUnmovableClass>::value);
   
   
   // check a class derived from lar::UncopiableAndUnmovableClass class
   struct Derived: protected lar::UncopiableAndUnmovableClass {};
   
   BOOST_CHECK(!std::is_copy_constructible<Derived>::value);
   BOOST_CHECK(!std::is_copy_assignable   <Derived>::value);
   BOOST_CHECK(!std::is_move_constructible<Derived>::value);
   BOOST_CHECK(!std::is_move_assignable   <Derived>::value);
   
   
   // check a class derived from lar::UncopiableAndUnmovableClass class
   // and made movable
   struct MovableDerived: protected lar::UncopiableAndUnmovableClass {
      MovableDerived(MovableDerived&&):
        lar::UncopiableAndUnmovableClass() {}
   };
   
   BOOST_CHECK(!std::is_copy_constructible<MovableDerived>::value);
   BOOST_CHECK(!std::is_copy_assignable   <MovableDerived>::value);
   BOOST_CHECK( std::is_move_constructible<MovableDerived>::value);
   BOOST_CHECK(!std::is_move_assignable   <MovableDerived>::value);
   
   
   // check a class derived from lar::UncopiableAndUnmovableClass class
   // and made both copy- and move-assignable
   struct AssignableDerived: protected lar::UncopiableAndUnmovableClass {
      AssignableDerived& operator=(AssignableDerived const&) { return *this; }
      AssignableDerived& operator=(AssignableDerived&&) { return *this; }
   };
   
   BOOST_CHECK(!std::is_copy_constructible<AssignableDerived>::value);
   BOOST_CHECK( std::is_copy_assignable   <AssignableDerived>::value);
   BOOST_CHECK(!std::is_move_constructible<AssignableDerived>::value);
   BOOST_CHECK( std::is_move_assignable   <AssignableDerived>::value);
   
} // BOOST_AUTO_TEST_CASE(larProviderFromTest)


//------------------------------------------------------------------------------
