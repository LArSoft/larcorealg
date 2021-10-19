/**
 * @file   UncopiableAndUnmovableClass_test.cc
 * @brief  Tests the content of UncopiableAndUnmovableClass.h
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   April 28, 2016
 * @see    larcorealg/CoreUtils/UncopiableAndUnmovableClass.h
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
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

// Boost libraries
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <type_traits>


//------------------------------------------------------------------------------
// static checks
static_assert( std::is_copy_constructible_v<lar::PolymorphicClass>);
static_assert( std::is_copy_assignable_v   <lar::PolymorphicClass>);
static_assert( std::is_move_constructible_v<lar::PolymorphicClass>);
static_assert( std::is_move_assignable_v   <lar::PolymorphicClass>);

static_assert(!std::is_copy_constructible_v<lar::PolymorphicUncopiableClass>);
static_assert(!std::is_copy_assignable_v   <lar::PolymorphicUncopiableClass>);
static_assert( std::is_move_constructible_v<lar::PolymorphicUncopiableClass>);
static_assert( std::is_move_assignable_v   <lar::PolymorphicUncopiableClass>);

static_assert( std::is_copy_constructible_v<lar::PolymorphicUnmovableClass>);
static_assert( std::is_copy_assignable_v   <lar::PolymorphicUnmovableClass>);
static_assert(!std::is_move_constructible_v<lar::PolymorphicUnmovableClass>);
static_assert(!std::is_move_assignable_v   <lar::PolymorphicUnmovableClass>);

static_assert(!std::is_copy_constructible_v<lar::PolymorphicUncopiableAndUnmovableClass>);
static_assert(!std::is_copy_assignable_v   <lar::PolymorphicUncopiableAndUnmovableClass>);
static_assert(!std::is_move_constructible_v<lar::PolymorphicUncopiableAndUnmovableClass>);
static_assert(!std::is_move_assignable_v   <lar::PolymorphicUncopiableAndUnmovableClass>);


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(UncopiableAndUnmovableClassTest) {

   // check lar::UncopiableAndUnmovableClass class itself
   BOOST_TEST
     (!std::is_copy_constructible<lar::UncopiableAndUnmovableClass>::value);
   BOOST_TEST
     (!std::is_copy_assignable<lar::UncopiableAndUnmovableClass>::value);
   BOOST_TEST
     (!std::is_move_constructible<lar::UncopiableAndUnmovableClass>::value);
   BOOST_TEST
     (!std::is_move_assignable<lar::UncopiableAndUnmovableClass>::value);


   // check a class derived from lar::UncopiableAndUnmovableClass class
   struct Derived: protected lar::UncopiableAndUnmovableClass {};

   BOOST_TEST(!std::is_copy_constructible<Derived>::value);
   BOOST_TEST(!std::is_copy_assignable   <Derived>::value);
   BOOST_TEST(!std::is_move_constructible<Derived>::value);
   BOOST_TEST(!std::is_move_assignable   <Derived>::value);


   // check a class derived from lar::UncopiableAndUnmovableClass class
   // and made movable
   struct MovableDerived: protected lar::UncopiableAndUnmovableClass {
      MovableDerived(MovableDerived&&):
        lar::UncopiableAndUnmovableClass() {}
   };

   BOOST_TEST(!std::is_copy_constructible<MovableDerived>::value);
   BOOST_TEST(!std::is_copy_assignable   <MovableDerived>::value);
   BOOST_TEST( std::is_move_constructible<MovableDerived>::value);
   BOOST_TEST(!std::is_move_assignable   <MovableDerived>::value);


   // check a class derived from lar::UncopiableAndUnmovableClass class
   // and made both copy- and move-assignable
   struct AssignableDerived: protected lar::UncopiableAndUnmovableClass {
      AssignableDerived& operator=(AssignableDerived const&) { return *this; }
      AssignableDerived& operator=(AssignableDerived&&) { return *this; }
   };

   BOOST_TEST(!std::is_copy_constructible<AssignableDerived>::value);
   BOOST_TEST( std::is_copy_assignable   <AssignableDerived>::value);
   BOOST_TEST(!std::is_move_constructible<AssignableDerived>::value);
   BOOST_TEST( std::is_move_assignable   <AssignableDerived>::value);

} // BOOST_AUTO_TEST_CASE(UncopiableAndUnmovableClassTest)


//------------------------------------------------------------------------------
