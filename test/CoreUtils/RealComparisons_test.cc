/**
 * @file    RealComparisons_test.cc
 * @brief   Unit test for `RealComparisons` class
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    November 4th, 2016
 * @see     larcorealg/CoreUtils/RealComparisons.h
 */

// Boost libraries
#define BOOST_TEST_MODULE ( RealComparisons_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_TEST()

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"

// C/C++ standard libraries
#include <ostream>
#include <cmath>


// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_RealComparisons) {

  constexpr lar::util::RealComparisons check(1e-5f);

  constexpr double const sqrt2 = std::sqrt(2);

  // check zero()
  constexpr double const epsilon = 2.0 - (sqrt2 * sqrt2);
  BOOST_TEST( check.zero(epsilon));
  BOOST_TEST( check.zero(0.0));
  BOOST_TEST( check.zero(1e-5));
  BOOST_TEST( check.zero(-1e-5));
  BOOST_TEST(!check.zero(1.01e-5));
  BOOST_TEST(!check.zero(-1.01e-5 * 1.01));

  // check nonZero()
  BOOST_TEST(!check.nonZero(epsilon));
  BOOST_TEST(!check.nonZero(0.0));
  BOOST_TEST(!check.nonZero(1e-5));
  BOOST_TEST(!check.nonZero(-1e-5));
  BOOST_TEST( check.nonZero(1.01e-5));
  BOOST_TEST( check.nonZero(-1.01e-5 * 1.01));

  // check equal()
  BOOST_TEST(!check.equal(sqrt2, 1.4142));
  BOOST_TEST( check.nonEqual(sqrt2, 1.4142));
  BOOST_TEST( check.equal(sqrt2, 1.414213));
  BOOST_TEST(!check.nonEqual(sqrt2, 1.414213));

  // check strictlyNegative()
  BOOST_TEST(!check.strictlyNegative(+1e-5 + 1e-7)); // outside tolerance
  BOOST_TEST(!check.strictlyNegative(+1e-5 - 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyNegative(0.0));          // zero
  BOOST_TEST(!check.strictlyNegative(-1e-5 + 1e-7)); // within tolerance
  BOOST_TEST( check.strictlyNegative(-1e-5 - 1e-7)); // outside tolerance

  // check strictlyPositive()
  BOOST_TEST( check.strictlyPositive(+1e-5 + 1e-7)); // outside tolerance
  BOOST_TEST(!check.strictlyPositive(+1e-5 - 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyPositive(0.0));          // zero
  BOOST_TEST(!check.strictlyPositive(-1e-5 + 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyPositive(-1e-5 - 1e-7)); // outside tolerance

  // check nonNegative()
  BOOST_TEST( check.nonNegative(+1e-5 + 1e-7)); // outside tolerance
  BOOST_TEST( check.nonNegative(+1e-5 - 1e-7)); // within tolerance
  BOOST_TEST( check.nonNegative(0.0));          // zero
  BOOST_TEST( check.nonNegative(-1e-5 + 1e-7)); // within tolerance
  BOOST_TEST(!check.nonNegative(-1e-5 - 1e-7)); // outside tolerance

  // check nonPositive()
  BOOST_TEST(!check.nonPositive(+1e-5 + 1e-7)); // outside tolerance
  BOOST_TEST( check.nonPositive(+1e-5 - 1e-7)); // within tolerance
  BOOST_TEST( check.nonPositive(0.0));          // zero
  BOOST_TEST( check.nonPositive(-1e-5 + 1e-7)); // within tolerance
  BOOST_TEST( check.nonPositive(-1e-5 - 1e-7)); // outside tolerance

  // check strictlySmaller()
  BOOST_TEST(!check.strictlySmaller(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_TEST(!check.strictlySmaller(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlySmaller(1.0, 1.0));        // equal
  BOOST_TEST(!check.strictlySmaller(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_TEST( check.strictlySmaller(1.0, 1.0 + 1e-4)); // outside tolerance

  // check nonSmaller()
  BOOST_TEST( check.nonSmaller(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_TEST( check.nonSmaller(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_TEST( check.nonSmaller(1.0, 1.0));        // equal
  BOOST_TEST( check.nonSmaller(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_TEST(!check.nonSmaller(1.0, 1.0 + 1e-4)); // outside tolerance

  // check strictlyGreater()
  BOOST_TEST( check.strictlyGreater(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_TEST(!check.strictlyGreater(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyGreater(1.0, 1.0));        // equal
  BOOST_TEST(!check.strictlyGreater(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_TEST(!check.strictlyGreater(1.0, 1.0 + 1e-4)); // outside tolerance

  // check nonGreater()
  BOOST_TEST(!check.nonGreater(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_TEST( check.nonGreater(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_TEST( check.nonGreater(1.0, 1.0));        // equal
  BOOST_TEST( check.nonGreater(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_TEST( check.nonGreater(1.0, 1.0 + 1e-4)); // outside tolerance

  // check within()
  BOOST_TEST(!check.within(sqrt2, 0., 1.41420));
  BOOST_TEST( check.within(sqrt2, 0., 1.41421));
  BOOST_TEST( check.within(sqrt2, 1.41422, 2.));
  BOOST_TEST(!check.within(sqrt2, 1.41423, 2.));

  // check inverted limits
  BOOST_TEST(!check.within(sqrt2, 1.41421, 0.));
  BOOST_TEST(check.withinSorted(sqrt2, 1.41421, 0.));

} // BOOST_AUTO_TEST_CASE(test_RealComparisons)


// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(Vector2DComparisons_testcase) {

  constexpr lar::util::Vector2DComparison check{ 1e-5f };

  struct VectorBase {
    float x { 0.0 }, y { 0.0 };
    constexpr float X() const { return x; }
    constexpr float Y() const { return y; }
  };
  struct Vector1_t: VectorBase {};
  struct Vector2_t: VectorBase {};
  
  constexpr float const epsilon { std::numeric_limits<float>::epsilon() };
  
  constexpr Vector1_t Orig   {  0.0,  0.0 + 4.0*epsilon };
  constexpr Vector1_t Xpos   { +1.0,  0.0 - 8.0*epsilon };
  constexpr Vector2_t Xneg   { -1.0 + 2.0*epsilon,  0.0 };
  constexpr Vector2_t Ypos   {  0.0, +1.0 + 8.0*epsilon };
  constexpr Vector1_t Yneg   {  0.0, -1.0 + 2.0*epsilon };
  constexpr Vector1_t XYpos  { +3.0 + 1.0*epsilon, +4.0 };
  constexpr Vector2_t Xpos2  { +2.0 + 1.0*epsilon,  0.0 };
  constexpr Vector1_t Xneg2  { -2.0,  0.0 + 1.0*epsilon };
  constexpr Vector1_t Ypos2  {  0.0, +2.0 };
  constexpr Vector2_t Yneg2  {  0.0 - 4.0*epsilon, +1.0 };
  
  // zero()
  BOOST_TEST( check.zero(Orig ));
  BOOST_TEST(!check.zero(Xpos ));
  BOOST_TEST(!check.zero(Xneg ));
  BOOST_TEST(!check.zero(Ypos ));
  BOOST_TEST(!check.zero(Yneg ));
  BOOST_TEST(!check.zero(XYpos));
  BOOST_TEST(!check.zero(Xpos2));
  BOOST_TEST(!check.zero(Xneg2));
  BOOST_TEST(!check.zero(Ypos2));
  BOOST_TEST(!check.zero(Yneg2));
  
  // nonZero()
  BOOST_TEST(!check.nonZero(Orig ));
  BOOST_TEST( check.nonZero(Xpos ));
  BOOST_TEST( check.nonZero(Xneg ));
  BOOST_TEST( check.nonZero(Ypos ));
  BOOST_TEST( check.nonZero(Yneg ));
  BOOST_TEST( check.nonZero(XYpos));
  BOOST_TEST( check.nonZero(Xpos2));
  BOOST_TEST( check.nonZero(Xneg2));
  BOOST_TEST( check.nonZero(Ypos2));
  BOOST_TEST( check.nonZero(Yneg2));
  
  // equal()
  BOOST_TEST( check.equal(Orig , Orig ));
  BOOST_TEST(!check.equal(Orig , Xpos ));
  BOOST_TEST(!check.equal(Orig , Xneg ));
  BOOST_TEST(!check.equal(Orig , XYpos));
  BOOST_TEST(!check.equal(Xpos , Orig ));
  BOOST_TEST( check.equal(Xpos , Xpos ));
  BOOST_TEST(!check.equal(Xpos , Xneg ));
  BOOST_TEST(!check.equal(Xpos , XYpos));
  BOOST_TEST(!check.equal(Xneg , Orig ));
  BOOST_TEST(!check.equal(Xneg , Xpos ));
  BOOST_TEST( check.equal(Xneg , Xneg ));
  BOOST_TEST(!check.equal(Xneg , XYpos));
  BOOST_TEST(!check.equal(XYpos, Orig ));
  BOOST_TEST(!check.equal(XYpos, Xpos ));
  BOOST_TEST(!check.equal(XYpos, Xneg ));
  BOOST_TEST( check.equal(XYpos, XYpos));
  
  // nonEqual()
  BOOST_TEST(!check.nonEqual(Orig , Orig ));
  BOOST_TEST( check.nonEqual(Orig , Xpos ));
  BOOST_TEST( check.nonEqual(Orig , Xneg ));
  BOOST_TEST( check.nonEqual(Orig , XYpos));
  BOOST_TEST( check.nonEqual(Xpos , Orig ));
  BOOST_TEST(!check.nonEqual(Xpos , Xpos ));
  BOOST_TEST( check.nonEqual(Xpos , Xneg ));
  BOOST_TEST( check.nonEqual(Xpos , XYpos));
  BOOST_TEST( check.nonEqual(Xneg , Orig ));
  BOOST_TEST( check.nonEqual(Xneg , Xpos ));
  BOOST_TEST(!check.nonEqual(Xneg , Xneg ));
  BOOST_TEST( check.nonEqual(Xneg , XYpos));
  BOOST_TEST( check.nonEqual(XYpos, Orig ));
  BOOST_TEST( check.nonEqual(XYpos, Xpos ));
  BOOST_TEST( check.nonEqual(XYpos, Xneg ));
  BOOST_TEST(!check.nonEqual(XYpos, XYpos));
  
  // parallel()
  BOOST_TEST(!check.parallel(Orig , Orig ));
  BOOST_TEST(!check.parallel(Orig , Xpos ));
  BOOST_TEST(!check.parallel(Xpos , Orig ));
  BOOST_TEST(!check.parallel(Orig , Xneg ));
  BOOST_TEST(!check.parallel(Xneg , Orig ));
  BOOST_TEST(!check.parallel(Orig , XYpos));
  BOOST_TEST( check.parallel(Xpos , Xpos ));
  BOOST_TEST( check.parallel(Xpos , Xneg ));
  BOOST_TEST(!check.parallel(Xpos , Ypos ));
  BOOST_TEST(!check.parallel(Xpos , Ypos2));
  BOOST_TEST( check.parallel(Ypos , Ypos ));
  BOOST_TEST( check.parallel(Ypos , Ypos2));
  BOOST_TEST( check.parallel(Yneg , Ypos2));
  BOOST_TEST(!check.parallel(Ypos , XYpos));
  BOOST_TEST(!check.parallel(Yneg , XYpos));
  
  // parallel()
  BOOST_TEST(!check.parallel(Orig ,  0.0, Orig ,  0.0));
  BOOST_TEST(!check.parallel(Orig ,  1.0, Xpos ,  1.0));
  BOOST_TEST(!check.parallel(Xpos ,  1.0, Orig ,  1.0));
  BOOST_TEST(!check.parallel(Orig ,  1.0, Xneg ,  1.0));
  BOOST_TEST(!check.parallel(Xneg ,  1.0, Orig ,  1.0));
  BOOST_TEST(!check.parallel(Orig ,  1.0, XYpos,  1.0));
  BOOST_TEST( check.parallel(Xpos ,  1.0, Xpos ,  1.0));
  BOOST_TEST( check.parallel(Xpos ,  1.0, Xneg ,  1.0));
  BOOST_TEST(!check.parallel(Xpos ,  1.0, Ypos ,  1.0));
  BOOST_TEST(!check.parallel(Xpos ,  1.0, Ypos2,  4.0));
  BOOST_TEST( check.parallel(Ypos ,  1.0, Ypos ,  1.0));
  BOOST_TEST( check.parallel(Ypos ,  1.0, Ypos2,  4.0));
  BOOST_TEST( check.parallel(Yneg ,  1.0, Ypos2,  4.0));
  BOOST_TEST(!check.parallel(Ypos ,  1.0, XYpos,  1.0));
  BOOST_TEST(!check.parallel(Yneg ,  1.0, XYpos,  1.0));
  
  // nonParallel()
  BOOST_TEST( check.nonParallel(Orig , Orig ));
  BOOST_TEST( check.nonParallel(Orig , Xpos ));
  BOOST_TEST( check.nonParallel(Xpos , Orig ));
  BOOST_TEST( check.nonParallel(Orig , Xneg ));
  BOOST_TEST( check.nonParallel(Xneg , Orig ));
  BOOST_TEST( check.nonParallel(Orig , XYpos));
  BOOST_TEST(!check.nonParallel(Xpos , Xpos ));
  BOOST_TEST(!check.nonParallel(Xpos , Xneg ));
  BOOST_TEST( check.nonParallel(Xpos , Ypos ));
  BOOST_TEST( check.nonParallel(Xpos , Ypos2));
  BOOST_TEST(!check.nonParallel(Ypos , Ypos ));
  BOOST_TEST(!check.nonParallel(Ypos , Ypos2));
  BOOST_TEST(!check.nonParallel(Yneg , Ypos2));
  BOOST_TEST( check.nonParallel(Ypos , XYpos));
  BOOST_TEST( check.nonParallel(Yneg , XYpos));
  
  // nonParallel()
  BOOST_TEST( check.nonParallel(Orig ,  0.0, Orig ,  0.0));
  BOOST_TEST( check.nonParallel(Orig ,  1.0, Xpos ,  1.0));
  BOOST_TEST( check.nonParallel(Xpos ,  1.0, Orig ,  1.0));
  BOOST_TEST( check.nonParallel(Orig ,  1.0, Xneg ,  1.0));
  BOOST_TEST( check.nonParallel(Xneg ,  1.0, Orig ,  1.0));
  BOOST_TEST( check.nonParallel(Orig ,  1.0, XYpos,  1.0));
  BOOST_TEST(!check.nonParallel(Xpos ,  1.0, Xpos ,  1.0));
  BOOST_TEST(!check.nonParallel(Xpos ,  1.0, Xneg ,  1.0));
  BOOST_TEST( check.nonParallel(Xpos ,  1.0, Ypos ,  1.0));
  BOOST_TEST( check.nonParallel(Xpos ,  1.0, Ypos2,  4.0));
  BOOST_TEST(!check.nonParallel(Ypos ,  1.0, Ypos ,  1.0));
  BOOST_TEST(!check.nonParallel(Ypos ,  1.0, Ypos2,  4.0));
  BOOST_TEST(!check.nonParallel(Yneg ,  1.0, Ypos2,  4.0));
  BOOST_TEST( check.nonParallel(Ypos ,  1.0, XYpos,  1.0));
  BOOST_TEST( check.nonParallel(Yneg ,  1.0, XYpos,  1.0));
  
  // orthogonal()
  BOOST_TEST(!check.orthogonal(Xpos , Xpos ));
  BOOST_TEST( check.orthogonal(Xpos , Ypos ));
  BOOST_TEST(!check.orthogonal(Xpos , Xneg ));
  BOOST_TEST( check.orthogonal(Xpos , Yneg ));
  BOOST_TEST(!check.orthogonal(Xpos , Xpos2));
  BOOST_TEST( check.orthogonal(Xpos , Ypos2));
  BOOST_TEST( check.orthogonal(Ypos , Xpos ));
  BOOST_TEST(!check.orthogonal(Ypos , Ypos ));
  BOOST_TEST( check.orthogonal(Ypos , Xneg ));
  BOOST_TEST(!check.orthogonal(Ypos , Yneg ));
  BOOST_TEST( check.orthogonal(Ypos , Xpos2));
  BOOST_TEST(!check.orthogonal(Ypos , Ypos2));
  BOOST_TEST(!check.orthogonal(XYpos, XYpos));
  BOOST_TEST(!check.orthogonal(XYpos, Xpos ));
  BOOST_TEST(!check.orthogonal(XYpos, Ypos ));
  BOOST_TEST(!check.orthogonal(XYpos, Ypos2));
  
  // nonOrthogonal()
  BOOST_TEST( check.nonOrthogonal(Xpos , Xpos ));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Ypos ));
  BOOST_TEST( check.nonOrthogonal(Xpos , Xneg ));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Yneg ));
  BOOST_TEST( check.nonOrthogonal(Xpos , Xpos2));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Ypos2));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Xpos ));
  BOOST_TEST( check.nonOrthogonal(Ypos , Ypos ));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Xneg ));
  BOOST_TEST( check.nonOrthogonal(Ypos , Yneg ));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Xpos2));
  BOOST_TEST( check.nonOrthogonal(Ypos , Ypos2));
  BOOST_TEST( check.nonOrthogonal(XYpos, XYpos));
  BOOST_TEST( check.nonOrthogonal(XYpos, Xpos ));
  BOOST_TEST( check.nonOrthogonal(XYpos, Ypos ));
  BOOST_TEST( check.nonOrthogonal(XYpos, Ypos2));
  
} // BOOST_AUTO_TEST_CASE(Vector2DComparisons_testcase)


// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(Vector3DComparisons_testcase) {

  constexpr lar::util::Vector3DComparison check{ 1e-5f };

  struct VectorBase {
    float x { 0.0 }, y { 0.0 }, z { 0.0 };
    constexpr float X() const { return x; }
    constexpr float Y() const { return y; }
    constexpr float Z() const { return z; }
  };
  struct Vector1_t: VectorBase {};
  struct Vector2_t: VectorBase {};
  
  constexpr float const epsilon { std::numeric_limits<float>::epsilon() };
  
  constexpr Vector1_t Orig   {  0.0,  0.0,  0.0 + 4.0*epsilon };
  constexpr Vector1_t Xpos   { +1.0,  0.0 + 4.0*epsilon,  0.0 };
  constexpr Vector2_t Xneg   { -1.0 + 2.0*epsilon,  0.0,  0.0 };
  constexpr Vector2_t Ypos   {  0.0, +1.0,  0.0 + 8.0*epsilon };
  constexpr Vector1_t Yneg   {  0.0, -1.0 + 2.0*epsilon,  0.0 };
  constexpr Vector2_t Zpos   {  0.0 + 1.0*epsilon,  0.0, +1.0 };
  constexpr Vector1_t Zneg   {  0.0,  0.0, -1.0 + 8.0*epsilon };
  constexpr Vector1_t XYpos  { +3.0 + 1.0*epsilon, +4.0,  0.0 + 1.0*epsilon };
  constexpr Vector2_t Xpos2  { +2.0 + 1.0*epsilon,  0.0,  0.0 - 2.0*epsilon };
  constexpr Vector1_t Xneg2  { -2.0,  0.0 + 1.0*epsilon,  0.0 + 8.0*epsilon };
  constexpr Vector1_t Ypos2  {  0.0, +2.0,  0.0 };
  constexpr Vector2_t Yneg2  {  0.0 - 4.0*epsilon, -2.0,  0.0 };
  constexpr Vector1_t Zpos2  {  0.0,  0.0 + 2.0*epsilon, +2.0 };
  constexpr Vector2_t Zneg2  {  0.0,  0.0, -2.0 + 2.0*epsilon };
  
  // zero()
  BOOST_TEST( check.zero(Orig ));
  BOOST_TEST(!check.zero(Xpos ));
  BOOST_TEST(!check.zero(Xneg ));
  BOOST_TEST(!check.zero(Ypos ));
  BOOST_TEST(!check.zero(Yneg ));
  BOOST_TEST(!check.zero(Zpos ));
  BOOST_TEST(!check.zero(Zneg ));
  BOOST_TEST(!check.zero(XYpos));
  BOOST_TEST(!check.zero(Xpos2));
  BOOST_TEST(!check.zero(Xneg2));
  BOOST_TEST(!check.zero(Ypos2));
  BOOST_TEST(!check.zero(Yneg2));
  BOOST_TEST(!check.zero(Zpos2));
  BOOST_TEST(!check.zero(Zneg2));
  
  // nonZero()
  BOOST_TEST(!check.nonZero(Orig ));
  BOOST_TEST( check.nonZero(Xpos ));
  BOOST_TEST( check.nonZero(Xneg ));
  BOOST_TEST( check.nonZero(Ypos ));
  BOOST_TEST( check.nonZero(Yneg ));
  BOOST_TEST( check.nonZero(Zpos ));
  BOOST_TEST( check.nonZero(Zneg ));
  BOOST_TEST( check.nonZero(XYpos));
  BOOST_TEST( check.nonZero(Xpos2));
  BOOST_TEST( check.nonZero(Xneg2));
  BOOST_TEST( check.nonZero(Ypos2));
  BOOST_TEST( check.nonZero(Yneg2));
  BOOST_TEST( check.nonZero(Zpos2));
  BOOST_TEST( check.nonZero(Zneg2));
  
  // equal()
  BOOST_TEST( check.equal(Orig , Orig ));
  BOOST_TEST(!check.equal(Orig , Xpos ));
  BOOST_TEST(!check.equal(Orig , Xneg ));
  BOOST_TEST(!check.equal(Orig , XYpos));
  BOOST_TEST(!check.equal(Xpos , Orig ));
  BOOST_TEST( check.equal(Xpos , Xpos ));
  BOOST_TEST(!check.equal(Xpos , Xneg ));
  BOOST_TEST(!check.equal(Xpos , XYpos));
  BOOST_TEST(!check.equal(Xneg , Orig ));
  BOOST_TEST(!check.equal(Xneg , Xpos ));
  BOOST_TEST( check.equal(Xneg , Xneg ));
  BOOST_TEST(!check.equal(Xneg , XYpos));
  BOOST_TEST(!check.equal(XYpos, Orig ));
  BOOST_TEST(!check.equal(XYpos, Xpos ));
  BOOST_TEST(!check.equal(XYpos, Xneg ));
  BOOST_TEST( check.equal(XYpos, XYpos));
  
  // nonEqual()
  BOOST_TEST(!check.nonEqual(Orig , Orig ));
  BOOST_TEST( check.nonEqual(Orig , Xpos ));
  BOOST_TEST( check.nonEqual(Orig , Xneg ));
  BOOST_TEST( check.nonEqual(Orig , XYpos));
  BOOST_TEST( check.nonEqual(Xpos , Orig ));
  BOOST_TEST(!check.nonEqual(Xpos , Xpos ));
  BOOST_TEST( check.nonEqual(Xpos , Xneg ));
  BOOST_TEST( check.nonEqual(Xpos , XYpos));
  BOOST_TEST( check.nonEqual(Xneg , Orig ));
  BOOST_TEST( check.nonEqual(Xneg , Xpos ));
  BOOST_TEST(!check.nonEqual(Xneg , Xneg ));
  BOOST_TEST( check.nonEqual(Xneg , XYpos));
  BOOST_TEST( check.nonEqual(XYpos, Orig ));
  BOOST_TEST( check.nonEqual(XYpos, Xpos ));
  BOOST_TEST( check.nonEqual(XYpos, Xneg ));
  BOOST_TEST(!check.nonEqual(XYpos, XYpos));
  
  // parallel()
  BOOST_TEST(!check.parallel(Orig , Orig ));
  BOOST_TEST(!check.parallel(Orig , Xpos ));
  BOOST_TEST(!check.parallel(Xpos , Orig ));
  BOOST_TEST(!check.parallel(Orig , Xneg ));
  BOOST_TEST(!check.parallel(Xneg , Orig ));
  BOOST_TEST(!check.parallel(Orig , XYpos));
  BOOST_TEST( check.parallel(Xpos , Xpos ));
  BOOST_TEST( check.parallel(Xpos , Xneg ));
  BOOST_TEST(!check.parallel(Xpos , Ypos ));
  BOOST_TEST(!check.parallel(Xpos , Ypos2));
  BOOST_TEST( check.parallel(Ypos , Ypos ));
  BOOST_TEST( check.parallel(Ypos , Ypos2));
  BOOST_TEST( check.parallel(Yneg , Ypos2));
  BOOST_TEST(!check.parallel(Ypos , XYpos));
  BOOST_TEST(!check.parallel(Yneg , XYpos));
  BOOST_TEST( check.parallel(Zpos , Zpos ));
  BOOST_TEST( check.parallel(Zpos , Zpos2));
  BOOST_TEST(!check.parallel(Zneg , Ypos2));
  BOOST_TEST(!check.parallel(Zpos , XYpos));
  
  // parallel()
  BOOST_TEST(!check.parallel(Orig ,  0.0, Orig ,  0.0));
  BOOST_TEST(!check.parallel(Orig ,  1.0, Xpos ,  1.0));
  BOOST_TEST(!check.parallel(Xpos ,  1.0, Orig ,  1.0));
  BOOST_TEST(!check.parallel(Orig ,  1.0, Xneg ,  1.0));
  BOOST_TEST(!check.parallel(Xneg ,  1.0, Orig ,  1.0));
  BOOST_TEST(!check.parallel(Orig ,  1.0, XYpos,  1.0));
  BOOST_TEST( check.parallel(Xpos ,  1.0, Xpos ,  1.0));
  BOOST_TEST( check.parallel(Xpos ,  1.0, Xneg ,  1.0));
  BOOST_TEST(!check.parallel(Xpos ,  1.0, Ypos ,  1.0));
  BOOST_TEST(!check.parallel(Xpos ,  1.0, Ypos2,  4.0));
  BOOST_TEST( check.parallel(Ypos ,  1.0, Ypos ,  1.0));
  BOOST_TEST( check.parallel(Ypos ,  1.0, Ypos2,  4.0));
  BOOST_TEST( check.parallel(Yneg ,  1.0, Ypos2,  4.0));
  BOOST_TEST(!check.parallel(Ypos ,  1.0, XYpos,  1.0));
  BOOST_TEST(!check.parallel(Yneg ,  1.0, XYpos,  1.0));
  BOOST_TEST( check.parallel(Zpos ,  1.0, Zpos ,  1.0));
  BOOST_TEST( check.parallel(Zpos ,  1.0, Zpos2,  4.0));
  BOOST_TEST(!check.parallel(Zneg ,  1.0, Ypos2,  4.0));
  BOOST_TEST(!check.parallel(Zpos ,  1.0, XYpos, 25.0));
  
  // nonParallel()
  BOOST_TEST( check.nonParallel(Orig , Orig ));
  BOOST_TEST( check.nonParallel(Orig , Xpos ));
  BOOST_TEST( check.nonParallel(Xpos , Orig ));
  BOOST_TEST( check.nonParallel(Orig , Xneg ));
  BOOST_TEST( check.nonParallel(Xneg , Orig ));
  BOOST_TEST( check.nonParallel(Orig , XYpos));
  BOOST_TEST(!check.nonParallel(Xpos , Xpos ));
  BOOST_TEST(!check.nonParallel(Xpos , Xneg ));
  BOOST_TEST( check.nonParallel(Xpos , Ypos ));
  BOOST_TEST( check.nonParallel(Xpos , Ypos2));
  BOOST_TEST(!check.nonParallel(Ypos , Ypos ));
  BOOST_TEST(!check.nonParallel(Ypos , Ypos2));
  BOOST_TEST(!check.nonParallel(Yneg , Ypos2));
  BOOST_TEST( check.nonParallel(Ypos , XYpos));
  BOOST_TEST( check.nonParallel(Yneg , XYpos));
  BOOST_TEST(!check.nonParallel(Zpos , Zpos ));
  BOOST_TEST(!check.nonParallel(Zpos , Zpos2));
  BOOST_TEST( check.nonParallel(Zneg , Ypos2));
  BOOST_TEST( check.nonParallel(Zpos , XYpos));
  
  // nonParallel()
  BOOST_TEST( check.nonParallel(Orig ,  0.0, Orig ,  0.0));
  BOOST_TEST( check.nonParallel(Orig ,  1.0, Xpos ,  1.0));
  BOOST_TEST( check.nonParallel(Xpos ,  1.0, Orig ,  1.0));
  BOOST_TEST( check.nonParallel(Orig ,  1.0, Xneg ,  1.0));
  BOOST_TEST( check.nonParallel(Xneg ,  1.0, Orig ,  1.0));
  BOOST_TEST( check.nonParallel(Orig ,  1.0, XYpos,  1.0));
  BOOST_TEST(!check.nonParallel(Xpos ,  1.0, Xpos ,  1.0));
  BOOST_TEST(!check.nonParallel(Xpos ,  1.0, Xneg ,  1.0));
  BOOST_TEST( check.nonParallel(Xpos ,  1.0, Ypos ,  1.0));
  BOOST_TEST( check.nonParallel(Xpos ,  1.0, Ypos2,  4.0));
  BOOST_TEST(!check.nonParallel(Ypos ,  1.0, Ypos ,  1.0));
  BOOST_TEST(!check.nonParallel(Ypos ,  1.0, Ypos2,  4.0));
  BOOST_TEST(!check.nonParallel(Yneg ,  1.0, Ypos2,  4.0));
  BOOST_TEST( check.nonParallel(Ypos ,  1.0, XYpos,  1.0));
  BOOST_TEST( check.nonParallel(Yneg ,  1.0, XYpos,  1.0));
  BOOST_TEST(!check.nonParallel(Zpos ,  1.0, Zpos ,  1.0));
  BOOST_TEST(!check.nonParallel(Zpos ,  1.0, Zpos2,  4.0));
  BOOST_TEST( check.nonParallel(Zneg ,  1.0, Ypos2,  4.0));
  BOOST_TEST( check.nonParallel(Zpos ,  1.0, XYpos, 25.0));
  
  // orthogonal()
  BOOST_TEST(!check.orthogonal(Xpos , Xpos ));
  BOOST_TEST( check.orthogonal(Xpos , Ypos ));
  BOOST_TEST( check.orthogonal(Xpos , Zpos ));
  BOOST_TEST(!check.orthogonal(Xpos , Xneg ));
  BOOST_TEST( check.orthogonal(Xpos , Yneg ));
  BOOST_TEST( check.orthogonal(Xpos , Zneg ));
  BOOST_TEST(!check.orthogonal(Xpos , Xpos2));
  BOOST_TEST( check.orthogonal(Xpos , Ypos2));
  BOOST_TEST( check.orthogonal(Xpos , Zpos2));
  BOOST_TEST( check.orthogonal(Ypos , Xpos ));
  BOOST_TEST(!check.orthogonal(Ypos , Ypos ));
  BOOST_TEST( check.orthogonal(Ypos , Zpos ));
  BOOST_TEST( check.orthogonal(Ypos , Xneg ));
  BOOST_TEST(!check.orthogonal(Ypos , Yneg ));
  BOOST_TEST( check.orthogonal(Ypos , Zneg ));
  BOOST_TEST( check.orthogonal(Ypos , Xpos2));
  BOOST_TEST(!check.orthogonal(Ypos , Ypos2));
  BOOST_TEST( check.orthogonal(Ypos , Zpos2));
  BOOST_TEST( check.orthogonal(Zpos , Xpos ));
  BOOST_TEST( check.orthogonal(Zpos , Ypos ));
  BOOST_TEST(!check.orthogonal(Zpos , Zpos ));
  BOOST_TEST( check.orthogonal(Zpos , Xpos2));
  BOOST_TEST( check.orthogonal(Zpos , Ypos2));
  BOOST_TEST(!check.orthogonal(Zpos , Zpos2));
  BOOST_TEST( check.orthogonal(Zpos , Xneg2));
  BOOST_TEST( check.orthogonal(Zpos , Yneg2));
  BOOST_TEST(!check.orthogonal(Zpos , Zneg2));
  BOOST_TEST(!check.orthogonal(XYpos, XYpos));
  BOOST_TEST(!check.orthogonal(XYpos, Xpos ));
  BOOST_TEST(!check.orthogonal(XYpos, Ypos ));
  BOOST_TEST(!check.orthogonal(XYpos, Ypos2));
  BOOST_TEST( check.orthogonal(XYpos, Zpos2));
  
  // nonOrthogonal()
  BOOST_TEST( check.nonOrthogonal(Xpos , Xpos ));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Ypos ));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Zpos ));
  BOOST_TEST( check.nonOrthogonal(Xpos , Xneg ));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Yneg ));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Zneg ));
  BOOST_TEST( check.nonOrthogonal(Xpos , Xpos2));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Ypos2));
  BOOST_TEST(!check.nonOrthogonal(Xpos , Zpos2));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Xpos ));
  BOOST_TEST( check.nonOrthogonal(Ypos , Ypos ));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Zpos ));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Xneg ));
  BOOST_TEST( check.nonOrthogonal(Ypos , Yneg ));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Zneg ));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Xpos2));
  BOOST_TEST( check.nonOrthogonal(Ypos , Ypos2));
  BOOST_TEST(!check.nonOrthogonal(Ypos , Zpos2));
  BOOST_TEST(!check.nonOrthogonal(Zpos , Xpos ));
  BOOST_TEST(!check.nonOrthogonal(Zpos , Ypos ));
  BOOST_TEST( check.nonOrthogonal(Zpos , Zpos ));
  BOOST_TEST(!check.nonOrthogonal(Zpos , Xpos2));
  BOOST_TEST(!check.nonOrthogonal(Zpos , Ypos2));
  BOOST_TEST( check.nonOrthogonal(Zpos , Zpos2));
  BOOST_TEST(!check.nonOrthogonal(Zpos , Xneg2));
  BOOST_TEST(!check.nonOrthogonal(Zpos , Yneg2));
  BOOST_TEST( check.nonOrthogonal(Zpos , Zneg2));
  BOOST_TEST( check.nonOrthogonal(XYpos, XYpos));
  BOOST_TEST( check.nonOrthogonal(XYpos, Xpos ));
  BOOST_TEST( check.nonOrthogonal(XYpos, Ypos ));
  BOOST_TEST( check.nonOrthogonal(XYpos, Ypos2));
  BOOST_TEST(!check.nonOrthogonal(XYpos, Zpos2));
  
} // BOOST_AUTO_TEST_CASE(Vector3DComparisons_testcase)


// -----------------------------------------------------------------------------
