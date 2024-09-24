/**
 * @file   RealComparisons.h
 * @brief  Class for approximate comparisons
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 2, 2016
 *
 * This is a header-only library.
 *
 */

#ifndef LARCORE_COREUTILS_REALCOMPARISONS_H
#define LARCORE_COREUTILS_REALCOMPARISONS_H

// C/C++ standard libraries
#include <cmath>   // std::abs()
#include <utility> // std::move()

namespace lar::util {

  /** ************************************************************************
   * @brief Provides simple real number checks
   * @tparam RealType type of value to operate on
   *
   * This class provides some members to perform comparisons between real
   * numbers, allowing for some tolerance for rounding errors.
   *
   * The tolerance parameter is fixed.
   * Example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * lar::util::RealComparisons<float> check(1e-5);
   *
   * double const sqrt2 = std::sqrt(2);
   *
   * // it should print the message
   * if (check.equal(sqrt2, 1.4142)) {
   *   std::cout
   *     << "Square root of 2 is not even close to 1.4142." << std::endl;
   * }
   *
   * // it should not print the message
   * if (check.equal(sqrt2, 1.414213)) {
   *   std::cout
   *     << "Square root of 2 is not even close to 1.414213." << std::endl;
   * }
   *
   * // this should print the message
   * if (check.within(std::sqrt(2), 0., 1.41421)) {
   *   std::cout
   *     << "Square root of 2 is between 0 and 1.41421 (tops)." << std::endl;
   * }
   *
   * // this will probably print the message
   * double const epsilon = 2.0 - (sqrt2 * sqrt2);
   * if (check.zero(epsilon)) {
   *   std::cout
   *     << "The square of the square root of 2 is roughly 2." << std::endl;
   * }
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  template <typename RealType>
  struct RealComparisons {
    using Value_t = RealType; /// type of values being compered

    /// Constructor: specify the threshold
    constexpr RealComparisons(Value_t threshold) : threshold(threshold) {}

    /// Returns whether the value is no farther from 0 than the threshold
    constexpr bool zero(Value_t value) const { return std::abs(value) <= threshold; }

    /// Returns whether the value is farther from 0 than the threshold
    constexpr bool nonZero(Value_t value) const { return !zero(value); }

    /// Returns whether a and b are no farther than the threshold
    constexpr bool equal(Value_t a, Value_t b) const { return zero(a - b); }

    /// Returns whether a and b are farther than the threshold
    constexpr bool nonEqual(Value_t a, Value_t b) const { return !equal(a, b); }

    /// Returns whether value is larger than zero beyond tolerance
    constexpr bool strictlyNegative(Value_t value) const { return value < -threshold; }

    /// Returns whether value is smaller than zero beyond tolerance
    constexpr bool strictlyPositive(Value_t value) const { return value > threshold; }

    /// Returns whether value is larger than or `equal()` to zero
    constexpr bool nonNegative(Value_t value) const { return value >= -threshold; }

    /// Returns whether value is smaller than or `equal()` to zero
    constexpr bool nonPositive(Value_t value) const { return value <= threshold; }

    /// Returns whether a is strictly smaller than b
    constexpr bool strictlySmaller(Value_t a, Value_t b) const { return strictlyNegative(a - b); }

    /// Returns whether a is greater than (or equal to) b
    constexpr bool nonSmaller(Value_t a, Value_t b) const { return nonNegative(a - b); }

    /// Returns whether a is strictly greater than b
    constexpr bool strictlyGreater(Value_t a, Value_t b) const { return strictlyPositive(a - b); }

    /// Returns whether a is smaller than (or equal to) b
    constexpr bool nonGreater(Value_t a, Value_t b) const { return nonPositive(a - b); }

    /// Returns whether value is between the bounds (included)
    constexpr bool within(Value_t value, Value_t lower, Value_t upper) const
    {
      return nonNegative(value - lower) && nonPositive(value - upper);
    }

    /// Returns whether value is between bounds (included); bounds are sorted
    constexpr bool withinSorted(Value_t value, Value_t lower, Value_t upper) const
    {
      return (lower < upper) ? within(value, lower, upper) : within(value, upper, lower);
    } // sortedWithin

    Value_t threshold; /// Threshold to compare the values to

  }; // struct RealComparisons<>

  //--------------------------------------------------------------------------
  /// Class comparing 2D vectors
  template <typename RealType>
  struct Vector2DComparison {

    using Comp_t = RealComparisons<RealType>;

    /// Copy the specified comparison.
    constexpr Vector2DComparison(Comp_t const& comparer) : comparer(comparer) {}

    /// Steal the specified comparison.
    Vector2DComparison(Comp_t&& comparer) : comparer(std::move(comparer)) {}

    /// Use the specified threshold.
    constexpr Vector2DComparison(RealType threshold) : comparer(threshold) {}

    /// Returns the basic value comparer.
    constexpr Comp_t comp() const { return comparer; }

    /// Returns whether the specified vector is null (within tolerance).
    template <typename Vect>
    constexpr bool zero(Vect const& v) const
    {
      return comp().zero(v.X()) && comp().zero(v.Y());
    }

    /// Returns whether the specified vector is not null (within tolerance).
    template <typename Vect>
    constexpr bool nonZero(Vect const& v) const
    {
      return !zero(v);
    }

    /// Returns whether the specified vectors match (within tolerance).
    template <typename VectA, typename VectB>
    constexpr bool equal(VectA const& a, VectB const& b) const
    {
      return comp().equal(a.X(), b.X()) && comp().equal(a.Y(), b.Y());
    }

    /// Returns whether the specified vectors do not match (within tolerance).
    template <typename VectA, typename VectB>
    constexpr bool nonEqual(VectA const& a, VectB const& b) const
    {
      return !equal(a, b);
    }

  private:
    Comp_t const comparer; ///< Comparison object.

  }; // struct Vector2DComparison

  //--------------------------------------------------------------------------
  /// Creates a `Vector2DComparison` from a `RealComparisons` object.
  template <typename RealType>
  auto makeVector2DComparison(RealType threshold)
  {
    return Vector2DComparison<RealType>(threshold);
  }

  /// Creates a `Vector2DComparison` from a `RealComparisons` object.
  template <typename RealType>
  auto makeVector2DComparison(lar::util::RealComparisons<RealType> const& comp)
  {
    return Vector2DComparison<RealType>(comp);
  }

  //--------------------------------------------------------------------------
  /// Class comparing 2D vectors.
  template <typename RealType>
  struct Vector3DComparison {

    /// Type of base value comparer.
    using Comp_t = RealComparisons<RealType>;

    /// Type of 2D vector comparer.
    using Comp2D_t = Vector2DComparison<typename Comp_t::Value_t>;

    /// Copy the specified comparison.
    constexpr Vector3DComparison(Comp_t const& comparer) : comparer(comparer) {}

    /// Steal the specified comparison.
    Vector3DComparison(Comp_t&& comparer) : comparer(std::move(comparer)) {}

    /// Use the specified threshold.
    constexpr Vector3DComparison(RealType threshold) : comparer(threshold) {}

    /// Returns the base value comparer.
    constexpr Comp_t comp() const { return comp2D().comp(); }

    /// Returns the 2D vector comparer.
    constexpr Comp2D_t comp2D() const { return comparer; }

    /// Returns whether the specified vector is null (within tolerance).
    template <typename Vect>
    constexpr bool zero(Vect const& v) const
    {
      return comp2D().zero(v) && comp().zero(v.Z());
    }

    /// Returns whether the specified vector is not null (within tolerance).
    template <typename Vect>
    constexpr bool nonZero(Vect const& v) const
    {
      return !zero(v);
    }

    /// Returns whether the specified vectors match (within tolerance).
    template <typename VectA, typename VectB>
    constexpr bool equal(VectA const& a, VectB const& b) const
    {
      return comp2D().equal(a, b) && comp().equal(a.Z(), b.Z());
    }

    /// Returns whether the specified vectors do not match (within tolerance).
    template <typename VectA, typename VectB>
    constexpr bool nonEqual(VectA const& a, VectB const& b) const
    {
      return !equal(a, b);
    }

  private:
    Comp2D_t comparer;

  }; // struct Vector3DComparison

  //--------------------------------------------------------------------------
  /// Creates a `Vector3DComparison` from a `RealComparisons` object.
  template <typename RealType>
  auto makeVector3DComparison(RealType threshold)
  {
    return Vector3DComparison<RealType>(threshold);
  }

  /// Creates a `Vector3DComparison` from a `RealComparisons` object.
  template <typename RealType>
  auto makeVector3DComparison(lar::util::RealComparisons<RealType> const& comp)
  {
    return Vector3DComparison<RealType>(comp);
  }

  //--------------------------------------------------------------------------

} // namespace lar::util

#endif // LARCORE_COREUTILS_REALCOMPARISONS_H
