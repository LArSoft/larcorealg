/**
 * @file   larcorealg/CoreUtils/counter.h
 * @brief  Test of `util::counter` and support utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 14, 2019
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_COREUTILS_COUNTER_H
#define LARCOREALG_COREUTILS_COUNTER_H

// External libraries
#include "range/v3/view.hpp"

// C/C++ libraries
#include <cstddef> // std::size_t
#include <utility> // std::move()

namespace util {

  // -- BEGIN -- Counted iterations --------------------------------------------
  /// @name Counted iterations
  /// @{

  /**
   * @brief An iterator dereferencing to a counter value.
   * @tparam T the type of counter returned on dereferenciation
   * @see `util::counter()`
   *
   * This iterator returns on dereferencing the net count of how many times it
   * has been incremented.
   * So for example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> data;
   * for (util::count_iterator it;; ++it) {
   *   if (*it >= 10) break;
   *   data.push_back(*it); // implicit conversion `std::size_t` to `int`
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * this infinite loop will push the value `0` into `data` on the first
   * iteration, `1` on the second and so forth, until at the eleventh iteration,
   * just before it can push `10`, the loop is forcibly broken leaving 10
   * elements in `data`.
   *
   * The iterator can be made to start from an index `n` different than `0`, in
   * which case it behaves like it had already been incremented `n` times.
   * End iterators can be build in that way too; see also `util::counter()`.
   *
   */
  template <typename T = std::size_t>
  class count_iterator {
  public:
    // --- BEGIN -- Traits and data types --------------------------------------
    /// @name Traits and data types
    /// @{

    using iterator_type = count_iterator<T>; ///< Type of this iterator.

    using difference_type = std::ptrdiff_t;
    using value_type = T;       ///< Type of index returned by this iterator.
    using reference = T const&; ///< Type returned by dereference operator.
    using pointer = T*;
    using iterator_category = std::bidirectional_iterator_tag;

    /// @}
    // --- END -- Traits and data types ----------------------------------------

    // --- BEGIN Constructors --------------------------------------------------
    /// @name Construction
    /// @{

    /**
     * @brief Initializes the iterator.
     *
     * The initial loop count is the default-constructed value of the counter
     * type, which is usually some variation on the concept of `0`.
     */
    count_iterator() = default;

    /**
     * @brief Initializes the iterator with the specified loop count.
     * @param count the initial loop count
     */
    count_iterator(value_type count) : fCount(count) {}

    /// @}
    // --- END Constructors ----------------------------------------------------

    // --- BEGIN -- Data access ------------------------------------------------
    /// @name Data access
    /// @{

    /// Returns the current loop count.
    reference operator*() const { return fCount; }

    /// @}
    // --- END -- Data access --------------------------------------------------

    // --- BEGIN -- Modification -----------------------------------------------
    /// @name Modification
    /// @{

    /// Increments the loop count of this iterator, which is then returned.
    iterator_type& operator++()
    {
      ++fCount;
      return *this;
    }

    /// Increments the loop count of this iterator, returning a copy with the
    /// value before the increment.
    iterator_type operator++(int) const
    {
      iterator_type const old = *this;
      operator++();
      return old;
    }

    /// Decrements the loop count of this iterator, which is then returned.
    iterator_type& operator--()
    {
      --fCount;
      return *this;
    }

    /// Decrements the loop count of this iterator, returning a copy with the
    /// value before the decrement.
    iterator_type operator--(int) const
    {
      iterator_type const old = *this;
      operator--();
      return old;
    }

    /// @}
    // --- END -- Modification -------------------------------------------------

    // --- BEGIN -- Comparisons ------------------------------------------------
    /// @name Comparisons
    /// @{

    /// Iterators are equal if their loop counts compare equal.
    template <typename U>
    bool operator==(count_iterator<U> const& other) const
    {
      return fCount == other.fCount;
    }

    /// Iterators are equal if their loop counts compare different.
    template <typename U>
    bool operator!=(count_iterator<U> const& other) const
    {
      return fCount != other.fCount;
    }

    /// @}
    // --- END -- Comparisons --------------------------------------------------

  private:
    value_type fCount{}; ///< Internal counter.

  }; // class count_iterator<>

  /**
   * @brief Returns an object to iterate values from `begin` to `end` in a
   *        range-for loop.
   * @tparam T type of counter value
   * @return a control object for range-for loop
   * @see `util::count_iterator`
   *
   * An example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<std::size_t> data;
   * for (auto i: util::counter(4, 8)) {
   *   data.push_back(i);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will insert in `data` the numbers from `4` to `7`, just like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * for (std::size_t i = 4; i < 8; ++i) {
   *   data.push_back(i);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * would.
   */
  template <typename T>
  auto counter(T begin, T end);

  /// Version of `util::counter()` starting at default-constructed `T`
  /// (usually some form of `0`).
  template <typename T>
  auto counter(T end)
  {
    return counter(T{}, end);
  }

  /**
   * @brief Version of `util::counter()` starting at `begin` and never ending.
   * @tparam T type of counter value
   * @param begin the count to start from (default-constructed, usually some
   *              form of `0`, by default)
   * @return a control object for range-for loop
   *
   * An example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<unsigned char> data;
   * for (auto ch: util::infinite_counter<unsigned char>()) {
   *   if (data.size() >= 512U) break;
   *   data.push_back(ch);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * This loop runs through the full range of a character (8 bits, supposedly)
   * twice before being broken.
   *
   */
  template <typename T = std::size_t>
  auto infinite_counter(T begin = T{});

  /// @}
  // -- END -- Counted iterations ----------------------------------------------

} // namespace util

//==============================================================================
//=== template implementation
//==============================================================================
//------------------------------------------------------------------------------
//--- util::counter()
//------------------------------------------------------------------------------
namespace util::details {

  /**
   * @brief Class used as end iterator (sentinel) for an infinite loop.
   * @tparam T the nominal count type
   *
   * This "iterator" offers a very limited interface, not at all like a real
   * iterator.
   * In fact, it can only be compared to other `count_iterator` objects on the
   * same counter type.
   */
  template <typename T>
  class infinite_endcount_iterator {
    using this_iterator_t = infinite_endcount_iterator<T>;
    using count_iterator_t = count_iterator<T>;

  public:
    // mock-up of stuff required by `ranges::subrange`
    using value_type = T;
    using reference = T const&;
    using pointer = T const*;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag; // this is a lie

    /// Never admit this iterator is equal to anything else (except the same).
    bool operator==(this_iterator_t const&) const { return true; }

    /// Never admit this iterator is equal to anything else (except the same).
    bool operator!=(this_iterator_t const&) const { return false; }

  }; // class infinite_endcount_iterator

  /// Never admit a `infinite_endcount_iterator` to be equal to anything else.
  template <typename T>
  bool operator!=(infinite_endcount_iterator<T> const&, count_iterator<T> const&)
  {
    return true;
  }

  template <typename T>
  bool operator!=(count_iterator<T> const&, infinite_endcount_iterator<T> const&)
  {
    return true;
  }

  template <typename T>
  bool operator==(infinite_endcount_iterator<T> const&, count_iterator<T> const&)
  {
    return false;
  }

  template <typename T>
  bool operator==(count_iterator<T> const&, infinite_endcount_iterator<T> const&)
  {
    return false;
  }

} // namespace util::details

//------------------------------------------------------------------------------
template <typename T>
auto util::counter(T begin, T end)
{
  return ranges::make_subrange(count_iterator(begin), count_iterator(end));
}

//------------------------------------------------------------------------------
template <typename T>
auto util::infinite_counter(T begin)
{
  return ranges::make_subrange(count_iterator(begin), details::infinite_endcount_iterator<T>());
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_COUNTER_H
