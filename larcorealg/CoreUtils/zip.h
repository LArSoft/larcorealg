/**
 * @file   larcorealg/CoreUtils/zip.h
 * @brief  Definition of `util::zip()`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 14, 2019
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_COREUTILS_ZIP_H
#define LARCOREALG_COREUTILS_ZIP_H

#include "larcorealg/CoreUtils/span.h"

// C/C++ libraries
#include <cstddef>  // std::size_t
#include <iterator> // std::begin(), std::end()
#include <tuple>
#include <type_traits> // std::remove_cv_t<>, ...
#include <utility>     // std::forward(), std::index_sequence_for(), ...

namespace util {

  // -- BEGIN -- Parallel iterations -------------------------------------------
  /// @name Parallel iterations
  /// @{

  /**
   * @brief Range-for loop helper iterating across many collections at the
   *        same time.
   * @tparam Lead index of the parameter driving the start and end of the loop
   * @tparam Iterables type of objects to be iterated together
   * @param iterables all iterable objects to be iterated together
   * @return an object suitable for range-for loop
   * @see `util::enumerate()`
   *
   * In the range-for loop, at each iteration this object yields a `tuple` of
   * values, each of the type returned by dereferencing `begin(iterable)`.
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * constexpr std::size_t N = 4;
   * std::array<int, N> twice;
   * std::vector<double> thrice(N + 1);
   *
   * unsigned int i = 0;
   * for (auto&& [ a, b]: util::zip(twice, thrice)) {
   *
   *   a = 2 * i;
   *   b = 3.0 * i;
   *
   *   ++i;
   *
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this example, `N` iterations will be run because that is the size of
   * the first iterable given to `enumerate`. If a different leading iterable
   * is needed, that has to be specified as an argument. The following loop
   * is completely equivalent to the former one:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * unsigned int i = 0;
   * for (auto&& [ b, a]: util::zip<1U>(thrice, twice)) {
   *
   *   a = 2 * i;
   *   b = 3.0 * i;
   *
   *   ++i;
   *
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (the index is zero-based, so `1U` refers to the second argument).
   *
   */
  template <std::size_t Lead, typename... Iterables>
  auto zip(Iterables&&... iterables);

  /// Version of `zip()` with first iterator implicitly leading the iteration.
  template <typename... Iterables>
  auto zip(Iterables&&... iterables)
  {
    return zip<0U>(std::forward<Iterables>(iterables)...);
  }

  /// @}
  // -- END -- Parallel iterations ---------------------------------------------

} // namespace util

//==============================================================================
//=== template implementation
//==============================================================================
//------------------------------------------------------------------------------
//--- util::zip()
//------------------------------------------------------------------------------
namespace util::details {

  //----------------------------------------------------------------------------
  template <std::size_t Lead, typename... Iters>
  class zip_iterator {

    static_assert(Lead < sizeof...(Iters), "The index (Lead) of the leading iterator is invalid.");

    /// Type of this object.
    using this_iterator_t = zip_iterator<Lead, Iters...>;

  public:
    // --- BEGIN -- Data types -------------------------------------------------
    /// @name Data types
    /// @{

    using difference_type = std::ptrdiff_t;
    using reference = std::tuple<typename std::iterator_traits<Iters>::reference...>;
    using value_type = std::remove_cv_t<reference>;
    using pointer = std::add_pointer_t<std::remove_reference_t<reference>>;
    using iterator_category = std::forward_iterator_tag;

    /// @}
    // --- END -- Data types ---------------------------------------------------

    // --- BEGIN -- Constructors -----------------------------------------------
    /// @name Constructors
    /// @{

    /// Constructor: copies all iterator values.
    zip_iterator(Iters&&... iterators) : fIterators(std::forward<Iters>(iterators)...) {}

    /// @}
    // --- END -- Constructors -------------------------------------------------

    // --- BEGIN -- Access -----------------------------------------------------
    /// @name Access
    /// @{

    /// Returns a tuple with values from all dereferenced iterators.
    auto operator*() const { return dereference_impl(std::index_sequence_for<Iters...>()); }

    /// Returns the iterator at the specified `Index`.
    template <std::size_t Index>
    decltype(auto) get() const
    {
      return std::get<Index>(fIterators);
    }

    /// @}
    // --- END -- Access -------------------------------------------------------

    // --- BEGIN -- Modification -----------------------------------------------
    /// @name Modification
    /// @{

    /// Increments all the iterators.
    this_iterator_t& operator++()
    {
      increment_impl(std::index_sequence_for<Iters...>());
      return *this;
    }

    /// Returns a copy of the current iterators, then increments all of the
    /// iterators in this object.
    this_iterator_t operator++(int)
    {
      this_iterator_t old(*this);
      operator++();
      return old;
    }

    /// @}
    // --- END -- Modification -------------------------------------------------

    // --- BEGIN -- Comparisons ------------------------------------------------
    /// @name Comparisons
    /// @{

    /// Comparison (based on the `Lead` iterator only).
    template <std::size_t OtherLead, typename... OtherIter>
    bool operator!=(zip_iterator<OtherLead, OtherIter...> const& other) const
    {
      return get<Lead>() != other.template get<OtherLead>();
    }

    /// Comparison (based on the `Lead` iterator only).
    template <std::size_t OtherLead, typename... OtherIter>
    bool operator==(zip_iterator<OtherLead, OtherIter...> const& other) const
    {
      return get<Lead>() == other.template get<OtherLead>();
    }

    /// @}
    // --- END -- Comparisons --------------------------------------------------

  private:
    std::tuple<Iters...> fIterators; ///< Tuple of all zipped iterators.

    /// Helper to trigger parameter pack expansion in expressions.
    template <typename... Args>
    static void expandStatements(Args&... args)
    {}

    template <std::size_t... Indices>
    void increment_impl(std::index_sequence<Indices...>)
    {
      expandStatements(++std::get<Indices>(fIterators)...);
    }

    template <std::size_t... Indices>
    auto dereference_impl(std::index_sequence<Indices...>) const
    {
      // std::tie returns a tuple of references, which is what we want.
      return std::tie(*std::get<Indices>(fIterators)...);
    }

  }; // class zip_iterator

  //----------------------------------------------------------------------------
  template <std::size_t Lead, typename... Iterables>
  auto make_zip_begin_iterator(Iterables&&... iterables)
  {
    using std::begin;
    return zip_iterator<Lead, decltype(begin(iterables))...>{begin(iterables)...};
  }

  //----------------------------------------------------------------------------
  template <std::size_t Lead, typename... Iterables>
  auto make_zip_end_iterator(Iterables&&... iterables)
  {
    using std::end;
    return zip_iterator<Lead, decltype(end(iterables))...>{end(iterables)...};
  }
} // namespace util::details

//------------------------------------------------------------------------------
template <std::size_t Lead /* = 0U */, typename... Iterables>
auto util::zip(Iterables&&... iterables)
{
  static_assert(sizeof...(Iterables) > 0, "Must zip at least one iterable.");
  return span{details::make_zip_begin_iterator<Lead>(iterables...),
              details::make_zip_end_iterator<Lead>(iterables...)};
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_ZIP_H
