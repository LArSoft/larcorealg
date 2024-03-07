/**
 * @file   larcorealg/CoreUtils/span.h
 * @brief  An object with a begin and end iterator.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 *
 * This library is header only.
 */

#ifndef LARCOREALG_COREUTILS_SPAN_H
#define LARCOREALG_COREUTILS_SPAN_H

// Boost libraries
#include "boost/iterator/transform_iterator.hpp"

// C/C++ standard library
#include <iterator>    // std::begin(), std::end()
#include <type_traits> // std::decay_t<>, std::declval(), std::invoke_result_t<>
#include <utility>     // std::as_const()

namespace util {

  /// Non-templated base class for `span`.
  struct span_base {

    /// Returns the begin iterator of the specified container.
    template <typename Cont>
    static decltype(auto) get_begin(Cont& cont)
    {
      using std::begin;
      return begin(cont);
    }

    /// Returns the end iterator of the specified container.
    template <typename Cont>
    static decltype(auto) get_end(Cont& cont)
    {
      using std::end;
      return end(cont);
    }

    /// Type of begin iterator of `Cont` type.
    template <typename Cont>
    using get_begin_iterator = std::decay_t<decltype(get_begin(std::declval<Cont>()))>;

    /// Type of end iterator of `Cont` type.
    template <typename Cont>
    using get_end_iterator = std::decay_t<decltype(get_end(std::declval<Cont>()))>;

    /// Returns the constant begin iterator of the specified container.
    template <typename Cont>
    static decltype(auto) get_cbegin(Cont& cont)
    {
      using std::cbegin;
      return cbegin(cont);
    }

    /// Returns the constant end iterator of the specified container.
    template <typename Cont>
    static decltype(auto) get_cend(Cont& cont)
    {
      using std::cend;
      return cend(cont);
    }

    /// Type of constant begin iterator of `Cont` type.
    template <typename Cont>
    using get_cbegin_iterator = std::decay_t<decltype(get_cbegin(std::declval<Cont>()))>;

    /// Type of constant end iterator of `Cont` type.
    template <typename Cont>
    using get_cend_iterator = std::decay_t<decltype(get_cend(std::declval<Cont>()))>;

  }; // struct span_base

  /**
   * @brief Simple class with a begin and an end.
   * @tparam BITER type of begin iterator
   * @tparam EITER type of end iterator (default: as `BITER`)
   *
   * This class is a glorified pair of iterators which can be used in a
   * range-for loop.
   * All input iterators are accepted.
   *
   * It is probably going to be redundant with the advent of C++20 and its span
   * and/or range libraries (if the latter is ever going to happen).
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using Data_t = std::vector<int>;
   *
   * int sum(util::span<Data_t::iterator> r) {
   *   decltype(r)::value_type s = 0;
   *   for (auto v: r) s += v;
   *   return s;
   * } // sum
   *
   * void analyse() {
   *
   *   Data_t v(10U);
   *   std::iota(v.begin(), v.end(), 1); // { 1, 2, 3, 4 ,5 ,6, 7, 8, 9, 10 }
   *
   *   auto span5 = util::span(v.begin(), v.begin() + 5);
   *   std::cout << "Sum of 5 numbers: " << sum(span5) << std::endl;
   *
   *   auto span8 = util::span(v.begin(), v.begin() + 8);
   *   std::cout << "Sum of 8 numbers: " << sum(span8) << std::endl;
   *
   *   auto full_span = util::make_span(v);
   *   std::cout << "Sum of all numbers: " << sum(full_span) << std::endl;
   *
   * } // analyse()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The call to `analyse()` will produce output like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Sum of 5 numbers: 15
   * Sum of 8 numbers: 36
   * Sum of all numbers: 55
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The convenience of this class, not evident from the example, is the
   * beginning and end of the span can be passed around as a single object,
   * with the added bonus of a (very reduced) collection interface.
   *
   * @note The constantness of the access is completely determined by the type
   *       of the iterators. A `const span` object does _not_ guarantee
   *       unmutable access to the iterated data.
   */
  template <typename BIter, typename EIter = BIter>
  struct span : public span_base, private std::pair<BIter, EIter> {

    using begin_iterator = BIter; ///< Type of begin iterator.
    using end_iterator = EIter;   ///< Type of end iterator.

    /// Type of this class.
    using span_t = span<begin_iterator, end_iterator>;

    /// Type of iterator pair.
    using pair_t = std::pair<begin_iterator, end_iterator>;

    /// Type of values pointed by the iterators.
    using value_type = typename begin_iterator::value_type;

    /// Type of reference pointed by the iterators.
    using reference = typename begin_iterator::reference;

    /// Constructor: specifies the begin and end iterator.
    span(begin_iterator b, end_iterator e) : pair_t(b, e) {}

    /// Constructor: specifies the begin and end iterator and an adapter.
    template <typename SrcIterB, typename SrcIterE, typename Adaptor>
    span(SrcIterB&& b, SrcIterE&& e, Adaptor&& adaptor)
      : span(adaptor(std::forward<SrcIterB>(b)), adaptor(std::forward<SrcIterE>(e)))
    {}

    /// Constructor: copies from another span (possibly with different types).
    template <typename OBIter, typename OEIter>
    span(span<OBIter, OEIter> const& from) : span(from.begin(), from.end())
    {}

    // C++ boilerplate
    span(span_t const&) = default;
    span(span_t&&) = default;
    span_t& operator=(span_t const&) = default;
    span_t& operator=(span_t&&) = default;

    /// Returns a copy of the begin iterator.
    begin_iterator begin() const { return pair_t::first; }

    /// Returns a copy of the end iterator.
    end_iterator end() const { return pair_t::second; }

    /// Returns the size between begin and end, converted to `std::size_t`.
    std::size_t size() const { return std::distance(begin(), end()); }

    /// Returns whether the span is empty (that is, no steps between them).
    bool empty() const { return begin() == end(); }

  }; // span

  // deduction guide for adapted span
  template <typename IterB, typename IterE, typename Adaptor>
  span(IterB&& b, IterE&& e, Adaptor&& adaptor)
    -> span<std::invoke_result_t<Adaptor, IterB>, std::invoke_result_t<Adaptor, IterE>>;

  // --- BEGIN -- Span helper functions ----------------------------------------
  /// @name Span helper functions
  /// @{

  /// Creates a span from specified iterators (can use constructor instead).
  template <typename BIter, typename EIter>
  auto make_span(BIter begin, EIter end)
  {
    return span{begin, end};
  }

  /// Creates a span from a container type.
  template <typename Cont>
  auto make_span(Cont& cont)
  {
    return span{span_base::get_begin(cont), span_base::get_end(cont)};
  }

  /// Creates a span with constant iterator access from a container type.
  template <typename Cont>
  auto make_const_span(Cont& cont)
  {
    return span{span_base::get_cbegin(cont), span_base::get_cend(cont)};
  }

  /// @}
  // --- END -- Span helper functions ------------------------------------------

  // --- BEGIN -- Adapted span helper functions --------------------------------
  /// @name Adapted span helper functions
  /// @{

  /// Creates a span from specified iterators via an adaptor.
  /// @see `make_adapted_span(Cont&, Adaptor&&)`
  template <typename BIter, typename EIter, typename Adaptor>
  auto make_adapted_span(BIter begin, EIter end, Adaptor&& adaptor)
  {
    return util::span(adaptor(begin), adaptor(end));
  }

  /**
   * @brief Creates a span from specified collection via an adaptor.
   * @param cont collection to be iterated through
   * @param adapter iterator transformation to be applied on `cont`
   * @return a `util::span` object iterating via an adapter to `cont`
   * @see `make_adapted_span(BIter, EIter, Adaptor&&)`,
   *      `make_adapted_const_span(Cont, Adaptor&&)`
   *
   * This interface is just a little help on using _iterator_ adapters.
   * Note that `adapter` transforms the iterators of `cont`, not its values.
   * The adapter needs to be written aside and it's in general not trivial.
   * An example of iterating through a collection of objects (actually, just
   * `float` for simplicity) stored as unique pointers in a collection:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * float accumulate(std::vector<std::unique_ptr<float>> const& v) {
   *
   *   using src_iterator = std::vector<std::unique_ptr<float>>::const_iterator;
   *
   *   float sum = 0.0F;
   *   for (float v: util::make_adapted_span(v, boost::make_indirect_iterator<src_iterator>))
   *     sum += v;
   *
   *   return sum;
   * } // accumulate()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * This example shows the usage of `util::make_adapted_span()`. In the
   * specific example, it would have been more explicit to use the constant
   * counterpart, `util::make_adapted_const_span()`.
   * The adaptor helper `boost::make_indirect_iterator()` is provided in
   * `boost/iterator/indirect_iterator.hpp`.
   */
  template <typename Cont, typename Adaptor>
  auto make_adapted_span(Cont& cont, Adaptor&& adaptor)
  {
    return make_adapted_span(
      span_base::get_begin(cont), span_base::get_end(cont), std::forward<Adaptor>(adaptor));
  }

  /// Creates constant iteration span from specified collection via an adaptor.
  // @see `make_adapted_span(Cont, Adaptor&&)`
  template <typename Cont, typename Adaptor>
  auto make_adapted_const_span(Cont& cont, Adaptor&& adaptor)
  {
    return make_adapted_span(
      span_base::get_cbegin(cont), span_base::get_cend(cont), std::forward<Adaptor>(adaptor));
  }

  /// @}
  // --- END -- Adapted span helper functions ----------------------------------

  // --- BEGIN -- Transformed span helper functions ----------------------------
  /**
   * @name Transformed span helper functions
   *
   *
   *
   */
  /// @{

  /// Creates a span from specified iterators via an adaptor.
  /// @see `make_transformed_span(Cont&, Op&&)`
  template <typename BIter, typename EIter, typename Op>
  auto make_transformed_span(BIter begin, EIter end, Op&& op)
  {
    auto adaptor = [&op](auto iter) { return boost::make_transform_iterator(iter, op); };
    return util::make_adapted_span(begin, end, adaptor);
  }

  /**
   * @brief Creates a span from specified collection via an adaptor.
   * @param cont collection to be iterated through
   * @param adapter iterator transformation to be applied on `cont`
   * @return a `util::span` object iterating via an adapter to `cont`
   * @see `make_transformed_span(BIter, EIter, Op&&)`,
   *      `make_transformed_const_span(Cont, Op&&)`,
   *      `make_adapted_span()`,
   *
   * This function transforms as directed by the unary operation `op` the result
   * of each iteration cycle, so that instead of `cont[i]`, a `op(cont[i])`
   * is assigned to the iteration control variable.
   * An example of iterating through a collection of objects (actually, just
   * `float` for simplicity) stored as unique pointers in a collection:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * float accumulate(std::vector<std::unique_ptr<float>> const& v) {
   *
   *   float sum = 0.0F;
   *   for (float v: util::make_transformed_span(v, [](auto& ptr){ return *ptr; }))
   *     sum += v;
   *
   *   return sum;
   * } // accumulate()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * This example shows the usage of `util::make_transformed_span()`. In the
   * specific example, it would have been more explicit to use the constant
   * counterpart, `util::make_transformed_const_span()`.
   * Note that in the example the dereference result is _copied_ into `v`.
   * For have it by reference, the operation `op` must be properly written:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void scale(std::vector<std::unique_ptr<float>>& v, float factor) {
   *
   *   for (float& v: util::make_transformed_span(v, [](auto& ptr) -> float& { return *ptr; }))
   *     v *= factor;
   *
   * } // scale()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename Cont, typename Op>
  auto make_transformed_span(Cont& cont, Op&& op)
  {
    return make_transformed_span(
      span_base::get_begin(cont), span_base::get_end(cont), std::forward<Op>(op));
  }

  /// Creates constant iteration span from specified collection via a
  /// transformation `op`.
  // @see `make_transformed_span(Cont, Op&&)`
  template <typename Cont, typename Op>
  auto make_transformed_const_span(Cont& cont, Op&& op)
  {
    return make_transformed_span(std::as_const(cont), std::forward<Op>(op));
  }

  /// @}
  // --- END -- Adapted span helper functions ----------------------------------

} // namespace util

#endif // LARCOREALG_COREUTILS_SPAN_H
