/**
 * @file   larcorealg/CoreUtils/operations.h
 * @brief  Provides a few simple operations for use in generic programming.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 5, 2019
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_COREUTILS_OPERATIONS_H
#define LARCOREALG_COREUTILS_OPERATIONS_H

// C++ standard library
#include <memory> // std::addressof()

namespace util {

  //----------------------------------------------------------------------------
  /**
   * @brief Functor returning the address in memory of the operand.
   * @see `util::takeAddress()`
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int*> ptrs(data.size());
   * std::transform
   *   (data.begin(), data.end(), ptrs.begin(), util::AddressTaker{});
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will fill the vector `ptrs` with pointers to the elements of `data`.
   *
   * @note The address is extracted via `std::addressof()` and it bypasses the
   *       `operator&()` of the operand.
   */
  struct AddressTaker {

    /// Returns the address of the argument.
    template <typename T>
    auto operator()(T& ref) const
    {
      return std::addressof(ref);
    }

  }; // struct AddressTaker

  /**
   * @brief Returns a functor that returns the address of its argument.
   * @see `util::AddressTaker`
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int*> ptrs(data.size());
   * std::transform
   *   (data.begin(), data.end(), ptrs.begin(), util::takeAddress());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will fill the vector `ptrs` with pointers to the elements of `data`.
   *
   * Why bother?
   * ------------
   *
   * C++ already provides a tool to effectively take an address,
   * `std::addressof`. The reason for `takeAddress()` is that `std::addressof()`
   * is a function, with many overloads, and to use it in a STL algorithm the
   * overload has to be resolved first. For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using addressof_t = int const*(*)(int const&);
   *
   * std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtr),
   *   ((addressof_t) &std::addressof));
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * One important limit is that the type of the argument (in this case
   * `int const&`) needs to be known or deduced in a quite precise way, in
   * particular regarding constantness and referenceness.
   * This is unconvenient and convoluted enough that one would rather create
   * a new function, like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto takeAddress = [](auto&& ref){ return std::addressof(ref); };
   *
   * std::vector<int const*> dataPtr;
   * std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtr),
   *   takeAddress);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * This `util::takeAddress()` operates in a very similar way to the lambda in
   * the last example.
   */
  decltype(auto) takeAddress()
  {
    return AddressTaker();
  }

  //----------------------------------------------------------------------------
  /**
   * @brief Functor dereferencing the operand.
   * @see `util::dereference()`
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> values(ptrs.size());
   * std::transform
   *   (ptrs.cbegin(), ptrs.cend(), values.begin(), util::Dereferencer{});
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will fill the vector `values` with the values pointed by the elements in
   * `ptrs`.
   */
  struct Dereferencer {

    /// Returns `*ptr`.
    template <typename T>
    decltype(auto) operator()(T&& ptr) const
    {
      return *ptr;
    }

  }; // struct Dereferencer

  /**
   * @brief Returns a functor that returns `*ptr` of its argument `ptr`.
   * @see `util::Dereferencer`
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> values(ptrs.size());
   * std::transform
   *   (ptrs.cbegin(), ptrs.cend(), values.begin(), util::dereference());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will fill the vector `values` with the values pointed by the elements in
   * `ptrs`.
   */
  decltype(auto) dereference()
  {
    return Dereferencer();
  }

  //----------------------------------------------------------------------------

} // namespace util

#endif // LARCOREALG_COREUTILS_OPERATIONS_H
