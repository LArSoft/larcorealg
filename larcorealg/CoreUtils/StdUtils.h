/**
 * @file   larcorealg/CoreUtils/StdUtils.h
 * @brief  Functions pulling in STL customization if available.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July , 2019
 */

#ifndef LARCOREALG_COREUTILS_STDUTILS_H
#define LARCOREALG_COREUTILS_STDUTILS_H

// LArSoft libraries
#include "larcorealg/CoreUtils/MetaUtils.h" // util::is_basic_string_type_v ...

// C/C++ standard libraries
#include <iterator> // std::begin(), std::end(), ...
#include <string>   // std::to_string()
#include <utility>  // std::get()

namespace util {

  /**
   * @name C++ standard library customization for user-defined classes.
   * @defgroup LArSoft_CoreUtils_StdUtils C++ STL customizations
   *
   * There are a number of functions that are provided by C++ standard library
   * for the data types and classes defined in the standard.
   * It is often desirable to have your class react to these standard functions
   * in a standard way, for example for a container to react to `std::begin()`
   * to return its `begin()` iterator. While sometimes this is easy (for example
   * `std::begin()` calls `begin()` member function if available), some other
   * times that is not possible. In that case, since overloading of functions in
   * the `std` namespace is not allowed by C++, the usual pattern is to rely on
   * the argument-dependent lookup (known as "ADL") to have the comnpiler find
   * the overloaded function that is defined in the same namespace as any of
   * the arguments. For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * userns::MyObj obj;
   * // ...
   * using std::to_string;
   * std::string objstr = to_string(obj);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will look in the namespace where the type of `obj` is defined (that is
   * `userns`) for a `userns::to_string`, then will consider `std::to_string`
   * it self.
   *
   * The utilities provided here provide a transparent way to do that, at the
   * cost of a new header and some non-standard call. The equivalent call of the
   * above would be:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * userns::MyObj obj;
   * // ...
   * std::string objstr = util::to_string(obj);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * @note For customization of templates, like `std::hash` or
   *       `std::numeric_limits`, specialization of classes in `std` is allowed
   *       by the standard, so no particular trick is required.
   */
  /// @{

  /// ADL-aware version of `std::to_string`.
  template <typename T>
  constexpr decltype(auto) to_string(T&& obj);
  //     { using std::to_string; return to_string(std::forward<T>(obj)); }

  // --- BEGIN --- Containers and iterators ------------------------------------
  /// ADL-aware version of `std::begin`.
  template <typename T>
  constexpr decltype(auto) begin(T&& obj)
  {
    using std::begin;
    return begin(std::forward<T>(obj));
  }

  /// ADL-aware version of `std::end`.
  template <typename T>
  constexpr decltype(auto) end(T&& obj)
  {
    using std::end;
    return end(std::forward<T>(obj));
  }

  /// ADL-aware version of `std::cbegin`.
  template <typename T>
  constexpr decltype(auto) cbegin(T&& obj)
  {
    using std::cbegin;
    return cbegin(std::forward<T>(obj));
  }

  /// ADL-aware version of `std::cend`.
  template <typename T>
  constexpr decltype(auto) cend(T&& obj)
  {
    using std::cend;
    return cend(std::forward<T>(obj));
  }

  /// ADL-aware version of `std::size`.
  template <typename T>
  constexpr decltype(auto) size(T&& obj)
  {
    using std::size;
    return size(std::forward<T>(obj));
  }

  /// ADL-aware version of `std::empty`.
  template <typename T>
  constexpr decltype(auto) empty(T&& obj)
  {
    using std::empty;
    return empty(std::forward<T>(obj));
  }

  // --- END --- Containers and iterators --------------------------------------

  // --- BEGIN --- tuples ------------------------------------------------------

  template <std::size_t I, typename T>
  decltype(auto) get(T&& obj)
  {
    using std::get;
    return get<I>(std::forward<T>(obj));
  }

  // --- END --- tuples --------------------------------------------------------

  /// @}

} // namespace util

// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
namespace util::details {

  // ---------------------------------------------------------------------------
  template <typename T, typename = void>
  struct ToStringImpl {

    template <typename U>
    static std::string to_string(U&& obj)
    {
      using std::to_string;
      return to_string(std::forward<U>(obj));
    }

  }; // struct ToStringImpl

  // ---------------------------------------------------------------------------
  template <typename T>
  struct ToStringImpl<T, std::enable_if_t<util::is_basic_string_type_v<T>>> {

    template <typename U>
    static std::string to_string(U&& obj)
    {
      return obj;
    }

  }; // struct ToStringImpl<string>

  // ---------------------------------------------------------------------------
  template <typename T>
  struct ToStringImpl<T, std::enable_if_t<util::is_basic_string_view_type_v<T>>> {

    template <typename U>
    static std::string to_string(U&& obj)
    {
      return {obj.begin(), obj.end()};
    }

  }; // struct ToStringImpl<string_view>

  // ---------------------------------------------------------------------------

} // namespace util::details

// -----------------------------------------------------------------------------
template <typename T>
constexpr decltype(auto) util::to_string(T&& obj)
{
  return util::details::ToStringImpl<T>::to_string(std::forward<T>(obj));
}
//   { using std::to_string; return to_string(std::forward<T>(obj)); }

// -----------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_STDUTILS_H
