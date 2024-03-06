/**
 * @file   larcorealg/CoreUtils/makeValueIndex.h
 * @brief  Provides `util::makeValueIndex()` helper function.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 24, 2019
 *
 * This is a header-only library.
 *
 * This utility belongs to the same cetegory as `util::MakeIndex()`, but that
 * one is definded in `lardataalg` and we need this one for geometry.
 */

#ifndef LARCOREALG_COREUTILS_MAKEVALUEINDEX_H
#define LARCOREALG_COREUTILS_MAKEVALUEINDEX_H

// LArSoft libraries
#include "larcorealg/CoreUtils/DebugUtils.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/fromFutureImport.h" // util::pre_std::identity<>

// C/C++ standard library
#include <cstddef> // std::size_t
#include <map>
#include <stdexcept>   // std::runtime_error
#include <string>      // std::to_string()
#include <type_traits> // std::invoke_result_t, std::remove_reference_t

namespace util {

  /**
   * @brief Returns a map of value to index.
   * @tparam Coll type of container
   * @tparam Extractor type of value extraction function
   * @param coll container to get the map of
   * @param getter function applied to each element to extract the value for map
   * @return a value-to-index associative container
   * @throw std::runtime_error if multiple elements yield the same value
   *
   * The collection `coll` is navigated in sequence from `begin()` to `end()`,
   * and a map is created where to each key, `getter(coll[i])`, the index `i`
   * is associated. The value returned by `getter()` is copied into the key.
   * Therefore that value needs to satisfy all the requirements of the key of
   * the STL associative container `std::map`.
   * Duplicate values will trigger an exception.
   *
   * Requirements
   * -------------
   *
   * The collection type `Coll` must have:
   *  * `value_type` type defining the content of the container
   *  * support `begin()` and `end()` free functions returning input iterators,
   *    that is support a ranged-for loop
   *
   */
  template <typename Coll, typename Extractor>
  decltype(auto) makeValueIndex(Coll const& coll, Extractor getter);

  template <typename Coll>
  auto makeValueIndex(Coll const& coll)
  {
    return makeValueIndex(coll, util::pre_std::identity());
  }

} // namespace util

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Coll, typename Extractor>
decltype(auto) util::makeValueIndex(Coll const& coll, Extractor getter)
{

  using Value_t = typename Coll::value_type;
  using Key_t
#if 0  // this is C++17...
    = std::remove_reference_t<std::invoke_result_t<Extractor, Value_t>>;
#else  // ... and this is what Clang 5.0 understands:
    = std::remove_reference_t<decltype(getter(std::declval<Value_t>()))>;
#endif // 0

    using Map_t = std::map<Key_t, std::size_t>;

  Map_t index;
  for (auto&& [iValue, collValue] : util::enumerate(coll)) {

    Key_t const& key = getter(collValue);
    auto const iKey = index.lower_bound(key);
    if ((iKey != index.end()) && (iKey->first == key)) {
      // no guarantee that `key` supports `std::to_string()`: print only indices
      throw std::runtime_error(std::string(__func__) + ": element #" + std::to_string(iValue) +
                               " has the same key as #" + std::to_string(iKey->second));
    }
    index.emplace_hint(iKey, key, iValue);
  } // for

  return index;

} // util::makeValueIndex()

//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_MAKEVALUEINDEX_H
