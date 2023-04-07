/**
 * @file   larcorealg/Geometry/GeometryIDmapper.h
 * @brief  Mapping between geometry/readout ID and flat index.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   October 26, 2019
 * @ingroup Geometry
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYIDMAPPER_H
#define LARCOREALG_GEOMETRY_GEOMETRYIDMAPPER_H

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// C/C++ standard libraries
#include <algorithm> // std::copy()
#include <array>
#include <cassert>
#include <cstdlib> // std::size_t
#include <initializer_list>
#include <type_traits> // std::index_sequence

namespace geo {

  template <typename IDType, typename Index = std::size_t>
  class GeoIDmapper;

  template <typename Index = std::size_t>
  class TPCIDmapper;

  template <typename Index = std::size_t>
  class PlaneIDmapper;

  // ---------------------------------------------------------------------------
  namespace details {

    /// Returns a STL array of size `N` filled with `values` from the argument.
    template <std::size_t N, typename T>
    auto initializerListToArray(std::initializer_list<T> values);

  } // namespace details
  // ---------------------------------------------------------------------------

} // namespace geo

// --- BEGIN Geometry ID mappers -----------------------------------------------
/// @name Geometry ID mappers
/// @ingroup Geometry
/// @{

/**
 * @brief Class managing the mapping between geometry/readout ID and flat index.
 * @tparam IDType the geometry or readout ID to be managed
 * @tparam Index (default: `std::size_t`) type of flat index
 *
 * This class maps multi-level ID's (e.g. `geo::WireID`, `readout::TPCsetID`)
 * into a single linear index.
 *
 * The current implementation guarantees that the indices are ordered like their
 * respective IDs, and that there are no gaps so that each index up to `size()`
 * (excluded) corresponds to a "valid" ID.
 * The range of valid IDs is a hyperbox so that all elements of a certain level
 * have the same number of sub-elements (for examples, all planes have the same
 * number of wires, all valid from the mapping point of view).
 *
 * Note that this mapping is typically not suitable for channel mapping, since
 * usually planes with different wire orientations may have different number of
 * wires, making some of the indices that are valid from the point of view of
 * this ID mapping invalid in that they match wires that do not exist and should
 * not be assigned a channel number.
 */
template <typename IDType, typename Index /* = std::size_t */>
class geo::GeoIDmapper {

public:
  using ID_t = IDType;      ///< Type used as ID for this mapping.
  using index_type = Index; ///< Type of flat index.

  /**
   * @brief Default constructor: all dimensions empty.
   * @see `resize()`
   *
   * The indexer must be resized before being of any use.
   */
  GeoIDmapper();

  /**
   * @brief Prepares the indexer.
   * @param dims number of elements on all levels of the mapping
   * @see `resize()`
   *
   * The size of each dimension is specified by the corresponding number,
   * starting from the size of the outer dimension (cryostat).
   */
  GeoIDmapper(std::initializer_list<unsigned int> dims);

  // --- BEGIN Indexer status query --------------------------------------------
  /// @name Indexer status query
  /// @{

  /// Returns the number of elements in the mapping.
  index_type size() const;

  /// Returns whether the mapping has no elements (`false` by assumptions).
  bool empty() const;

  /// Dimensions of the `Level` dimension of this mapping.
  template <std::size_t Level>
  unsigned int dimSize() const;

  /// Dimensions of the ID of this mapping.
  static constexpr unsigned int dimensions();

  /// Returns whether this mapping hosts data for the specified ID.
  template <typename GeoID = ID_t>
  bool hasElement(GeoID const& id) const;

  /// Returns the ID of the first element with `GeoID` type.
  template <typename GeoID = ID_t>
  GeoID firstID() const;

  /// Returns the ID of the last covered element with `GeoID` type.
  template <typename GeoID = ID_t>
  GeoID lastID() const;

  /// @}
  // --- END Indexer status query ----------------------------------------------

  // --- BEGIN Mapping transformations -----------------------------------------
  /// @name Mapping transformations
  /// @{

  /// Returns the linear index corresponding to the specified ID.
  index_type index(ID_t const& id) const;

  /// Returns the ID corresponding to the specified linear `index`.
  ID_t ID(index_type const index) const;

  /// Returns the linear index corresponding to the specified ID.
  index_type operator()(ID_t const& id) const;

  /// Returns the ID corresponding to the specified linear `index`.
  ID_t operator()(index_type const index) const;

  /// @}
  // --- END Mapping transformations -------------------------------------------

  // --- BEGIN Mapping modification --------------------------------------------
  /// @name Mapping modification
  /// @{

  /**
   * @brief Resizes the mapping to accommodate the specified dimension sizes.
   * @param dims number of elements on all levels of the mapping
   * @see `resizeAs()`
   *
   * The size of each dimension is specified by the corresponding number,
   * starting from the size of the outer dimension (cryostat).
   */
  void resize(std::initializer_list<unsigned int> dims);

  /**
   * @brief Resizes the mapping to reflect the one from another mapping.
   * @param other ID mapping to take dimensions from
   *
   * The size of each dimension is taken by the matching one in `other`.
   */
  template <typename OIDType, typename OIndex>
  void resizeAs(GeoIDmapper<OIDType, OIndex> const& other);

  /**
   * @brief Sets all dimension sizes to `0`.
   * @see `resize()`
   *
   * The mapping needs to be resized before it is useful again.
   */
  void clear();

  /// @}
  // --- END Mapping modification ----------------------------------------------

private:
  ///< Type of dimension sizes.
  using Dimensions_t = std::array<unsigned int, dimensions()>;

  /// Number of maximum entries per ID level.
  Dimensions_t fN = zeroDimensions();

  template <std::size_t Level, typename GeoID>
  index_type indexLevel(GeoID const& id) const;

  /**
   * @brief Fills the specified ID with its index.
   * @tparam Level the level of the index to fill (`0` for cryostat level, etc._
   * @tparam GeoID type of ID to be filled
   * @param id the ID to be filled
   * @param index the index corresponding to the ID
   * @see `indexLevel()`, `index()`
   *
   * Fills the specified ID with its index. This can be considered the inverse
   * operation of the `index()` method.
   * The `index` argument is local to the level to be filled, e.g. for cryostats
   * it goes from `0` to the number of cryostats.
   */
  template <std::size_t Level, std::size_t N>
  void fillID(std::array<unsigned int, N>& id, index_type index) const;

  /// Returns whether all levels of `id` up to `Level` are within range.
  template <std::size_t Level, typename GeoID>
  bool hasElementLevel(GeoID const& id) const;

  /// Computes the expected size of this mapping.
  index_type computeSize() const;

  /// Implementation for `resizeAs()`.
  template <typename OIDType, typename OIndex, std::size_t... Indices>
  void resizeAsImpl(GeoIDmapper<OIDType, OIndex> const& other, std::index_sequence<Indices...>);

  /// Returns the number of elements at the specified `Level`.
  /// @param dimSizes the sizes of each of the levels
  template <std::size_t Level, typename Dims>
  static index_type sizeLevel(Dims const& dimSizes);

  /// Initializer with zero size for each of the dimensions.
  static Dimensions_t zeroDimensions()
  {
    Dimensions_t dims;
    dims.fill(0);
    return dims;
  }

  /// Returns whether the specified value is between `0` and the upper limit.
  template <typename Value, typename Upper>
  static bool bounded(Value v, Upper upper)
  {
    return (v >= 0) && (static_cast<index_type>(v) < upper);
  }

}; // class geo::GeoIDmapper<>

/** ****************************************************************************
 * @brief Mapping for TPC identifiers.
 * @tparam Index (default: `std::size_t`) type of flat index
 * @see `geo::GeoIDmapper`
 *
 * A customized version of `geo::GeoIDmapper` offering TPC ID-specific
 * interface.
 */
template <typename Index /* = std::size_t */>
class geo::TPCIDmapper : public geo::GeoIDmapper<geo::TPCID, Index> {

  /// Base class.
  using BaseMapper_t = GeoIDmapper<TPCID, Index>;

public:
  // import types
  using ID_t = typename BaseMapper_t::ID_t;
  using index_type = typename BaseMapper_t::index_type;

  // import all constructors from `geo::GeoIDmapper`
  using BaseMapper_t::BaseMapper_t;

  /**
   * @brief Prepares the mapping with the specified sizes.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs per cryostat
   *
   * The mapping is sized to map `nCryo` cryostats, each with `nTPCs` TPCs.
   */
  TPCIDmapper(unsigned int nCryo, unsigned int nTPCs) : BaseMapper_t({nCryo, nTPCs}) {}

  // --- BEGIN Mapping modification --------------------------------------------
  /// @name Mapping modification
  /// @{

  using BaseMapper_t::resize;

  /**
   * @brief Prepares the mapping for the specified sizes.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @see `resizeAs()`
   *
   * The mapping is sized to map `nCryo` cryostats, each with `nTPCs` TPCs.
   */
  void resize(unsigned int nCryo, unsigned int nTPCs) { BaseMapper_t::resize({nCryo, nTPCs}); }

  /// @}
  // --- END Mapping modification ----------------------------------------------

  // --- BEGIN Mapping status query --------------------------------------------
  /// @name Mapping status query
  /// @{

  /// Returns whether this mapping covers the specified cryostat.
  bool hasCryostat(CryostatID const& cryoid) const { return BaseMapper_t::hasElement(cryoid); }

  /// Returns whether this mapping covers the specified TPC.
  bool hasTPC(TPCID const& tpcid) const { return BaseMapper_t::hasElement(tpcid); }

  /// @}
  // --- END Mapping status query ----------------------------------------------

}; // geo::TPCIDmapper<>

/** ****************************************************************************
 * @brief Mapping for sensitive plane identifiers.
 * @tparam Index (default: `std::size_t`) type of flat index
 * @see `geo::GeoIDmapper`
 *
 * A customized version of `geo::GeoIDmapper` offering
 * sensitive plane ID-specific interface.
 */
template <typename Index /* = std::size_t */>
class geo::PlaneIDmapper : public geo::GeoIDmapper<geo::PlaneID, Index> {

  /// Base class.
  using BaseMapper_t = GeoIDmapper<PlaneID, Index>;

public:
  // import types
  using ID_t = typename BaseMapper_t::ID_t;
  using index_type = typename BaseMapper_t::index_type;

  // import all constructors from `geo::GeoIDmapper`
  using BaseMapper_t::BaseMapper_t;

  /**
   * @brief Prepares the mapping with the specified sizes.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs per cryostat
   * @param nPlanes number of planes per TPC
   *
   * The mapping is sized to map `nCryo` cryostats, each with `nTPCs` TPCs,
   * each one with `nPlanes` wire planes.
   */
  PlaneIDmapper(unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes)
    : BaseMapper_t({nCryo, nTPCs, nPlanes})
  {}

  // --- BEGIN Mapping modification --------------------------------------------
  /// @name Mapping modification
  /// @{

  using BaseMapper_t::resize;

  /**
   * @brief Prepares the mapping for the specified sizes.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @param nPlanes number of planes per TPC
   * @see `resizeAs()`
   *
   * The mapping is sized to map `nCryo` cryostats, each with `nTPCs` TPCs,
   * each one with `nPlanes` wire planes.
   */
  void resize(unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes)
  {
    BaseMapper_t::resize({nCryo, nTPCs, nPlanes});
  }

  /// @}
  // --- END Mapping modification ----------------------------------------------

  // --- BEGIN Mapping status query --------------------------------------------
  /// @name Mapping status query
  /// @{

  /// Returns whether this mapping covers the specified cryostat.
  bool hasCryostat(CryostatID const& cryoid) const { return BaseMapper_t::hasElement(cryoid); }

  /// Returns whether this mapping covers the specified TPC.
  bool hasTPC(TPCID const& tpcid) const { return BaseMapper_t::hasElement(tpcid); }

  /// Returns whether this mapping covers the specified plane.
  bool hasPlane(PlaneID const& planeid) const { return BaseMapper_t::hasElement(planeid); }

  /// @}
  // --- END Mapping status query ----------------------------------------------

}; // geo::PlaneIDmapper<>

/// @}
// --- END Geometry ID mappers -------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--- Template implementation
//------------------------------------------------------------------------------
template <std::size_t N, typename T>
auto geo::details::initializerListToArray(std::initializer_list<T> values)
{
  std::array<T, N> data;
  std::copy(values.begin(), values.end(), data.begin());
  return data;
}

//------------------------------------------------------------------------------
//--- geo::GeoIDmapper
//------------------------------------------------------------------------------
template <typename IDType, typename Index>
geo::GeoIDmapper<IDType, Index>::GeoIDmapper()
{
  fN.fill(0U);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
geo::GeoIDmapper<IDType, Index>::GeoIDmapper(std::initializer_list<unsigned int> dims)
  : fN(details::initializerListToArray<dimensions()>(dims))
{
  assert(dims.size() == dimensions()); // can't be static
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
auto geo::GeoIDmapper<IDType, Index>::size() const -> index_type
{
  return computeSize();
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
bool geo::GeoIDmapper<IDType, Index>::empty() const
{
  return size() == index_type{0};
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <std::size_t Level>
unsigned int geo::GeoIDmapper<IDType, Index>::dimSize() const
{
  if constexpr (Level >= dimensions())
    return 0U; // technically it would be 1...
  else
    return fN[Level];
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
constexpr unsigned int geo::GeoIDmapper<IDType, Index>::dimensions()
{
  return IDType::Level + 1;
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <typename GeoID>
bool geo::GeoIDmapper<IDType, Index>::hasElement(GeoID const& id) const
{
  return hasElementLevel<GeoID::Level>(id);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <typename GeoID /* = ID_t */>
GeoID geo::GeoIDmapper<IDType, Index>::firstID() const
{
  if constexpr (GeoID::Level == 0)
    return GeoID(0U);
  else
    return GeoID(firstID<typename GeoID::ParentID_t>(), 0U);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <typename GeoID /* = ID_t */>
GeoID geo::GeoIDmapper<IDType, Index>::lastID() const
{
  if constexpr (GeoID::Level == 0)
    return GeoID(fN[GeoID::Level] - 1U);
  else
    return GeoID(lastID<typename GeoID::ParentID_t>(), fN[GeoID::Level] - 1U);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
auto geo::GeoIDmapper<IDType, Index>::index(ID_t const& id) const -> index_type
{
  return indexLevel<ID_t::Level>(id);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
auto geo::GeoIDmapper<IDType, Index>::ID(index_type const index) const -> ID_t
{
  std::array<unsigned int, ID_t::Level + 1u> numbers;
  fillID<ID_t::Level>(numbers, index);
  auto id = std::make_from_tuple<ID_t>(numbers);
  id.setValidity(id.Cryostat < fN[0u]);
  return id;
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
auto geo::GeoIDmapper<IDType, Index>::operator()(ID_t const& id) const -> index_type
{
  return index(id);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
auto geo::GeoIDmapper<IDType, Index>::operator()(index_type const index) const -> ID_t
{
  return ID(index);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
void geo::GeoIDmapper<IDType, Index>::resize(std::initializer_list<unsigned int> dims)
{
  fN = details::initializerListToArray<dimensions()>(dims);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <typename OIDType, typename OIndex>
void geo::GeoIDmapper<IDType, Index>::resizeAs(GeoIDmapper<OIDType, OIndex> const& other)
{
  resizeAsImpl(other, std::make_index_sequence<dimensions()>{});
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
void geo::GeoIDmapper<IDType, Index>::clear()
{
  fN.fill(0U);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <std::size_t Level, typename GeoID>
auto geo::GeoIDmapper<IDType, Index>::indexLevel(GeoID const& id) const -> index_type
{
  if constexpr (Level == 0)
    return id.template getIndex<0U>();
  else {
    return indexLevel<(Level - 1U)>(id) * fN[Level] + id.template getIndex<Level>();
  }
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <std::size_t Level, std::size_t N>
void geo::GeoIDmapper<IDType, Index>::fillID(std::array<unsigned int, N>& numbers,
                                             index_type index) const
{
  if constexpr (Level == 0) { numbers[0] = index; }
  else {
    numbers[Level] = index % fN[Level];
    fillID<Level - 1U>(numbers, index / fN[Level]);
  }
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <std::size_t Level, typename GeoID>
bool geo::GeoIDmapper<IDType, Index>::hasElementLevel(GeoID const& id) const
{
  if (!bounded(id.template getIndex<Level>(), fN[Level])) return false;
  if constexpr (Level == 0U)
    return true;
  else
    return hasElementLevel<(Level - 1U)>(id);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
auto geo::GeoIDmapper<IDType, Index>::computeSize() const -> index_type
{
  return sizeLevel<0U>(fN);
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <typename OIDType, typename OIndex, std::size_t... Indices>
void geo::GeoIDmapper<IDType, Index>::resizeAsImpl(GeoIDmapper<OIDType, OIndex> const& other,
                                                   std::index_sequence<Indices...>)
{
  // Clang 5.0.1 does not understand `other.dimensions()` is constexpr
  static_assert(GeoIDmapper<OIDType, OIndex>::dimensions() >= dimensions(),
                "Can't resize a deeper mapping to a shallower one.");
  resize({other.template dimSize<Indices>()...});
}

//------------------------------------------------------------------------------
template <typename IDType, typename Index>
template <std::size_t Level, typename Dims>
auto geo::GeoIDmapper<IDType, Index>::sizeLevel(Dims const& dimSizes) -> index_type
{
  if constexpr (Level >= dimensions())
    return 1U;
  else
    return sizeLevel<(Level + 1U)>(dimSizes) * dimSizes[Level];
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_GEOMETRYIDMAPPER_H
