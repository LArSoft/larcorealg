/**
 * @file   larcorealg/Geometry/GeometryDataContainers.h
 * @brief  Containers to hold one datum per TPC or plane.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 2nd, 2018
 * @ingroup Geometry
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYDATACONTAINERS_H
#define LARCOREALG_GEOMETRY_GEOMETRYDATACONTAINERS_H

// LArSoft libraries
#include "larcorealg/CoreUtils/span.h"
#include "larcorealg/Geometry/GeometryIDmapper.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Boost libraries
#include <boost/iterator/iterator_adaptor.hpp>

// C/C++ standard libraries
#include <algorithm> // std::fill(), std::for_each()
#include <cassert>
#include <initializer_list>
#include <stdexcept> // std::out_of_range
#include <string>
#include <utility> // std::forward()
#include <vector>

namespace geo {

  template <typename T, typename Mapper>
  class GeoIDdataContainer;

  template <typename T>
  class TPCDataContainer;

  template <typename T>
  class PlaneDataContainer;

  // ---------------------------------------------------------------------------
  namespace details {

    template <typename T>
    class GeoContainerData;

    template <typename GeoIDdataContainerClass, typename BaseIterator>
    class GeoIDdataContainerIterator;

    template <typename GeoIDIteratorClass>
    class GeoIDdataContainerItemIterator;

  } // namespace details
  // ---------------------------------------------------------------------------

} // namespace geo

// --- BEGIN Geometry data containers ----------------------------------------

/// @name Geometry data containers
/// @ingroup Geometry
/// @{

/** **************************************************************************
 * @brief Container with one element per geometry TPC.
 * @tparam T type of the contained datum
 * @see `geo::GeometryCore::makeTPCData`
 *
 * The container is of fixed size and can't be neither resized nor freed
 * before destruction.
 *
 * This example creates a "map" of tracks starting on each TPC:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto const* geom = lar::providerFrom<geo::GeometryCore>();
 * geo::TPCDataContainer<std::vector<recob::Track const*>> TracksPerTPC
 *   (geom->NCryostats(), geom->MaxTPCs());
 *
 * for (recob::Track const& track: tracks) {
 *   geo::TPCGeo const* tpc = geom->PositionToTPCptr(track.Start());
 *   if (tpc) TracksPerTPC[tpc->ID()].push_back(tpc);
 * } // for
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 * Assumptions
 * ============
 *
 * The following assumptions should be considered unchecked, and the behavior
 * when they are violated undefined (but note that in debug mode some of them
 * might be actually checked):
 * * the container assumes the same number of TPCs in each cryostat. While
 *   this is not effectively a necessary condition, keep in mind that this
 *   container has no notion whether a given TPC actually exists in the
 *   geometry or not
 * * at least one element is expected to be present
 *
 */
template <typename T, typename Mapper>
class geo::GeoIDdataContainer {

  using This_t = GeoIDdataContainer<T, Mapper>; ///< Type of this class.

  /// Type of data container helper.
  using Container_t = details::GeoContainerData<T>;

  /// Type of iterator to the data.
  using BaseIter_t = typename Container_t::iterator;

  /// Type of constant iterator to the data.
  using BaseConstIter_t = typename Container_t::const_iterator;

  /// Functor to extract an ID data member.
  struct IDextractor {
    template <typename Obj>
    decltype(auto) operator()(Obj&& obj) const
    {
      return obj.ID();
    }
  }; // struct IDextractor

public:
  /// Type of mapper between IDs and index.
  using Mapper_t = Mapper;

  using ID_t = typename Mapper_t::ID_t; ///< Type used as ID for this container.

  /// @{
  /// @name STL container types.

  using value_type = typename Container_t::value_type;
  using reference = typename Container_t::reference;
  using const_reference = typename Container_t::const_reference;
  using pointer = typename Container_t::pointer;
  using const_pointer = typename Container_t::const_pointer;
  using iterator = details::GeoIDdataContainerIterator<Mapper_t, BaseIter_t>;
  using const_iterator = details::GeoIDdataContainerIterator<Mapper_t, BaseConstIter_t>;
  using difference_type = typename Container_t::difference_type;
  using size_type = typename Container_t::size_type;

  /// Special iterator dereferencing to pairs ( ID, value ) (see `items()`).
  using item_iterator = details::GeoIDdataContainerItemIterator<iterator>;

  /// Special iterator dereferencing to pairs ( ID, value ) (see `items()`).
  using item_const_iterator = details::GeoIDdataContainerItemIterator<const_iterator>;

  /// @}

  /**
   * @brief Default constructor: container has no room at all.
   * @see `resize()`
   *
   * The object *must* be resized before being of any use.
   */
  GeoIDdataContainer() = default;

  /**
   * @brief Prepares the container with default-constructed data.
   * @param dims number of elements on all levels of the container
   * @see `resize()`
   *
   * The size of each dimension is specified by the corresponding number,
   * starting from the size of the outer dimension (cryostat).
   *
   * The container is sized to host data for all the elements.
   * Each element in the container is default-constructed.
   */
  GeoIDdataContainer(std::initializer_list<unsigned int> dims);

  /**
   * @brief Prepares the container initializing all its data.
   * @param dims number of elements on all levels of the container
   * @param defValue the value copied to fill all entries in the container
   * @see `resize()`
   *
   * The size of each dimension is specified by the corresponding number,
   * starting from the size of the outer dimension (cryostat).
   *
   * The container is sized to host data for all the elements.
   * Each element in the container is constructed as copy of `defValue`.
   */
  GeoIDdataContainer(std::initializer_list<unsigned int> dims, value_type const& defValue);

  // --- BEGIN Container status query ----------------------------------------
  /// @name Container status query
  /// @{

  /// Returns the number of elements in the container.
  size_type size() const;

  /// Returns the number of elements the container has memory for.
  size_type capacity() const;

  /// Returns whether the container has no elements (`false` by assumptions).
  bool empty() const;

  /// Dimensions of the `Level` dimension of this container.
  template <std::size_t Level>
  unsigned int dimSize() const;

  /// Dimensions of the ID of this container.
  static constexpr unsigned int dimensions();

  /// Returns whether this container hosts data for the specified ID.
  template <typename GeoID>
  bool hasElement(GeoID const& id) const;

  /// Returns the ID of the first element with GeoID type.
  template <typename GeoID = ID_t>
  GeoID firstID() const;

  /// Returns the ID of the last covered element with GeoID type.
  template <typename GeoID = ID_t>
  GeoID lastID() const;

  /// Returns the mapper object used to convert ID's and container positions.
  Mapper_t const& mapper() const;

  /// @}
  // --- END Container status query ------------------------------------------

  // --- BEGIN Element access ------------------------------------------------
  /// @name Element access
  /// @{

  /// Returns the element for the specified geometry element.

  reference operator[](ID_t const& id);

  /// Returns the element for the specified geometry element (read-only).
  const_reference operator[](ID_t const& id) const;

  /// Returns the element for the specified geometry element.
  /// @throw std::out_of_range if element `id` is not within the container range
  reference at(ID_t const& id);

  /// Returns the element for the specified geometry element (read-only).
  /// @throw std::out_of_range if element `id` is not within the container range
  const_reference at(ID_t const& id) const;

  /// Returns the element for the first ID (unchecked).
  reference first();

  /// Returns the element for the first ID (unchecked).
  const_reference first() const;

  /// Returns the element for the last ID (unchecked).
  reference last();

  /// Returns the element for the last ID (unchecked).
  const_reference last() const;

  /// @}
  // --- END Element access --------------------------------------------------

  // --- BEGIN Iterators -----------------------------------------------------
  /**
   * @name Iterators
   *
   * Two types of iterators are provided:
   *
   * 1. "standard" iterators pointing to data values
   * 2. "item" pseudo-iterators dereferencing to a (ID, value) pair
   *
   * Reverse iterators are not supported (yet?).
   *
   * Standard iterators
   * -------------------
   *
   * The STL-like interface provides iterators that go through the entire range
   * of allowed data, i.e. all the `size()` elements that are also reached
   * via random access (`operator[]()`).
   *
   * These iterators have an interface extension: the data member `ID()` returns
   * the ID of the element the iterator is pointing to. For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto iData = data.begin();
   * auto const dend = data.end();
   * while (iData != dend) {
   *   std::cout << "data[" << iData.ID() << "] = " << *iData << std::endl;
   *   ++iData;
   * } // while
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Note that using the range-for loop, you don't get access to the iterator
   * and therefore not even to the ID.
   *
   *
   * Item iterators
   * ---------------
   *
   * The item iterators are iterators adapted from the standard ones, which
   * when dereferenced return a pair ( `ID_t`, `reference` ).
   * They can be accessed with `item_begin()`, `item_end()` etc, and range-for
   * loop can be obtained via `items()` member function:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * for (auto&& [ ID, value ]: data.items()) {
   *   std::cout << "data[" << ID << "] = " << value << std::endl;
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (this loop has the same effect as the one in the example of the standard
   * iterators, but it's more compact).
   *
   */
  /// @{

  /// Returns an iterator to the beginning of the data.
  iterator begin();

  /// Returns an iterator to past the end of the data.
  iterator end();

  /// Returns a constant iterator to the beginning of the data.
  const_iterator begin() const;

  /// Returns a constant iterator to past the end of the data.
  const_iterator end() const;

  /// Returns a constant iterator to the beginning of the data.
  const_iterator cbegin() const;

  /// Returns a constant iterator to past the end of the data.
  const_iterator cend() const;

  /// Returns an item iterator to the beginning of the data.
  item_iterator item_begin();

  /// Returns an item iterator to past the end of the data.
  item_iterator item_end();

  /// Returns a item constant iterator to the beginning of the data.
  item_const_iterator item_begin() const;

  /// Returns a item constant iterator to past the end of the data.
  item_const_iterator item_end() const;

  /// Returns a item constant iterator to the beginning of the data.
  item_const_iterator item_cbegin() const;

  /// Returns a item constant iterator to past the end of the data.
  item_const_iterator item_cend() const;

  /// Returns an object suitable for a range-for loop with `item_iterator`.
  auto items();

  /// Returns an object suitable for a range-for loop with `item_const_iterator`.
  auto items() const;

  /// @}
  // --- END Iterators -------------------------------------------------------

  // --- BEGIN Data modification ---------------------------------------------
  /**
   * @name Data modification
   *
   * In general, each single element can be accessed and changed.
   * In addition, this section includes methods acting on multiple elements
   * at once.
   */
  /// @{

  /// Sets all elements to the specified `value` (copied).
  void fill(value_type value);

  /// Sets all the elements to a default-constructed `value_type`.
  void reset();

  /**
   * @brief Applies an operation on all elements.
   * @tparam Op type of operation
   * @param op Operation
   * @return the operation object after operations took place
   *
   * The operation `op` is a unary functor, i.e. an object that supports the
   * call to `op(value_type&)`.
   *
   * The return values of `op` calls are discarded.
   */
  template <typename Op>
  Op apply(Op&& op);

  /**
   * @brief Applies an operation on all elements.
   * @tparam Op type of operation
   * @param op Operation
   * @return the operation object after operations took place
   *
   * The operation `op` is a unary functor, i.e. an object that supports the
   * call to `op(value_type const&)`.
   *
   * The return values of `op` calls are discarded.
   * Note that while the elements of this container can't be modified, the
   * operation itself still can if not constant, and it is returned.
   */
  template <typename Op>
  decltype(auto) apply(Op&& op) const;

  /// @}
  // --- END Data modification -------------------------------------------------

  // --- BEGIN Container modification ------------------------------------------
  /// @name Container modification
  /// @{

  /**
   * @brief Prepares the container with default-constructed data.
   * @param dims number of elements on all levels of the container
   * @see `clear()`, `fill()`
   *
   * The size of each dimension is specified by the corresponding number,
   * starting from the size of the outer dimension (cryostat).
   *
   * The container is sized to host data for all the elements.
   * Each new element in the container is default-constructed.
   * Existing data is not touched, but it may be rearranged in a
   * non-straightforward way.
   */
  void resize(std::initializer_list<unsigned int> dims);

  /**
   * @brief Prepares the container initializing all its data.
   * @param dims number of elements on all levels of the container
   * @param defValue the value copied to fill all entries in the container
   * @see `clear()`, `fill()`
   *
   * The size of each dimension is specified by the corresponding number,
   * starting from the size of the outer dimension (cryostat).
   *
   * The container is sized to host data for all the elements.
   * Each new element in the container is constructed as copy of `defValue`.
   * Existing data is not touched, but it may be rearranged in a
   * non-straightforward way.
   */
  void resize(std::initializer_list<unsigned int> dims, value_type const& defValue);

  /**
   * @brief Prepares the container with default-constructed data.
   * @param other data collection to take dimensions from
   *
   * The size of each dimension is taken by the matching one in `other`.
   *
   * The container is sized to host data for all the elements.
   * Each new element in the container is default-constructed.
   * Existing data is not touched, but it may be rearranged in a
   * non-straightforward way.
   */
  template <typename OT>
  void resizeAs(GeoIDdataContainer<OT, Mapper_t> const& other);

  /**
   * @brief Prepares the container initializing all its data.
   * @param other data collection to take dimensions from
   * @param defValue the value copied to fill all entries in the container
   *
   * The size of each dimension is taken by the matching one in `other`.
   *
   * The container is sized to host data for all the elements.
   * Each new element in the container is constructed as copy of `defValue`.
   * Existing data is not touched, but it may be rearranged in a
   * non-straightforward way.
   */
  template <typename OT>
  void resizeAs(GeoIDdataContainer<OT, Mapper_t> const& other, value_type const& defValue);

  /**
   * @brief Makes the container empty, with no usable storage space.
   * @see `resize()`
   *
   * The container needs to be resized before it is useful again.
   */
  void clear();

  /// @}
  // --- END Container modification --------------------------------------------

private:
  Mapper_t fMapper; ///< Mapping of IDs to indices.

  Container_t fData; ///< Data storage.

  /// Returns the internal index of the specified ID in the storage area.
  size_type index(ID_t const& id) const;

  /// Returns the ID corresponding to a internal index in the storage area.
  ID_t ID(size_type const index) const;

}; // class geo::GeoIDdataContainer<>

//------------------------------------------------------------------------------
/**
 * @brief Container with one element per geometry TPC.
 * @tparam T type of the contained datum
 * @see `geo::GeometryCore::makeTPCData`
 *
 * The container is of fixed size and can't be neither resized nor freed
 * before destruction.
 *
 * This example creates a "map" of tracks starting on each TPC:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto const* geom = lar::providerFrom<geo::GeometryCore>();
 * geo::TPCDataContainer<std::vector<recob::Track const*>> TracksPerTPC
 *   (geom->NCryostats(), geom->MaxTPCs());
 *
 * for (recob::Track const& track: tracks) {
 *   geo::TPCGeo const* tpc = geom->PositionToTPCptr(track.Start());
 *   if (tpc) TracksPerTPC[tpc->ID()].push_back(tpc);
 * } // for
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 * Assumptions
 * ============
 *
 * The following assumptions should be considered unchecked, and the behavior
 * when they are violated undefined (but note that in debug mode some of them
 * might be actually checked):
 * * the container assumes the same number of TPCs in each cryostat. While
 *   this is not effectively a necessary condition, keep in mind that this
 *   container has no notion whether a given TPC actually exists in the
 *   geometry or not
 * * at least one element is expected to be present
 *
 */
template <typename T>
class geo::TPCDataContainer : public geo::GeoIDdataContainer<T, geo::TPCIDmapper<>> {

  using BaseContainer_t = GeoIDdataContainer<T, TPCIDmapper<>>;

public:
  using value_type = typename BaseContainer_t::value_type;

  /**
   * @brief Default constructor: empty container.
   * @see `resize()`
   *
   * The container starts with no room for any data.
   * The only guarantee is that `empty()` is `true` and `size()` is `0`.
   * Use `resize()` before anything else.
   */
  TPCDataContainer() = default;

  /**
   * @brief Prepares the container with default-constructed data.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs. Each element in the container is default-constructed.
   */
  TPCDataContainer(unsigned int nCryo, unsigned int nTPCs) : BaseContainer_t({nCryo, nTPCs}) {}

  /**
   * @brief Prepares the container with copies of the specified default value.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @param defValue the value to be replicated
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs. Each element in the container is a copy of defValue.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const* geom = lar::providerFrom<geo::GeometryCore>();
   * geo::TPCDataContainer<unsigned int> PlanesPerTPC
   *   (geom->NCryostats(), geom->MaxTPCs(), 3U);
   * for (geo::TPCGeo const& TPC: geom->IterateTPC())
   *   assert(PlanesPerTPC[TPC.ID()] == TPC.Nplanes());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  TPCDataContainer(unsigned int nCryo, unsigned int nTPCs, value_type const& defValue)
    : BaseContainer_t({nCryo, nTPCs}, defValue)
  {}

  // --- BEGIN Container modification ------------------------------------------
  /// @name Container modification
  /// @{

  /**
   * @brief Prepares the container with default-constructed data.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @see `clear()`, `fill()`
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs. Each element in the container is default-constructed.
   *
   * Existing data is not touched, but it may be rearranged in a
   * non-straightforward way.
   */
  void resize(unsigned int nCryo, unsigned int nTPCs) { BaseContainer_t::resize({nCryo, nTPCs}); }

  /**
   * @brief Prepares the container initializing all its data.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @param defValue the value copied to fill all entries in the container
   * @see `clear()`, `fill()`
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs. Each element in the container is a copy of `defValue`.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const* geom = lar::providerFrom<geo::GeometryCore>();
   * geo::TPCDataContainer<unsigned int> countPerPlane
   * countPerPlane.resize(geom->NCryostats(), geom->MaxTPCs(), 0U);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * Existing data is not touched, but it may be rearranged in a
   * non-straightforward way.
   */
  void resize(unsigned int nCryo, unsigned int nTPCs, T const& defValue)
  {
    BaseContainer_t::resize({nCryo, nTPCs}, defValue);
  }

  /// @}
  // --- END Container modification --------------------------------------------

  // --- BEGIN Container status query ------------------------------------------
  /// @name Container status query
  /// @{

  /// Returns whether this container hosts data for the specified cryostat.
  bool hasCryostat(CryostatID const& cryoid) const { return BaseContainer_t::hasElement(cryoid); }

  /// Returns whether this container hosts data for the specified TPC.
  bool hasTPC(TPCID const& tpcid) const { return BaseContainer_t::hasElement(tpcid); }

  /// @}
  // --- END Container status query --------------------------------------------

}; // class geo::details::TPCDataContainer<>

//------------------------------------------------------------------------------
/**
 * @brief Container with one element per geometry wire plane.
 * @tparam T type of the contained datum
 * @see `geo::GeometryCore::makePlaneData`
 *
 * The container is of fixed size and can't be neither resized nor freed
 * before destruction.
 *
 * This example creates a "map" of hits on each wire plane:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto const* geom = lar::providerFrom<geo::GeometryCore>();
 * geo::PlaneDataContainer<std::vector<recob::Hit const*>> hitsPerPlane
 *   (geom->NCryostats(), geom->MaxTPCs(), geom->MaxPlanes());
 *
 * for (recob::Hit const& hit: hits) {
 *   if (hit.WireID()) hitsPerPlane[hit.WireID().planeID()].push_back(&hit);
 * } // for
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 * Assumptions
 * ============
 *
 * The following assumptions should be considered unchecked, and the behavior
 * when they are violated undefined (but note that in debug mode some of them
 * might be actually checked):
 * * the container assumes the same number of planes in all TPCs, and of TPCs
 *   on all cryostats. While this is not effectively a necessary condition,
 *   keep in mind that this container has no notion whether a given plane or
 *   TPC actually exists in the geometry or not
 * * at least one element is expected to be present
 *
 */
template <typename T>
class geo::PlaneDataContainer : public geo::GeoIDdataContainer<T, geo::PlaneIDmapper<>> {

  /// Base class.
  using BaseContainer_t = GeoIDdataContainer<T, PlaneIDmapper<>>;

public:
  /**
   * @brief Default constructor: empty container.
   * @see `resize()`
   *
   * The container starts with no room for any data.
   * The only guarantee is that `empty()` is `true` and `size()` is `0`.
   * Use `resize()` before anything else.
   */
  PlaneDataContainer() = default;

  /**
   * @brief Prepares the container with default-constructed data.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs per cryostat
   * @param nPlanes number of planes per TPC
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs, each one with `nPlanes` wire planes. Each element in the
   * container is default-constructed.
   */
  PlaneDataContainer(unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes)
    : BaseContainer_t({nCryo, nTPCs, nPlanes})
  {}

  /**
   * @brief Prepares the container with copies of the specified default value.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @param nPlanes number of planes per TPC
   * @param defValue the value to be replicated
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs, and each of them with `nPlanes` planes. Each element in
   * the container is a copy of `defValue`.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const* geom = lar::providerFrom<geo::GeometryCore>();
   * geo::PlaneDataContainer<unsigned int> countPerPlane
   *   (geom->NCryostats(), geom->MaxTPCs(), geom->MaxPlanes(), 0U);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  PlaneDataContainer(unsigned int nCryo,
                     unsigned int nTPCs,
                     unsigned int nPlanes,
                     T const& defValue)
    : BaseContainer_t{{nCryo, nTPCs, nPlanes}, defValue}
  {}

  // --- BEGIN Container modification ------------------------------------------
  /// @name Container modification
  /// @{

  /**
   * @brief Prepares the container with default-constructed data.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @param nPlanes number of planes per TPC
   * @see `clear()`, `fill()`
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs, each one with `nPlanes` wire planes. Each element in the
   * container is default-constructed.
   *
   * Existing data is not touched, but it may be rearranged in a
   * non-straightforward way.
   */
  void resize(unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes)
  {
    BaseContainer_t::resize({nCryo, nTPCs, nPlanes});
  }

  /**
   * @brief Prepares the container initializing all its data.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @param nPlanes number of planes per TPC
   * @param defValue the value copied to fill all entries in the container
   * @see `clear()`, `fill()`
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs, and each of them with `nPlanes` planes. Each element in
   * the container is a copy of `defValue`.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const* geom = lar::providerFrom<geo::GeometryCore>();
   * geo::PlaneDataContainer<unsigned int> countPerPlane
   * countPerPlane.resize
   *   (geom->NCryostats(), geom->MaxTPCs(), geom->MaxPlanes(), 0U);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * Existing data is not touched, but it may be rearranged in a
   * non-straightforward way.
   */
  void resize(unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes, T const& defValue)
  {
    BaseContainer_t::resize({nCryo, nTPCs, nPlanes}, defValue);
  }

  /// @}
  // --- END Container modification --------------------------------------------

  // --- BEGIN Container status query ------------------------------------------
  /// @name Container status query
  /// @{

  /// Returns whether this container hosts data for the specified cryostat.
  bool hasCryostat(CryostatID const& cryoid) const { return BaseContainer_t::hasElement(cryoid); }

  /// Returns whether this container hosts data for the specified TPC.
  bool hasTPC(TPCID const& tpcid) const { return BaseContainer_t::hasElement(tpcid); }

  /// Returns whether this container hosts data for the specified plane.
  bool hasPlane(PlaneID const& planeid) const { return BaseContainer_t::hasElement(planeid); }

  /// @}
  // --- END Container status query --------------------------------------------

}; // class geo::PlaneDataContainer

//------------------------------------------------------------------------------
/**
 * @brief Iterator for `geo::GeoIDdataContainer` class.
 * @tparam GeoIDdataContainerClass type of the class being iterated
 * @tparam BaseIterator type of iterator to the actual data
 *
 * @note These iterators haven't been extensively tested. Caveat emptor...
 */
template <typename GeoIDmapperClass, typename BaseIterator>
class geo::details::GeoIDdataContainerIterator
  : public boost::iterator_adaptor<
      geo::details::GeoIDdataContainerIterator<GeoIDmapperClass, BaseIterator>,
      BaseIterator> {

  ///< Type of mapping of the container this class iterates.
  using Mapper_t = GeoIDmapperClass;

  using BaseIterator_t = BaseIterator; ///< Type of iterator to the actual data.

  /// Type of index in the container mapping.
  using Index_t = typename Mapper_t::index_type;

  struct ExtraData_t {

    /// Mapping of the container being iterated.
    Mapper_t const* mapper = nullptr;

    BaseIterator_t start; ///< Iterator to the first element.

  }; // struct ExtraData_t

  ExtraData_t fData; ///< Data for extended features of this iterator.

  /// Returns the iterator to the current element.
  BaseIterator_t const& current() const
  {
    return GeoIDdataContainerIterator::iterator_adaptor_::base();
  }
  BaseIterator_t& current() { return GeoIDdataContainerIterator::iterator_adaptor_::base(); }

  /// Returns the iterator to the begin element.
  BaseIterator_t const& start() const { return fData.start; }

  /// Returns the mapping of the container being iterated.
  Mapper_t const& mapper() const
  {
    assert(fData.mapper);
    return *(fData.mapper);
  }

  /// Returns the index of the current element
  Index_t index() const { return static_cast<Index_t>(current() - start()); }

public:
  using ID_t = typename Mapper_t::ID_t; ///< Type of the ID in this iterator.

  /// Default constructor: undefined status.
  GeoIDdataContainerIterator() = default;

  /// Constructor: points to data pointed by `current`.
  GeoIDdataContainerIterator(Mapper_t const& mapper,
                             BaseIterator_t const& start,
                             BaseIterator_t const& current)
    : GeoIDdataContainerIterator::iterator_adaptor_(current), fData{&mapper, start}
  {}

  /// Generalized copy constructor, only if argument iterator can be converted.
  template <typename OBaseIterator>
  GeoIDdataContainerIterator(
    GeoIDdataContainerIterator<Mapper_t, OBaseIterator> const& other,
    std::enable_if_t<std::is_convertible_v<OBaseIterator, BaseIterator_t>> = nullptr)
    : GeoIDdataContainerIterator::iterator_adaptor_(other.base()), fData(other.fData)
  {}

  /// Returns the ID corresponding to the current element.
  ID_t ID() const { return mapper().ID(index()); }

private:
  //   friend class boost::iterator_core_access;

  // there might be some need for other interoperability methods (comparison?)

}; // geo::details::GeoIDdataContainerIterator

/**
 * @brief Item iterator for `geo::GeoIDdataContainer` class.
 * @tparam GeoIDIteratorClass type of iterator being wrapped
 *
 * This iterator is just a wrapper.
 *
 * @note These iterators haven't been extensively tested. Caveat emptor...
 */
template <typename GeoIDIteratorClass>
class geo::details::GeoIDdataContainerItemIterator
  : public boost::iterator_adaptor<geo::details::GeoIDdataContainerItemIterator<GeoIDIteratorClass>,
                                   GeoIDIteratorClass,
                                   std::pair<typename GeoIDIteratorClass::ID_t,
                                             typename GeoIDIteratorClass::reference> // Value
                                   ,
                                   boost::use_default // Category
                                   ,
                                   std::pair<typename GeoIDIteratorClass::ID_t,
                                             typename GeoIDIteratorClass::reference> // Reference
                                   > {
  /*
   * Implementation notes:
   *  * boost::transform_iterator can't be used to wrap
   *      `geo::GeoIDdataContainerIterator` because it acts on the value of
   *      the wrapped iterator, while we need to extract information from
   *      the iterator itself (its `ID()` data member)
   *  * we need to specify the type of reference the iterator returns,
   *      because... it's not a reference; this is not really a STL random
   *      access iterator since it dereferences to a temporary value (rvalue)
   *      like an input iterator can; but for the rest it's a full blown random
   *      access iterator (same stuff as the infamous `std::vector<bool>`)
   *
   */

  using GeoIDiterator_t = GeoIDIteratorClass; ///< Type of wrapped iterator.

  using iterator_adaptor_ = typename GeoIDdataContainerItemIterator::iterator_adaptor_;

public:
  /// Type of the ID in this iterator.
  using ID_t = typename GeoIDiterator_t::ID_t;

  /// Default constructor: undefined status.
  GeoIDdataContainerItemIterator() = default;

  /// Constructor: points to data pointed by `current`.
  GeoIDdataContainerItemIterator(GeoIDiterator_t const& iter) : iterator_adaptor_(iter) {}

  /// Generalized copy constructor, only if argument iterator can be converted.
  template <typename OGeoIDIteratorClass>
  GeoIDdataContainerItemIterator(
    GeoIDdataContainerItemIterator<OGeoIDIteratorClass> const& other,
    std::enable_if_t<std::is_convertible_v<OGeoIDIteratorClass, GeoIDiterator_t>> = nullptr)
    : iterator_adaptor_(other.base())
  {}

private:
  friend class boost::iterator_core_access;

  using iterator_adaptor_::base;

  typename iterator_adaptor_::reference dereference() const { return {base().ID(), *base()}; }

  // there might be some need for other interoperability methods (comparison?)

}; // geo::details::GeoIDdataContainerItemIterator

/// @}
// --- END Geometry data containers --------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--- Template implementation
//------------------------------------------------------------------------------
//--- geo::details::GeoContainerData
//------------------------------------------------------------------------------
template <typename T>
class geo::details::GeoContainerData {

  using Container_t = std::vector<T>;

public:
  // --- BEGIN STL container types ---------------------------------------------
  /// @{
  /// @name STL container types.

  using value_type = typename Container_t::value_type;
  using reference = typename Container_t::reference;
  using const_reference = typename Container_t::const_reference;
  using pointer = typename Container_t::pointer;
  using const_pointer = typename Container_t::const_pointer;
  using iterator = typename Container_t::iterator;
  using const_iterator = typename Container_t::const_iterator;
  using reverse_iterator = typename Container_t::reverse_iterator;
  using const_reverse_iterator = typename Container_t::const_reverse_iterator;
  using difference_type = typename Container_t::difference_type;
  using size_type = typename Container_t::size_type;

  /// @}
  // --- END STL container types -----------------------------------------------

  using index_type = size_type; ///< Type used internally (so far) for indexing.

  // --- BEGIN Constructors ----------------------------------------------------

  /// Default constructor with empty container. Good for nothing.
  GeoContainerData() = default;

  /// Prepares the container with default-constructed data.
  GeoContainerData(size_type size) : fData(size) {}

  // Prepares the container with copies of the specified default value.
  GeoContainerData(size_type size, value_type const& defValue) : fData(size, defValue) {}

  // --- END Constructors ------------------------------------------------------

  // --- BEGIN Container status query ------------------------------------------
  /// @name Container status query
  /// @{

  /// Returns the number of elements in the container.
  size_type size() const { return fData.size(); }

  /// Returns the number of elements the container has memory for.
  size_type capacity() const { return fData.capacity(); }

  /// Returns whether the container has no elements (`false` by assumptions).
  bool empty() const { return fData.empty(); }

  /// @}
  // --- END Container status query --------------------------------------------

  // --- BEGIN Data modification -----------------------------------------------
  /**
   * @name Data modification
   *
   * In general, each single element can be accessed and changed.
   * In addition, this section includes methods acting on multiple elements
   * at once.
   */
  /// @{

  /// Sets all elements to the specified `value` (copied).
  void fill(value_type value) { std::fill(fData.begin(), fData.end(), value); }

  /// Sets all the elements to a default-constructed `value_type`.
  void reset() { fill(value_type{}); }

  /**
   * @brief Applies an operation on all elements.
   * @tparam Op type of operation
   * @param op Operation
   * @return the operation object after operations took place
   *
   * The operation `op` is a unary functor, i.e. an object that supports the
   * call to `op(value_type&)`.
   *
   * The return values of `op` calls are discarded.
   */
  template <typename Op>
  Op apply(Op&& op)
  {
    for (auto& data : fData)
      op(data);
    return op;
  }

  /**
   * @brief Applies an operation on all elements.
   * @tparam Op type of operation
   * @param op Operation
   * @return the operation object after operations took place
   *
   * The operation `op` is a unary functor, i.e. an object that supports the
   * call to `op(value_type const&)`.
   *
   * The return values of `op` calls are discarded.
   */
  template <typename Op>
  Op apply(Op&& op) const
  {
    for (auto const& data : fData)
      op(data);
    return op;
  }

  /// @}
  // --- END Element access ----------------------------------------------------

  // --- BEGIN Container modification ------------------------------------------
  /// @name Container modification
  /// @{

  /**
   * @brief Prepares the container with default-constructed data.
   * @param size number of elements in the container
   * @see `clear()`, `fill()`
   *
   * The container is sized to host data for all the elements.
   * Each new element in the container is default-constructed.
   * Existing data is not touched.
   */
  void resize(size_type size) { fData.resize(size); }

  /**
   * @brief Prepares the container with copies of the specified default value.
   * @param size number of elements in the container
   * @param defValue the value copied to fill all entries in the container
   * @see `clear()`, `fill()`
   *
   * The container is sized to host data for all the elements.
   * Each new element in the container is constructed as copy of `defValue`.
   * Existing data is not touched.
   */
  void resize(size_type size, value_type const& defValue) { fData.resize(size, defValue); }

  /**
   * @brief Makes the container empty, with no usable storage space.
   * @see `resize()`
   *
   * The container needs to be resized before it is useful again.
   */
  void clear() { fData.clear(); }

  /// @}
  // --- END Container modification --------------------------------------------

  // --- BEGIN Element access --------------------------------------------------
  /// @name Element access
  /// @{

  /// Returns the element for the specified index.
  reference operator[](index_type index) { return fData[index]; }

  /// Returns the element for the specified index (read-only).
  const_reference operator[](index_type index) const { return fData[index]; }

  /// @}
  // --- END Element access ----------------------------------------------------

  // --- BEGIN Iterators -------------------------------------------------------
  /// @name Iterators
  /// @{

  iterator begin() { return fData.begin(); }
  iterator end() { return fData.end(); }
  const_iterator begin() const { return fData.begin(); }
  const_iterator end() const { return fData.end(); }
  const_iterator cbegin() const { return fData.cbegin(); }
  const_iterator cend() const { return fData.cend(); }

  reverse_iterator rbegin() { return fData.rbegin(); }
  reverse_iterator rend() { return fData.rend(); }
  const_reverse_iterator rbegin() const { return fData.rbegin(); }
  const_reverse_iterator rend() const { return fData.rend(); }
  const_reverse_iterator crbegin() const { return fData.crbegin(); }
  const_reverse_iterator crend() const { return fData.crend(); }

  /// @}
  // --- END Iterators ---------------------------------------------------------

  /// Returns whether the specified value is between `0` and the upper limit.
  template <typename Value, typename Upper>
  static bool bounded(Value v, Upper upper)
  {
    return (v >= 0) && (static_cast<size_type>(v) < upper);
  }

private:
  Container_t fData; ///< Data storage area.

}; // class geo::details::GeoContainerData

//------------------------------------------------------------------------------
//--- geo::GeoIDdataContainer
//------------------------------------------------------------------------------
template <typename T, typename Mapper>
geo::GeoIDdataContainer<T, Mapper>::GeoIDdataContainer(std::initializer_list<unsigned int> dims)
  : fMapper(dims), fData(fMapper.size())
{
  assert(!fData.empty());
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
geo::GeoIDdataContainer<T, Mapper>::GeoIDdataContainer(std::initializer_list<unsigned int> dims,
                                                       value_type const& defValue)
  : fMapper(dims), fData(fMapper.computeSize(), defValue)
{
  assert(!fData.empty());
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::size() const -> size_type
{
  return fData.size();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::capacity() const -> size_type
{
  return fData.capacity();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
bool geo::GeoIDdataContainer<T, Mapper>::empty() const
{
  return fData.empty();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
template <std::size_t Level>
unsigned int geo::GeoIDdataContainer<T, Mapper>::dimSize() const
{
  return mapper().template dimSize<Level>();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
constexpr unsigned int geo::GeoIDdataContainer<T, Mapper>::dimensions()
{
  return Mapper_t::dimensions();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
template <typename GeoID>
bool geo::GeoIDdataContainer<T, Mapper>::hasElement(GeoID const& id) const
{
  return mapper().template hasElement<GeoID>(id);
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
template <typename GeoID /* = ID_t */>
GeoID geo::GeoIDdataContainer<T, Mapper>::firstID() const
{
  return mapper().template firstID<GeoID>();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
template <typename GeoID /* = ID_t */>
GeoID geo::GeoIDdataContainer<T, Mapper>::lastID() const
{
  return mapper().template lastID<GeoID>();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::mapper() const -> Mapper_t const&
{
  return fMapper;
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::operator[](ID_t const& id) -> reference
{
  return fData[mapper().index(id)];
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::operator[](ID_t const& id) const -> const_reference
{
  return fData[mapper().index(id)];
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::at(ID_t const& id) -> reference
{
  if (hasElement(id)) return operator[](id);
  throw std::out_of_range("No data for " + std::string(id));
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::at(ID_t const& id) const -> const_reference
{
  if (hasElement(id)) return operator[](id);
  throw std::out_of_range("No data for " + std::string(id));
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::first() -> reference
{
  return operator[](firstID());
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::first() const -> const_reference
{
  return operator[](firstID());
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::last() -> reference
{
  return operator[](lastID());
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::last() const -> const_reference
{
  return operator[](lastID());
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::begin() -> iterator
{
  return {mapper(), fData.begin(), fData.begin()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::end() -> iterator
{
  return {mapper(), fData.begin(), fData.end()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::begin() const -> const_iterator
{
  return {mapper(), fData.begin(), fData.begin()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::end() const -> const_iterator
{
  return {mapper(), fData.begin(), fData.end()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::cbegin() const -> const_iterator
{
  return begin();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::cend() const -> const_iterator
{
  return end();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::item_begin() -> item_iterator
{
  return {begin()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::item_end() -> item_iterator
{
  return {end()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::item_begin() const -> item_const_iterator
{
  return {begin()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::item_end() const -> item_const_iterator
{
  return {end()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::item_cbegin() const -> item_const_iterator
{
  return item_begin();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::item_cend() const -> item_const_iterator
{
  return item_end();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::items()
{
  return util::span{item_begin(), item_end()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::items() const
{
  return util::span{item_begin(), item_end()};
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
void geo::GeoIDdataContainer<T, Mapper>::fill(value_type value)
{
  fData.fill(value);
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
void geo::GeoIDdataContainer<T, Mapper>::reset()
{
  fData.reset();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
template <typename Op>
Op geo::GeoIDdataContainer<T, Mapper>::apply(Op&& op)
{
  return fData.apply(std::forward<Op>(op));
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
void geo::GeoIDdataContainer<T, Mapper>::resize(std::initializer_list<unsigned int> dims)
{
  fMapper.resize(dims);
  fData.resize(mapper().size());
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
void geo::GeoIDdataContainer<T, Mapper>::resize(std::initializer_list<unsigned int> dims,
                                                value_type const& defValue)
{
  fMapper.resize(dims);
  fData.resize(mapper().size(), defValue);
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
template <typename OT>
void geo::GeoIDdataContainer<T, Mapper>::resizeAs(GeoIDdataContainer<OT, Mapper_t> const& other)
{
  fMapper.resizeAs(other.mapper());
  fData.resize(mapper().size());
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
template <typename OT>
void geo::GeoIDdataContainer<T, Mapper>::resizeAs(GeoIDdataContainer<OT, Mapper_t> const& other,
                                                  value_type const& defValue)
{
  fMapper.resizeAs(other.mapper());
  fData.resize(mapper().size(), defValue);
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
void geo::GeoIDdataContainer<T, Mapper>::clear()
{
  fMapper.clear();
  fData.clear();
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
template <typename Op>
decltype(auto) geo::GeoIDdataContainer<T, Mapper>::apply(Op&& op) const
{
  return fData.apply(std::forward<Op>(op));
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::index(ID_t const& id) const -> size_type
{
  return mapper().index(id);
}

//------------------------------------------------------------------------------
template <typename T, typename Mapper>
auto geo::GeoIDdataContainer<T, Mapper>::ID(size_type const index) const -> ID_t
{
  return mapper().ID(index);
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_GEOMETRYDATACONTAINERS_H
