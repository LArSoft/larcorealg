#ifndef LARCOREALG_GEOMETRY_DETAILS_GEOMETRY_ITERATORS_H
#define LARCOREALG_GEOMETRY_DETAILS_GEOMETRY_ITERATORS_H

// LArSoft libraries
#include "larcorealg/Geometry/details/helpers.h"
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

// C/C++ standard libraries
#include <iterator> // std::forward_iterator_tag
#include <string>
#include <utility>

namespace geo::details {

  /// Base class for geometry iterators (note: this is not an iterator)
  class geometry_iterator_base {
  public:
    /// Constructor: associates with the specified geometry
    geometry_iterator_base(GeometryCore const* geom) : pGeo(geom) {}

  protected:
    /// Returns a pointer to the geometry
    GeometryCore const* geometry() const { return pGeo; }

    /// Default constructor; do not use a default-constructed iterator as-is!
    geometry_iterator_base() = default;

  private:
    GeometryCore const* pGeo = nullptr; ///< pointer to the geometry

  }; // class geometry_iterator_base

  /**
     * @brief Base forward iterator browsing all cryostat IDs in the detector
     * @tparam GEOID ID type to be used
     *
     * This iterator assumes that GEOID is derived from geo::CryostatID.
     * Note that no polymorphic behaviour is required, or expected, from GEOID.
     *
     * This iterator is designed to carry on, untouched, anything else that the
     * GEOID type defines beyond the required CryostatID data.
     *
     * Currently, backward iterations are not supported.
     */
  template <typename LocalID, typename GEOID>
  class id_iterator_base;

  template <typename GEOID>
  class id_iterator_base<CryostatID, GEOID> : public geometry_iterator_base {
  public:
    using GeoID_t = GEOID; ///< type of the actual ID stored in the iterator

    using LocalID_t = CryostatID;                          ///< type of the ID we change
    using iterator = id_iterator_base<LocalID_t, GeoID_t>; ///< this iterator

    static_assert(std::is_base_of<LocalID_t, GEOID>{}, "template type GEOID is not a LocalID_t");

    /// @name Iterator traits
    /// @{
    using difference_type = std::ptrdiff_t;
    using value_type = LocalID_t;
    using reference = value_type const&;
    using pointer = value_type const*;
    using iterator_category = std::forward_iterator_tag;
    /// @}

    /// Default constructor; effect not defined: assign to it before using!
    id_iterator_base() = default;

    /// Constructor: points to the specified cryostat
    id_iterator_base(GeometryCore const* geom, GeoID_t const& start_from)
      : geometry_iterator_base{geom}, id{start_from}, limit{NSiblings(geom, localID())}
    {}

    /// Returns true if the two iterators point to the same cryostat
    bool operator==(id_iterator_base const& as) const { return localID() == as.localID(); }

    bool operator!=(id_iterator_base const& as) const { return localID() != as.localID(); }

    /// Returns the ID the iterator points to
    reference operator*() const { return localID(); }

    /// Returns a pointer to the ID the iterator points to
    pointer operator->() const { return &(localID()); }

    /// Prefix increment: returns this iterator pointing to the next cryostat
    iterator& operator++()
    {
      next();
      return *this;
    }

    /// Postfix increment: returns the current iterator, then increments it
    iterator operator++(int)
    {
      iterator old(*this);
      next();
      return old;
    }

  protected:
    using ID_t = typename LocalID_t::CryostatID_t;

    //@{
    /// Returns the actual type of ID we store
    GeoID_t const& ID() const { return id; }
    GeoID_t& ID() { return id; }
    //@}

    /// Skips to the next cryostat
    void next()
    {
      if (at_end()) return;
      if (++local_index() < limit) return;
      localID().isValid = false;
    }

    /// Returns whether this iterator has reached the end
    bool at_end() const { return local_index() == limit; }

  private:
    GeoID_t id{};                      ///< ID of the current cryostat
    ID_t limit = LocalID_t::InvalidID; ///< maximum number of cryostats

    //@{
    /// Returns the type of ID we act on
    LocalID_t const& localID() const { return ID(); }
    LocalID_t& localID() { return ID(); }
    //@}

    //@{
    /// Returns the index (part if the ID) this iterator runs on
    ID_t const& local_index() const { return localID().deepestIndex(); }
    ID_t& local_index() { return localID().deepestIndex(); }
    //@}
  }; // class id_iterator_base<CryostatID>

  /**
     * @brief Base forward iterator browsing all TPC IDs in the detector
     * @tparam GEOID ID type to be used
     *
     * This iterator requires that GEOID is derived from geo::TPCID.
     * Note that no polymorphic behaviour is required, or expected, from GEOID.
     *
     * This iterator is designed to carry on, untouched, anything else that the
     * GEOID type defines beyond the required TPCID data.
     *
     * @note A number of "local" methods are overloaded: since there is no
     * polymorphism here and they are not virtual functions, these are designed
     * not to replace the inherited methods except within the non-inherited and
     * explicitly redefined methods.
     *
     * Currently, backward iterations are not supported.
     */
  template <typename LocalID, typename GEOID>
  class id_iterator_base : protected id_iterator_base<typename LocalID::ParentID_t, GEOID> {
    using upper_iterator = id_iterator_base<typename LocalID::ParentID_t, GEOID>;

  public:
    using GeoID_t = typename upper_iterator::GeoID_t;

    using LocalID_t = LocalID; ///< type of the ID we change
    static_assert(std::is_base_of<LocalID_t, GEOID>{}, "template type GEOID is not a LocalID_t");

    using iterator = id_iterator_base<LocalID_t, GeoID_t>; ///< type of this iterator

    /// @name Iterator traits
    /// @{
    using difference_type = std::ptrdiff_t;
    using value_type = LocalID_t;
    using reference = value_type const&;
    using pointer = value_type const*;
    using iterator_category = std::forward_iterator_tag;
    /// @}

    /// Default constructor; effect not defined: assign to it before using!
    id_iterator_base() = default;

    /// Constructor: points to the specified element
    id_iterator_base(GeometryCore const* geo, GeoID_t const& start_from)
      : upper_iterator{geo, start_from}, geom{geo}, limit(NSiblings(geom, localID()))
    {}

    /// Returns true if the two iterators point to the same element
    bool operator==(id_iterator_base const& as) const { return localID() == as.localID(); }

    /// Returns true if the two iterators point to different elements
    bool operator!=(id_iterator_base const& as) const { return !operator==(as); }

    /// Returns the element ID the iterator points to
    reference operator*() const { return localID(); }

    /// Returns the element ID the iterator points to
    pointer operator->() const { return &(localID()); }

    /// Prefix increment: returns this iterator pointing to the next element
    iterator& operator++()
    {
      next();
      return *this;
    }

    /// Postfix increment: returns the current iterator, then increments it
    iterator operator++(int)
    {
      iterator old(*this);
      next();
      return old;
    }

  protected:
    /// Returns the type of ID we act on
    LocalID_t const& localID() const { return upper_iterator::ID(); }

    using ID_t = std::remove_reference_t<
      decltype(std::declval<LocalID_t>().deepestIndex())>; ///< specific type for element ID

    /// Skips to the next element
    void next()
    {
      // if at end (checked in the inherited context), do nothing
      if (upper_iterator::at_end()) return;

      // if after incrementing we haven't reached the limit, we are done
      if (++local_index() < limit) return;

      // we reached the end of the current elements list, we need to escalate:
      // - go to the next parent; if that becomes invalid, too bad, but we go on
      upper_iterator::next();
      // - set the index to the first element of the new parent
      local_index() = 0;
      // - update how many elements there are
      //   (expect 0 if it is now at_end() -- and it does not even matter)
      limit = NSiblings(geom, localID());
    }

    /// Returns the index (part if the ID) this iterator runs on
    ID_t const& local_index() const { return localID().deepestIndex(); }

  private:
    geo::GeometryCore const* geom;
    /// maximum number of elements in the current cryostat
    ID_t limit = LocalID_t::InvalidID;

    /// Returns the type of ID we act on (non-const version)
    LocalID_t& localID() { return upper_iterator::ID(); }

    /// Returns the index (part if the ID) this iterator runs on  (non-const)
    ID_t& local_index() { return localID().deepestIndex(); }
  }; // class id_iterator_base

  template <typename LocalID>
  using id_iterator = id_iterator_base<LocalID, LocalID>;

  /// Stream output for all geometry ID iterator types: prints the pointed ID.
  template <typename GEOIT>
  std::enable_if_t<std::is_base_of_v<geometry_iterator_base, GEOIT>, std::ostream&> operator<<(
    std::ostream& out,
    GEOIT const& it)
  {
    return out << "geometry_iterator{ " << *it << " }";
  }

  /**
     * @brief Forward iterator browsing all geometry elements in the detector
     * @tparam GEOITER type of geometry ID iterator
     *
     * This iterator works as the corresponding ID iterator in the template
     * argument. The difference is the dereferenciation operator: this one
     * obtains the geometry element directly, or throws on failure.
     * The boolean conversion operator checks that it can obtain a pointer to
     * the geometry element.
     *
     * In particular, get() and ID() methods still return the pointer to the
     * geometry element and its ID, respectively.
     *
     * It can also be initialized and compare with the corresponding ID
     * iterator.
     */
  template <typename Element, typename GEOIDITER>
  class geometry_element_iterator {
  public:
    using id_iterator_t = GEOIDITER;

    static_assert(std::is_base_of<geometry_iterator_base, id_iterator_t>{},
                  "template class for geometry_element_iterator"
                  " must be a geometry iterator");

    using iterator = geometry_element_iterator<Element, id_iterator_t>; ///< this type

    /// @{
    /// @name Types mirrored from the ID iterator
    using LocalID_t = typename id_iterator_t::LocalID_t;
    using GeoID_t = typename id_iterator_t::GeoID_t;
    using ElementPtr_t = Element const*;
    /// @}

    /// @name Iterator traits
    /// @{
    using difference_type = std::ptrdiff_t;
    using value_type = Element;
    using reference = value_type const&;
    using pointer = value_type const*;
    using iterator_category = std::forward_iterator_tag;
    /// @}

    /// Default constructor; effect not defined: assign to it before using!
    geometry_element_iterator() = default;

    /// Constructor: points to the same element as the specified ID iterator.
    geometry_element_iterator(GeometryCore const* geo, id_iterator_t const& iter)
      : geom{geo}, id_iterator(iter)
    {}

    /// Constructor: points to the same element as the specified ID iterator.
    geometry_element_iterator(GeometryCore const* geo, id_iterator_t&& iter)
      : geom{geo}, id_iterator(iter)
    {}

    /// Constructor: points to the specified geometry element
    geometry_element_iterator(GeometryCore const* geo, GeoID_t const& start_from)
      : geom{geo}, id_iterator(geom, start_from)
    {}

    /// Returns true if the two iterators point to the same object
    bool operator==(iterator const& as) const { return id_iterator == as.id_iterator; }

    /// Returns true if the two iterators point to different objects
    bool operator!=(iterator const& as) const { return id_iterator != as.id_iterator; }

    /**
       * @brief Returns the geometry element the iterator points to
       * @return a constant reference to the element the iterator points to
       * @throw cet::exception (category "geometry_iterator") if no valid
       *   geometry element is currently pointed by the iterator
       */
    reference operator*() const
    {
      ElementPtr_t ptr = get();
      if (ptr) return *ptr;
      throw cet::exception("geometry_iterator")
        << "iterator attempted to obtain geometry element " << std::string(ID());
    } // operator*()

    /// Returns a pointer to the element the iterator points to (or nullptr)
    pointer operator->() const { return get(); }

    /// Prefix increment: returns this iterator pointing to the next element
    iterator& operator++()
    {
      ++id_iterator;
      return *this;
    }

    /// Postfix increment: returns the current iterator, then increments it
    iterator operator++(int)
    {
      iterator old(*this);
      ++id_iterator;
      return old;
    }

    /// Returns whether the iterator is pointing to a valid geometry element
    operator bool() const { return validElement(geom, *id_iterator); }

    /// Returns a pointer to the geometry element, or nullptr if invalid
    ElementPtr_t get() const { return getElementPtr(geom, *id_iterator); }

    /// Returns the ID of the pointed geometry element
    LocalID_t const& ID() const { return *id_iterator; }

  private:
    geo::GeometryCore const* geom;
    id_iterator_t id_iterator;
  }; // class geometry_element_iterator<>

  // Element here supports types like CryostatGeo, etc.
  template <typename Element>
  using element_iterator_for =
    details::geometry_element_iterator<Element, details::id_iterator<typename Element::ID_t>>;

} // namespace geo::details

#endif // LARCOREALG_GEOMETRY_DETAILS_GEOMETRY_ITERATORS_H
