/**
 * @file   larcorealg/Geometry/GeometryCore.h
 * @brief  Access the description of detector geometry
 * @author brebel@fnal.gov
 * @see    larcorealg/Geometry/GeometryCore.cxx
 * @ingroup Geometry
 *
 * Structure of the header:
 *
 *   namespace geo {
 *
 *     // forward class declarations
 *
 *     namespace details {
 *
 *       // geometry iterator base class
 *
 *     }
 *
 *     // geometry iterators declaration
 *     //  - cryostat_id_iterator
 *     //  - TPC_id_iterator
 *     //  - plane_id_iterator
 *     //  - wire_id_iterator
 *     //  - TPCset_id_iterator
 *     //  - ROP_id_iterator
 *
 *     // GeometryData_t definition (part of GeometryCore)
 *
 *     // GeometryCore declaration
 *   }
 *
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYCORE_H
#define LARCOREALG_GEOMETRY_GEOMETRYCORE_H

// LArSoft libraries
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryData.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/ReadoutDataContainers.h" // readout::ROPDataContainer
#include "larcorealg/Geometry/GeometryDataContainers.h" // geo::TPCDataContainer
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect namespace
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/CoreUtils/span.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

// Framework and infrastructure libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT libraries
#include "TVector3.h"

// C/C++ standard libraries
#include <cstddef> // size_t
#include <string>
#include <vector>
#include <set>
#include <memory> // std::shared_ptr<>
#include <iterator> // std::forward_iterator_tag
#include <type_traits> // std::is_base_of<>

// ROOT class prototypes
class TGeoManager;
class TGeoNode;
class TGeoVolume;
class TGeoMaterial;

/// Namespace collecting geometry-related classes utilities
namespace geo {

  // Forward declarations within namespace.
  class AuxDetGeo;
  class AuxDetSensitiveGeo;
  class OpDetGeo;
  class GeometryCore;

  //
  // iterators
  //

  namespace details {

    /// Base class for geometry iterators, containing some type definitions
    class geometry_iterator_types {
        public:

      //@{
      /// Structures to distinguish the constructors.
      struct BeginPos_t {};
      struct EndPos_t {};
      struct UndefinedPos_t {};

      static constexpr BeginPos_t begin_pos = {};
      static constexpr EndPos_t end_pos = {};
      static constexpr UndefinedPos_t undefined_pos = {};
      //@}

    }; // class geometry_iterator_types

    /// Base class for geometry iterators (note: this is not an iterator)
    class geometry_iterator_base: public geometry_iterator_types {
        public:

      /// Constructor: associates with the specified geometry
      geometry_iterator_base(geo::GeometryCore const* geom): pGeo(geom) {}

        protected:
      /// Returns a pointer to the geometry
      geo::GeometryCore const* geometry() const { return pGeo; }

      /// Default constructor; do not use a default-constructed iterator as-is!
      geometry_iterator_base() {}

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
    template <typename GEOID>
    class cryostat_id_iterator_base:
      virtual public std::forward_iterator_tag, public geometry_iterator_base
    {
        public:
      using ElementPtr_t = geo::CryostatGeo const*;
      using GeoID_t = GEOID; ///< type of the actual ID stored in the iterator

      using iterator = cryostat_id_iterator_base<GeoID_t>; ///< this iterator

      using LocalID_t = geo::CryostatID; ///< type of the ID we change
      static_assert(std::is_base_of<LocalID_t, GEOID>::value,
        "template type GEOID is not a LocalID_t");


      /// @name Iterator traits
      /// @{
      using difference_type = std::ptrdiff_t;
      using value_type = LocalID_t;
      using reference = value_type const&;
      using pointer  = value_type const*;
      using iterator_category = std::input_iterator_tag;
      /// @}


      /// Default constructor; effect not defined: assign to it before using!
      cryostat_id_iterator_base() {}

      /// Constructor: points to begin
      cryostat_id_iterator_base(geo::GeometryCore const* geom):
        cryostat_id_iterator_base(geom, begin_pos) {}

      /// Constructor: points to the specified cryostat
      cryostat_id_iterator_base
        (geo::GeometryCore const* geom, GeoID_t const& start_from):
        cryostat_id_iterator_base(geom, undefined_pos)
        { id = start_from; }

      /// Constructor: points to begin
      cryostat_id_iterator_base
        (geo::GeometryCore const* geom, BeginPos_t const):
        cryostat_id_iterator_base(geom, undefined_pos)
        { set_begin(); }

      /// Constructor: points to end
      cryostat_id_iterator_base(geo::GeometryCore const* geom, EndPos_t):
        cryostat_id_iterator_base(geom, undefined_pos)
        { set_end(); }

      // TODO reconsider if the additional template is indeed needed
      /// Returns true if the two iterators point to the same cryostat
      template <typename OTHERID>
      bool operator== (cryostat_id_iterator_base<OTHERID> const& as) const
        { return localID() == as.localID(); }

      /// Returns true if the two iterators point to different cryostats
      template <typename OTHERID>
      bool operator!= (cryostat_id_iterator_base<OTHERID> const& as) const
        { return localID() != as.localID(); }

      /// Returns the ID the iterator points to
      reference operator* () const { return localID(); }

      /// Returns a pointer to the ID the iterator points to
      pointer operator-> () const { return &(localID()); }

      /// Prefix increment: returns this iterator pointing to the next cryostat
      iterator& operator++ () { next(); return *this; }

      /// Postfix increment: returns the current iterator, then increments it
      iterator operator++ (int) { iterator old(*this); next(); return old; }

      /// Returns whether the iterator is pointing to a valid cryostat
      operator bool() const;

      /// Returns a pointer to cryostat, or nullptr if invalid
      ElementPtr_t get() const;

        protected:
      using ID_t = typename LocalID_t::CryostatID_t;

      /// Constructor: does not set the current ID
      cryostat_id_iterator_base(geo::GeometryCore const* geom, UndefinedPos_t):
        geometry_iterator_base(geom), id()
        { set_local_limits(); }

      //@{
      /// Returns the actual type of ID we store
      GeoID_t const& ID() const { return id; }
      GeoID_t& ID() { return id; }
      //@}

      /// Skips to the next cryostat
      void next();

      /// Returns whether this iterator has reached the end
      bool at_end() const { return local_index() == limit; }

        private:
      GeoID_t id; ///< ID of the current cryostat
      ID_t limit = LocalID_t::InvalidID; ///< maximum number of cryostats

      /// Sets the limit member to the past-the-end cryostat number
      void set_local_limits();

      /// Sets the iterator to the begin position
      void set_begin();

      /// Sets the iterator to the end position
      void set_end();

      //@{
      /// Returns the type of ID we act on
      LocalID_t const& localID() const
        { return static_cast<LocalID_t const&>(ID()); }
      LocalID_t& localID() { return static_cast<LocalID_t&>(ID()); }
      //@}

      //@{
      /// Returns the index (part if the ID) this iterator runs on
      ID_t const& local_index() const { return localID().Cryostat; }
      ID_t& local_index() { return localID().Cryostat; }
      //@}

    }; // class cryostat_id_iterator_base<>


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
    template <typename GEOID>
    class TPC_id_iterator_base:
      virtual public std::forward_iterator_tag,
      protected cryostat_id_iterator_base<GEOID>
    {
      using upper_iterator = cryostat_id_iterator_base<GEOID>;

        public:
      using ElementPtr_t = geo::TPCGeo const*;
      using GeoID_t = typename upper_iterator::GeoID_t;

      using LocalID_t = geo::TPCID; ///< type of the ID we change
      static_assert(std::is_base_of<LocalID_t, GEOID>::value,
        "template type GEOID is not a LocalID_t");

      using iterator = TPC_id_iterator_base<GeoID_t>; ///< type of this iterator

      // import all the useful types from the base templated class
      using typename upper_iterator::UndefinedPos_t;
      using typename upper_iterator::BeginPos_t;
      using typename upper_iterator::EndPos_t;

      // import all the useful members from the base templated class
      using upper_iterator::undefined_pos;
      using upper_iterator::begin_pos;
      using upper_iterator::end_pos;


      /// @name Iterator traits
      /// @{
      using difference_type = std::ptrdiff_t;
      using value_type = LocalID_t;
      using reference = value_type const&;
      using pointer  = value_type const*;
      using iterator_category = std::input_iterator_tag;
      /// @}


      /// Default constructor; effect not defined: assign to it before using!
      TPC_id_iterator_base() {}

      /// Constructor: points to begin
      TPC_id_iterator_base(geo::GeometryCore const* geom):
        TPC_id_iterator_base(geom, begin_pos) {}

      /// Constructor: points to the specified TPC
      TPC_id_iterator_base
        (geo::GeometryCore const* geom, GeoID_t const& start_from):
        upper_iterator(geom, start_from)
        { set_local_limits(); }

      /// Constructor: points to begin
      TPC_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t const):
        upper_iterator(geom, begin_pos)
        { set_local_limits(); }

      /// Constructor: points to end
      TPC_id_iterator_base(geo::GeometryCore const* geom, EndPos_t):
        upper_iterator(geom, end_pos)
        {} // the local limit is ill-defined and left invalid

      // TODO reconsider if the additional template is indeed needed
      /// Returns true if the two iterators point to the same TPC
      template <typename OTHERID>
      bool operator== (TPC_id_iterator_base<OTHERID> const& as) const
        { return localID() == as.localID(); }

      /// Returns true if the two iterators point to different TPCs
      template <typename OTHERID>
      bool operator!= (TPC_id_iterator_base<OTHERID> const& as) const
        { return localID() != as.localID(); }

      /// Returns the TPCID the iterator points to
      reference operator* () const { return localID(); }

      /// Returns the TPCID the iterator points to
      pointer operator-> () const { return &(localID()); }

      /// Prefix increment: returns this iterator pointing to the next TPC
      iterator& operator++ () { next(); return *this; }

      /// Postfix increment: returns the current iterator, then increments it
      iterator operator++ (int) { iterator old(*this); next(); return old; }

      /// Returns whether the iterator is pointing to a valid TPC
      operator bool() const;

      /// Returns a pointer to TPC, or nullptr if invalid
      ElementPtr_t get() const;

        protected:

      using ID_t = typename LocalID_t::TPCID_t; ///< specific type for TPC ID

      /// Constructor: position undefined (meaning undefined local limits too)
      TPC_id_iterator_base(geo::GeometryCore const* geom, UndefinedPos_t):
        upper_iterator(geom, undefined_pos)
        {}

      using upper_iterator::ID; // to be explicit; this is NOT overloaded

      /// Returns the type of ID we act on
      LocalID_t const& localID() const
        { return static_cast<LocalID_t const&>(upper_iterator::ID()); }

      using upper_iterator::at_end; // to be explicit; this is NOT overloaded

      /// Skips to the next TPC
      void next();

      /// Returns the index (part if the ID) this iterator runs on
      ID_t const& local_index() const { return localID().TPC; }

        private:

      /// maximum number of TPCs in the current cryostat
      ID_t limit = LocalID_t::InvalidID;

      /// Sets limit to the past-the-end TPC number of current croystat
      void set_local_limits();

      /// Returns the type of ID we act on (non-const version)
      LocalID_t& localID() { return static_cast<LocalID_t&>(ID()); }

      /// Returns the index (part if the ID) this iterator runs on  (non-const)
      ID_t& local_index() { return localID().TPC; }

    }; // class TPC_id_iterator_base


    /**
     * @brief Base forward iterator browsing all plane IDs in the detector
     * @tparam GEOID ID type to be used
     *
     * This iterator requires that GEOID is derived from geo::PlaneID.
     * Note that no polymorphic behaviour is required, or expected, from GEOID.
     *
     * This iterator is designed to carry on, untouched, anything else that the
     * GEOID type defines beyond the required PlaneID data.
     *
     * @note A number of "local" methods are overloaded: since there is no
     * polymorphism here and they are not virtual functions, these are designed
     * not to replace the inherited methods except within the non-inherited and
     * explicitly redefined methods.
     *
     * Currently, backward iterations are not supported.
     */
    template <typename GEOID>
    class plane_id_iterator_base:
      virtual public std::forward_iterator_tag,
      protected TPC_id_iterator_base<GEOID>
    {
      using upper_iterator = TPC_id_iterator_base<GEOID>;

        public:
      using ElementPtr_t = geo::PlaneGeo const*;
      using GeoID_t = typename upper_iterator::GeoID_t;

      using LocalID_t = geo::PlaneID; ///< type of the ID we change
      static_assert(std::is_base_of<LocalID_t, GEOID>::value,
        "template type GEOID is not a LocalID_t");

      /// type of this iterator
      using iterator = plane_id_iterator_base<GeoID_t>;

      // import all the useful types from the base templated class
      using typename upper_iterator::UndefinedPos_t;
      using typename upper_iterator::BeginPos_t;
      using typename upper_iterator::EndPos_t;

      // import all the useful members from the base templated class
      using upper_iterator::undefined_pos;
      using upper_iterator::begin_pos;
      using upper_iterator::end_pos;


      /// @name Iterator traits
      /// @{
      using difference_type = std::ptrdiff_t;
      using value_type = LocalID_t;
      using reference = value_type const&;
      using pointer  = value_type const*;
      using iterator_category = std::input_iterator_tag;
      /// @}


      /// Default constructor; effect not defined: assign to it before using!
      plane_id_iterator_base() {}

      /// Constructor: points to begin
      plane_id_iterator_base(geo::GeometryCore const* geom):
        plane_id_iterator_base(geom, begin_pos) {}

      /// Constructor: points to the specified plane
      plane_id_iterator_base
        (geo::GeometryCore const* geom, GeoID_t const& start_from):
        upper_iterator(geom, start_from)
        { set_local_limits(); }

      /// Constructor: points to begin
      plane_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t const):
        upper_iterator(geom, begin_pos)
        { set_local_limits(); }

      /// Constructor: points to end
      plane_id_iterator_base(geo::GeometryCore const* geom, EndPos_t):
        upper_iterator(geom, end_pos)
        {} // the local limit is ill-defined and left invalid

      // TODO reconsider if the additional template is indeed needed
      /// Returns true if the two iterators point to the same plane
      template <typename OTHERID>
      bool operator== (plane_id_iterator_base<OTHERID> const& as) const
        { return localID() == as.localID(); }

      /// Returns true if the two iterators point to different planes
      template <typename OTHERID>
      bool operator!= (plane_id_iterator_base<OTHERID> const& as) const
        { return localID() != as.localID(); }

      /// Returns the PlaneID the iterator points to
      reference operator* () const { return localID(); }

      /// Returns the PlaneID the iterator points to
      pointer operator-> () const { return &(localID()); }

      /// Prefix increment: returns this iterator pointing to the next plane
      iterator& operator++ () { next(); return *this; }

      /// Postfix increment: returns the current iterator, then increments it
      iterator operator++ (int) { iterator old(*this); next(); return old; }

      /// Returns whether the iterator is pointing to a valid plane
      operator bool() const;

      /// Returns a pointer to plane, or nullptr if invalid
      ElementPtr_t get() const;

        protected:

      using ID_t = typename LocalID_t::PlaneID_t; ///< specific type for plane ID

      /// Constructor: position undefined (meaning undefined local limits too)
      plane_id_iterator_base(geo::GeometryCore const* geom, UndefinedPos_t):
        upper_iterator(geom, undefined_pos)
        {}

      using upper_iterator::ID; // to be explicit; this is NOT overloaded

      /// Returns the type of ID we act on
      LocalID_t const& localID() const
        { return static_cast<LocalID_t const&>(upper_iterator::ID()); }

      using upper_iterator::at_end; // to be explicit; this is NOT overloaded

      /// Skips to the next plane
      void next();

      /// Returns the index (part if the ID) this iterator runs on
      ID_t const& local_index() const { return localID().Plane; }

        private:

      /// maximum number of planes in the current TPC
      ID_t limit = LocalID_t::InvalidID;

      /// Sets limit to the past-the-end plane number of current TPC
      void set_local_limits();

      /// Returns the type of ID we act on (non-const version)
      LocalID_t& localID() { return static_cast<LocalID_t&>(ID()); }

      /// Returns the index (part if the ID) this iterator runs on  (non-const)
      ID_t& local_index() { return localID().Plane; }

    }; // class plane_id_iterator_base


    /**
     * @brief Base forward iterator browsing all wire IDs in the detector
     * @tparam GEOID ID type to be used
     *
     * This iterator requires that GEOID is derived from geo::WireID.
     * Note that no polymorphic behaviour is required, or expected, from GEOID.
     *
     * This iterator is designed to carry on, untouched, anything else that the
     * GEOID type defines beyond the required WireID data.
     *
     * @note A number of "local" methods are overloaded: since there is no
     * polymorphism here and they are not virtual functions, these are designed
     * not to replace the inherited methods except within the non-inherited and
     * explicitly redefined methods.
     *
     * Currently, backward iterations are not supported.
     */
    template <typename GEOID>
    class wire_id_iterator_base:
      virtual public std::forward_iterator_tag,
      protected plane_id_iterator_base<GEOID>
    {
      using upper_iterator = plane_id_iterator_base<GEOID>;

        public:
      using ElementPtr_t = geo::WireGeo const*;
      using GeoID_t = typename upper_iterator::GeoID_t;

      using LocalID_t = geo::WireID; ///< type of the ID we change
      static_assert(std::is_base_of<LocalID_t, GEOID>::value,
        "template type GEOID is not a LocalID_t");

      /// type of this iterator
      using iterator = wire_id_iterator_base<GeoID_t>;

      // import all the useful types from the base templated class
      using typename upper_iterator::UndefinedPos_t;
      using typename upper_iterator::BeginPos_t;
      using typename upper_iterator::EndPos_t;

      // import all the useful members from the base templated class
      using upper_iterator::undefined_pos;
      using upper_iterator::begin_pos;
      using upper_iterator::end_pos;


      /// @name Iterator traits
      /// @{
      using difference_type = std::ptrdiff_t;
      using value_type = LocalID_t;
      using reference = value_type const&;
      using pointer  = value_type const*;
      using iterator_category = std::input_iterator_tag;
      /// @}


      /// Default constructor; effect not defined: assign to it before using!
      wire_id_iterator_base() {}

      /// Constructor: points to begin
      wire_id_iterator_base(geo::GeometryCore const* geom):
        wire_id_iterator_base(geom, begin_pos) {}

      /// Constructor: points to the specified wire
      wire_id_iterator_base
        (geo::GeometryCore const* geom, GeoID_t const& start_from):
        upper_iterator(geom, start_from)
        { set_local_limits(); }

      /// Constructor: points to begin
      wire_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t const):
        upper_iterator(geom, begin_pos)
        { set_local_limits(); }

      /// Constructor: points to end
      wire_id_iterator_base(geo::GeometryCore const* geom, EndPos_t):
        upper_iterator(geom, end_pos)
        {} // the local limit is ill-defined and left invalid

      // TODO reconsider if the additional template is indeed needed
      /// Returns true if the two iterators point to the same wire
      template <typename OTHERID>
      bool operator== (wire_id_iterator_base<OTHERID> const& as) const
        { return localID() == as.localID(); }

      /// Returns true if the two iterators point to different wires
      template <typename OTHERID>
      bool operator!= (wire_id_iterator_base<OTHERID> const& as) const
        { return localID() != as.localID(); }

      /// Returns the WireID the iterator points to
      reference operator* () const { return localID(); }

      /// Returns the WireID the iterator points to
      pointer operator-> () const { return &(localID()); }

      /// Prefix increment: returns this iterator pointing to the next wire
      iterator& operator++ () { next(); return *this; }

      /// Postfix increment: returns the current iterator, then increments it
      iterator operator++ (int) { iterator old(*this); next(); return old; }

      /// Returns whether the iterator is pointing to a valid wire
      operator bool() const;

      /// Returns a pointer to wire, or nullptr if invalid
      ElementPtr_t get() const;

        protected:

      using ID_t = typename LocalID_t::WireID_t; ///< specific type for wire ID

      /// Constructor: position undefined (meaning undefined local limits too)
      wire_id_iterator_base(geo::GeometryCore const* geom, UndefinedPos_t):
        upper_iterator(geom, undefined_pos)
        {}

      using upper_iterator::ID; // to be explicit; this is NOT overloaded

      /// Returns the type of ID we act on
      LocalID_t const& localID() const
        { return static_cast<LocalID_t const&>(upper_iterator::ID()); }

      using upper_iterator::at_end; // to be explicit; this is NOT overloaded

      /// Skips to the next wire
      void next();

      /// Returns the index (part if the ID) this iterator runs on
      ID_t const& local_index() const { return localID().Wire; }

        private:

      /// maximum number of wires in the current plane
      ID_t limit = LocalID_t::InvalidID;

      /// Sets limit to the past-the-end wire number of current plane
      void set_local_limits();

      /// Returns the type of ID we act on (non-const version)
      LocalID_t& localID() { return static_cast<LocalID_t&>(ID()); }

      /// Returns the index (part if the ID) this iterator runs on (non-const)
      ID_t& local_index() { return localID().Wire; }

    }; // class wire_id_iterator_base


    // forward declarations:
    template <typename GEOIDITER>
    class geometry_element_iterator;

    /// Comparison operator: geometry ID and element point to the same ID.
    template <typename GEOIDITER>
    bool operator== (
      geometry_element_iterator<GEOIDITER> const& iter,
      GEOIDITER const& id_iter
      );
    /// Comparison operator: geometry ID and element point to the same ID.
    template <typename GEOIDITER>
    inline bool operator== (
      GEOIDITER const& id_iter,
      geometry_element_iterator<GEOIDITER> const& iter
      )
      { return iter == id_iter; }

    /// Comparison operator: geometry ID and element point to different IDs.
    template <typename GEOIDITER>
    bool operator!= (
      geometry_element_iterator<GEOIDITER> const& iter,
      GEOIDITER const& id_iter
      );
    /// Comparison operator: geometry ID and element point to different IDs.
    template <typename GEOIDITER>
    inline bool operator!= (
      GEOIDITER const& id_iter,
      geometry_element_iterator<GEOIDITER> const& iter
      )
      { return iter != id_iter; }

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
    template <typename GEOIDITER>
    class geometry_element_iterator:
      public std::forward_iterator_tag, public geometry_iterator_types
    {
        public:
      using id_iterator_t = GEOIDITER;

      static_assert(
        std::is_base_of<geometry_iterator_base, id_iterator_t>::value,
        "template class for geometry_element_iterator"
        " must be a geometry iterator"
        );

      using iterator = geometry_element_iterator<id_iterator_t>; ///< this type

      /// @{
      /// @name Types mirrored from the ID iterator
      using LocalID_t = typename id_iterator_t::LocalID_t;
      using GeoID_t = typename id_iterator_t::GeoID_t;
      using UndefinedPos_t = typename id_iterator_t::UndefinedPos_t;
      using BeginPos_t = typename id_iterator_t::BeginPos_t;
      using EndPos_t = typename id_iterator_t::EndPos_t;
      using ElementPtr_t = typename id_iterator_t::ElementPtr_t;
      /// @}

      /// @{
      /// @name Constants inherited from the ID iterator
      using geometry_iterator_types::undefined_pos;
      using geometry_iterator_types::begin_pos;
      using geometry_iterator_types::end_pos;
      /// @}

      /// Geometry class pointed by the iterator
      using Element_t = typename std::remove_pointer<ElementPtr_t>::type;


      /// @name Iterator traits
      /// @{
      using difference_type = std::ptrdiff_t;
      using value_type = Element_t;
      using reference = value_type const&;
      using pointer  = value_type const*;
      using iterator_category = std::forward_iterator_tag;
      /// @}


      /// Default constructor; effect not defined: assign to it before using!
      geometry_element_iterator() = default;

      /// Constructor: points to begin
      geometry_element_iterator(geo::GeometryCore const* geom):
        id_iter(geom) {}

      /// Constructor: points to the same element as the specified ID iterator.
      geometry_element_iterator(id_iterator_t const& iter): id_iter(iter) {}

      /// Constructor: points to the same element as the specified ID iterator.
      geometry_element_iterator(id_iterator_t&& iter): id_iter(iter) {}

      /// Constructor: points to the specified geometry element
      geometry_element_iterator
        (geo::GeometryCore const* geom, GeoID_t const& start_from):
        id_iter(geom, start_from)
        {}

      /// Constructor: points to beginning
      geometry_element_iterator
        (geo::GeometryCore const* geom, BeginPos_t const pos):
        id_iter(geom, pos)
        {}

      /// Constructor: points to end
      geometry_element_iterator
        (geo::GeometryCore const* geom, EndPos_t const pos):
        id_iter(geom, pos)
        {}

      /// Returns true if the two iterators point to the same object
      bool operator== (iterator const& as) const
        { return id_iterator() == as.id_iterator(); }

      /// Returns true if the two iterators point to different objects
      bool operator!= (iterator const& as) const
        { return id_iterator() != as.id_iterator(); }

      /**
       * @brief Returns the geometry element the iterator points to
       * @return a constant reference to the element the iterator points to
       * @throw cet::exception (category "geometry_iterator") if no valid
       *   geometry element is currently pointed by the iterator
       */
      reference operator* () const
        {
          ElementPtr_t ptr = get();
          if (ptr) return *ptr;
          throw cet::exception("geometry_iterator")
            << "iterator attempted to obtain geometry element "
            << std::string(ID());
        } // operator*()

      /// Returns a pointer to the element the iterator points to (or nullptr)
      pointer operator-> () const { return get(); }

      /// Prefix increment: returns this iterator pointing to the next element
      iterator& operator++ () { ++id_iterator(); return *this; }

      /// Postfix increment: returns the current iterator, then increments it
      iterator operator++ (int)
        { iterator old(*this); ++id_iterator(); return old; }

      /// Returns whether the iterator is pointing to a valid geometry element
      operator bool() const
        { return bool(id_iterator()) && (id_iterator().get() != nullptr); }

      /// Returns a pointer to the geometry element, or nullptr if invalid
      ElementPtr_t get() const { return id_iterator().get(); }

      /// Returns the ID of the pointed geometry element
      LocalID_t const& ID() const { return *(id_iterator()); }

        protected:
      friend bool geo::details::operator== <id_iterator_t>
        (iterator const& iter, id_iterator_t const& id_iter);
      friend bool geo::details::operator== <id_iterator_t>
        (id_iterator_t const& id_iter, iterator const& iter);
      friend bool geo::details::operator!= <id_iterator_t>
        (iterator const& iter, id_iterator_t const& id_iter);
      friend bool geo::details::operator!= <id_iterator_t>
        (id_iterator_t const& id_iter, iterator const& iter);

      //@{
      /// Access to the base ID iterator
      id_iterator_t const& id_iterator() const { return id_iter; }
      id_iterator_t& id_iterator() { return id_iter; }
      //@}

        private:
      id_iterator_t id_iter; ///< iterator performing the job

    }; // class geometry_element_iterator<>


    /**
     * @brief Base forward iterator browsing all TPC set IDs in the detector.
     * @tparam GEOID ID type to be used
     *
     * This iterator requires that GEOID is derived from geo::TPCSetID.
     * Note that no polymorphic behaviour is required, or expected, from GEOID.
     *
     * This iterator is designed to carry on, untouched, anything else that the
     * GEOID type defines beyond the required TPCsetID data.
     *
     * @note A number of "local" methods are overloaded: since there is no
     * polymorphism here and they are not virtual functions, these are designed
     * not to replace the inherited methods except within the non-inherited and
     * explicitly redefined methods.
     *
     * Currently, backward iterations are not supported.
     */
    template <typename GEOID>
    class TPCset_id_iterator_base:
      virtual public std::forward_iterator_tag,
      protected cryostat_id_iterator_base<GEOID>
    {
      using upper_iterator = cryostat_id_iterator_base<GEOID>;

        public:
      using GeoID_t = typename upper_iterator::GeoID_t;

      using LocalID_t = readout::TPCsetID; ///< Type of the ID we change.
      static_assert(std::is_base_of<LocalID_t, GEOID>::value,
        "template type GEOID is not a LocalID_t");

      ///< Type of this iterator.
      using iterator = TPCset_id_iterator_base<GeoID_t>;

      // import all the useful types from the base templated class
      using typename upper_iterator::UndefinedPos_t;
      using typename upper_iterator::BeginPos_t;
      using typename upper_iterator::EndPos_t;

      // import all the useful members from the base templated class
      using upper_iterator::undefined_pos;
      using upper_iterator::begin_pos;
      using upper_iterator::end_pos;


      /// @name Iterator traits
      /// @{
      using difference_type = std::ptrdiff_t;
      using value_type = LocalID_t;
      using reference = value_type const&;
      using pointer  = value_type const*;
      using iterator_category = std::input_iterator_tag;
      /// @}


      /// Default constructor; effect not defined: assign to it before using!
      TPCset_id_iterator_base() {}

      /// Constructor: points to begin.
      TPCset_id_iterator_base(geo::GeometryCore const* geom)
        : TPCset_id_iterator_base(geom, begin_pos)
        {}

      /// Constructor: points to the specified TPC set.
      TPCset_id_iterator_base
        (geo::GeometryCore const* geom, GeoID_t const& start_from)
        : upper_iterator(geom, start_from)
        { set_local_limits(); }

      /// Constructor: points to begin.
      TPCset_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t const)
        : upper_iterator(geom, begin_pos)
        { set_local_limits(); }

      /// Constructor: points to end.
      TPCset_id_iterator_base(geo::GeometryCore const* geom, EndPos_t)
        : upper_iterator(geom, end_pos)
        {} // the local limit is ill-defined and left invalid

      // TODO reconsider if the additional template is indeed needed
      /// Returns true if the two iterators point to the same TPC set.
      template <typename OTHERID>
      bool operator== (TPCset_id_iterator_base<OTHERID> const& as) const
        { return localID() == as.localID(); }

      /// Returns true if the two iterators point to different TPC sets.
      template <typename OTHERID>
      bool operator!= (TPCset_id_iterator_base<OTHERID> const& as) const
        { return localID() != as.localID(); }

      /// Returns the TPCsetID the iterator points to.
      reference operator* () const { return localID(); }

      /// Returns the TPCsetID the iterator points to.
      pointer operator-> () const { return &(localID()); }

      /// Prefix increment: returns this iterator pointing to the next TPC set.
      iterator& operator++ () { next(); return *this; }

      /// Postfix increment: returns the current iterator, then increments it.
      iterator operator++ (int) { iterator old(*this); next(); return old; }

      /// Returns whether the iterator is pointing to a valid TPC set.
      operator bool() const;

        protected:

      /// Specific type for TPC set ID.
      using ID_t = typename LocalID_t::TPCsetID_t;

      /// Constructor: position undefined (meaning undefined local limits too).
      TPCset_id_iterator_base(geo::GeometryCore const* geom, UndefinedPos_t)
        : upper_iterator(geom, undefined_pos)
        {}

      using upper_iterator::ID; // to be explicit; this is NOT overloaded

      /// Returns the type of ID we act on.
      LocalID_t const& localID() const
        { return static_cast<LocalID_t const&>(upper_iterator::ID()); }

      using upper_iterator::at_end; // to be explicit; this is NOT overloaded

      /// Skips to the next TPC set.
      void next();

      /// Returns the index (part if the ID) this iterator runs on.
      ID_t const& local_index() const { return localID().TPCset; }

        private:

      /// maximum number of TPC sets in the current cryostat.
      ID_t limit = LocalID_t::InvalidID;

      /// Sets limit to the past-the-end TPC set number of current croystat.
      void set_local_limits();

      /// Returns the type of ID we act on (non-const version).
      LocalID_t& localID() { return static_cast<LocalID_t&>(ID()); }

      /// Returns the index (part if the ID) this iterator runs on  (non-const).
      ID_t& local_index() { return localID().TPCset; }

      // no object is currently implemented for TPC sets
      typename upper_iterator::ElementPtr_t get() const = delete;


    }; // class TPCset_id_iterator_base


    /**
     * @brief Base forward iterator browsing all readout plane IDs in the
     *        detector
     * @tparam GEOID ID type to be used
     *
     * This iterator requires that GEOID is derived from geo::ROPID.
     * Note that no polymorphic behaviour is required, or expected, from GEOID.
     *
     * This iterator is designed to carry on, untouched, anything else that the
     * GEOID type defines beyond the required ROPID data.
     *
     * @note A number of "local" methods are overloaded: since there is no
     * polymorphism here and they are not virtual functions, these are designed
     * not to replace the inherited methods except within the non-inherited and
     * explicitly redefined methods.
     *
     * Currently, backward iterations are not supported.
     */
    template <typename GEOID>
    class ROP_id_iterator_base:
      virtual public std::forward_iterator_tag,
      protected TPCset_id_iterator_base<GEOID>
    {
      using upper_iterator = TPCset_id_iterator_base<GEOID>;

        public:
      using GeoID_t = typename upper_iterator::GeoID_t;

      using LocalID_t = readout::ROPID; ///< type of the ID we change
      static_assert(std::is_base_of<LocalID_t, GEOID>::value,
        "template type GEOID is not a LocalID_t");

      /// Type of this iterator.
      using iterator = ROP_id_iterator_base<GeoID_t>;

      // import all the useful types from the base templated class
      using typename upper_iterator::UndefinedPos_t;
      using typename upper_iterator::BeginPos_t;
      using typename upper_iterator::EndPos_t;

      // import all the useful members from the base templated class
      using upper_iterator::undefined_pos;
      using upper_iterator::begin_pos;
      using upper_iterator::end_pos;


      /// @name Iterator traits
      /// @{
      using difference_type = std::ptrdiff_t;
      using value_type = LocalID_t;
      using reference = value_type const&;
      using pointer  = value_type const*;
      using iterator_category = std::input_iterator_tag;
      /// @}


      /// Default constructor; effect not defined: assign to it before using!
      ROP_id_iterator_base() = default;

      /// Constructor: points to begin.
      ROP_id_iterator_base(geo::GeometryCore const* geom)
        : ROP_id_iterator_base(geom, begin_pos) {}

      /// Constructor: points to the specified readout plane.
      ROP_id_iterator_base
        (geo::GeometryCore const* geom, GeoID_t const& start_from)
        : upper_iterator(geom, start_from)
        { set_local_limits(); }

      /// Constructor: points to begin.
      ROP_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t const)
        : upper_iterator(geom, begin_pos)
        { set_local_limits(); }

      /// Constructor: points to end.
      ROP_id_iterator_base(geo::GeometryCore const* geom, EndPos_t)
        : upper_iterator(geom, end_pos)
        {} // the local limit is ill-defined and left invalid

      // TODO reconsider if the additional template is indeed needed
      /// Returns true if the two iterators point to the same readout plane.
      template <typename OTHERID>
      bool operator== (ROP_id_iterator_base<OTHERID> const& as) const
        { return localID() == as.localID(); }

      /// Returns true if the two iterators point to different readout planes.
      template <typename OTHERID>
      bool operator!= (ROP_id_iterator_base<OTHERID> const& as) const
        { return localID() != as.localID(); }

      /// Returns the PlaneID the iterator points to
      reference operator* () const { return localID(); }

      /// Returns the PlaneID the iterator points to
      pointer operator-> () const { return &(localID()); }

      /// Prefix increment: returns this iterator pointing to the next plane
      iterator& operator++ () { next(); return *this; }

      /// Postfix increment: returns the current iterator, then increments it.
      iterator operator++ (int) { iterator old(*this); next(); return old; }

      /// Returns whether the iterator is pointing to a valid plane.
      operator bool() const;

        protected:

      using ID_t = typename LocalID_t::ROPID_t; ///< Specific type for plane ID.

      /// Constructor: position undefined (meaning undefined local limits too).
      ROP_id_iterator_base(geo::GeometryCore const* geom, UndefinedPos_t)
        : upper_iterator(geom, undefined_pos)
        {}

      using upper_iterator::ID; // to be explicit; this is NOT overloaded

      /// Returns the type of ID we act on.
      LocalID_t const& localID() const
        { return static_cast<LocalID_t const&>(upper_iterator::ID()); }

      using upper_iterator::at_end; // to be explicit; this is NOT overloaded

      /// Skips to the next readout plane.
      void next();

      /// Returns the index (part if the ID) this iterator runs on.
      ID_t const& local_index() const { return localID().ROP; }

        private:

      /// Maximum number of readout planes in the current TPC set.
      ID_t limit = LocalID_t::InvalidID;

      /// Sets limit to the past-the-end readout plane number of current TPC
      /// set.
      void set_local_limits();

      /// Returns the type of ID we act on (non-const version).
      LocalID_t& localID() { return static_cast<LocalID_t&>(ID()); }

      /// Returns the index (part if the ID) this iterator runs on  (non-const).
      ID_t& local_index() { return localID().ROP; }

    }; // class ROP_id_iterator_base


  } // namespace details

  // BEGIN Geometry group ------------------------------------------------------
  /// @ingroup Geometry
  /// @{
  /**
   * @brief Forward iterator browsing all cryostats in the detector
   *
   * Prefer asking GeometryCore object for iterators rather than constructing
   * them anew: see geo::GeometryCore::cryostat_id_iterator for the recommended
   * usage.
   * Stand-alone example (not recommended):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * geo::GeometryCore::cryostat_id_iterator iCryostat,
   *   cbegin(geom, geo::cryostat_id_iterator::begin_pos),
   *   cend(geom, geo::cryostat_id_iterator::end_pos);
   * for (iCryostat = cbegin; iCryostat != cend; ++iCryostat) {
   *   geo::CryostatID const& cid = *iCryostat;
   *   geo::CryostatGeo const* pCryo = iCryostat.get();
   *   std::cout << "We are at: " << cid << std::endl;
   *   // ...
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  using cryostat_id_iterator
    = details::cryostat_id_iterator_base<geo::CryostatID>;

  /**
   * @brief Forward iterator browsing all cryostats in the detector
   *
   * The comments from cryostat_id_iterator are valid here as well.
   * This object has a different dereferenciation operator that obtains
   * the plane directly, or throws on failure.
   */
  using cryostat_iterator
    = details::geometry_element_iterator<cryostat_id_iterator>;


  /**
   * @brief Forward iterator browsing all TPCs in the detector
   *
   * Prefer asking the geometry object for iterators rather than constructing
   * them anew: see geo::GeometryCore::TPC_id_iterator for the recommended
   * usage.
   * Stand-alone example (not recommended):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * geo::GeometryCore::TPC_id_iterator iTPC,
   *   tbegin(geom, geo::TPC_id_iterator::begin_pos),
   *   tend(geom, geo::TPC_id_iterator::end_pos);
   * for (iTPC = tbegin; iTPC != tend; ++iTPC) {
   *   geo::TPCID const& tid = *iTPC;
   *   geo::TPCGeo const* pTPC = iTPC.get();
   *   std::cout << "We are at: " << tid << std::endl;
   *   // ...
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  using TPC_id_iterator = details::TPC_id_iterator_base<geo::TPCID>;

  /**
   * @brief Forward iterator browsing all TPCs in the detector
   *
   * The comments from TPC_id_iterator are valid here as well.
   * This object has a different dereferenciation operator that obtains
   * the TPC directly, or throws on failure.
   */
  using TPC_iterator = details::geometry_element_iterator<TPC_id_iterator>;


  /**
   * @brief Forward iterator browsing all planes in the detector
   *
   * Prefer asking the geometry object for iterators rather than constructing
   * them anew: see geo::GeometryCore::plane_id_iterator for the recommended
   * usage.
   * Stand-alone example (not recommended):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * geo::GeometryCore::plane_id_iterator iPlane,
   *   pbegin(geom, geo::plane_id_iterator::begin_pos),
   *   pend(geom, geo::plane_id_iterator::end_pos);
   * for (iPlane = pbegin; iPlane != pend; ++iPlane) {
   *   geo::PlaneID const& pid = *iPlane;
   *   geo::PlaneGeo const* pPlane = iPlane.get();
   *   std::cout << "We are at: " << pid << std::endl;
   *   // ...
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  using plane_id_iterator = details::plane_id_iterator_base<geo::PlaneID>;

  /**
   * @brief Forward iterator browsing all planes in the detector
   *
   * The comments from plane_id_iterator are valid here as well.
   * This object has a different dereferenciation operator that obtains
   * the plane directly, or throws on failure.
   */
  using plane_iterator = details::geometry_element_iterator<plane_id_iterator>;


  /**
   * @brief Forward iterator browsing all wires in the detector
   *
   * Prefer asking the geometry object for iterators rather than constructing
   * them anew: see geo::GeometryCore::wire_id_iterator for the recommended usage.
   * Stand-alone example (not recommended):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * geo::GeometryCore::wire_id_iterator iWire,
   *   wbegin(geom, geo::wire_id_iterator::begin_pos),
   *   wend(geom, geo::wire_id_iterator::end_pos);
   * for (iWire = wbegin; iWire != wend; ++iWire) {
   *   geo::WireID const& wid = *iWire;
   *   geo::WireGeo const* pWire = iWire.get();
   *   std::cout << "We are at: " << wid << std::endl;
   *   // ...
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  using wire_id_iterator = details::wire_id_iterator_base<geo::WireID>;

  /**
   * @brief Forward iterator browsing all wires in the detector
   *
   * The comments from wire_id_iterator are valid here as well.
   * This object has a different dereferenciation operator that obtains
   * the wire directly, or throws on failure.
   */
  using wire_iterator = details::geometry_element_iterator<wire_id_iterator>;


  /**
   * @brief Forward iterator browsing all TPC sets in the detector.
   *
   * Prefer asking the geometry object for iterators rather than constructing
   * them anew: see `geo::GeometryCore::TPCset_id_iterator` for the recommended
   * usage.
   * Stand-alone example (not recommended):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * geo::GeometryCore::TPCset_id_iterator iTPCset,
   *   tbegin(geom, geo::iterators::begin_pos),
   *   tend(geom, geo::iterators::end_pos);
   * for (iTPCset = tbegin; iTPCset != tend; ++iTPCset) {
   *   readout::TPCsetID const& tid = *iTPCset;
   *   std::cout << "We are at: " << tid << std::endl;
   *   // ...
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  using TPCset_id_iterator
    = details::TPCset_id_iterator_base<readout::TPCsetID>;


  /**
   * @brief Forward iterator browsing all readout planes in the detector.
   *
   * Prefer asking the geometry object for iterators rather than constructing
   * them anew: see geo::GeometryCore::ROP_id_iterator for the recommended
   * usage.
   * Stand-alone example (not recommended):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * geo::GeometryCore::ROP_id_iterator iROP,
   *   rbegin(geom, geo::iterators::begin_pos),
   *   rend(geom, geo::iterators::end_pos);
   * for (iROP = rbegin; iROP != rend; ++iROP) {
   *   readout::ROPID const& rid = *iROP;
   *   std::cout << "We are at: " << rid << std::endl;
   *   // ...
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  using ROP_id_iterator = details::ROP_id_iterator_base<readout::ROPID>;

  template <
    typename Iter,
    Iter (GeometryCore::*BeginFunc)() const,
    Iter (GeometryCore::*EndFunc)() const
    >
  class IteratorBox: public util::span<Iter> {
      public:

    IteratorBox(GeometryCore const* geom)
      : util::span<Iter>((geom->*BeginFunc)(), (geom->*EndFunc)())
      {}

  }; // IteratorBox<>


  template <
    typename Iter,
    typename GeoID,
    Iter (GeometryCore::*BeginFunc)(GeoID const&) const,
    Iter (GeometryCore::*EndFunc)(GeoID const&) const
    >
  class LocalIteratorBox: public util::span<Iter> {
      public:

    LocalIteratorBox(GeometryCore const* geom, GeoID const& ID)
      : util::span<Iter>((geom->*BeginFunc)(ID), (geom->*EndFunc)(ID))
      {}

  }; // LocalIteratorBox<>


  /// Namespace for geometry iterators.
  /// (currently quite depleted)
  namespace iterators {

    using BeginPos_t = details::geometry_iterator_types::BeginPos_t;
    using EndPos_t = details::geometry_iterator_types::EndPos_t;
    using UndefinedPos_t = details::geometry_iterator_types::UndefinedPos_t;

    constexpr auto begin_pos = details::geometry_iterator_types::begin_pos;
    constexpr auto end_pos = details::geometry_iterator_types::end_pos;
    constexpr auto undefined_pos = details::geometry_iterator_types::undefined_pos;

  } // namespace iterators


  //
  // GeometryCore
  //


  /** **************************************************************************
   * @brief Description of geometry of one entire detector
   *
   * @note All lengths are specified in centimetres
   *
   *
   * How to correctly instantiate a GeometryCore object
   * ---------------------------------------------------
   *
   * Instantiation is a multi-step procedure:
   * 1. construct a GeometryCore object (the "service provider"),
   *    with the full configuration; at this step, configuration is just stored
   * 2. load a geometry with GeometryCore::LoadGeometryFile();
   *    this loads the detector geometry information
   * 3. prepare a channel map algorithm object (might use for example
   *    GeometryCore::DetectorName() or the detector geometry from the
   *    newly created object, but any use of channel mapping related functions
   *    is forbidden and it would yield undefined behaviour (expected to be
   *    catastrophic)
   * 4. acquire the channel mapping algorithm with
   *    GeometryCore::ApplyChannelMap().
   *
   * Step 3 (creation of the channel mapping algorithm object) can be performed
   * at any time before step 4, provided that no GeometryCore instance is needed
   * for it.
   *
   *
   * Configuration parameters
   * -------------------------
   *
   * - *Name* (string; mandatory): string identifying the detector; it can be
   *   different from the base name of the file used to initialize the geometry;
   *   standard names are recommended by each experiment.
   *   This name can be used, for example, to select which channel mapping
   *   algorithm to use.
   * - *SurfaceY* (real; mandatory): depth of the detector, in centimetrs;
   *   see SurfaceY() for details
   * - *MinWireZDist* (real; default: 3)
   * - *PositionEpsilon* (real; default: 0.01%) set the default tolerance
   *   (see DefaultWiggle())
   *
   */
  class GeometryCore {

    using DefaultVector_t = TVector3; ///< Default template argument.
    using DefaultPoint_t = TVector3; ///< Default template argument.

  public:

    /// Type used for expressing coordinates
    /// @deprecated Use directly `geo::Length_t`
    using Coord_t [[deprecated("Use geo::Point_t instead")]] = geo::Length_t;

    /// Type used to represent a point in global coordinates
    /// @deprecated Use directly `geo::Point_t`
    using Point3D_t [[deprecated("Convert the code to use geo::Point_t")]]
      = DefaultPoint_t;


    /// Simple class with two points (a pair with aliases).
    template <typename Point>
    struct Segment: public std::pair<Point, Point> {

      // use the base class constructors
      using std::pair<Point, Point>::pair;

      Point const& start() const { return this->first; }
      Point& start() { return this->first; }

      Point const& end() const { return this->second; }
      Point& end() { return this->second; }

    }; // struct Segment_t

    using Segment_t = Segment<DefaultPoint_t>;

    /// Type of list of cryostats
    using CryostatList_t = GeometryData_t::CryostatList_t;
    /// Type of list of auxiliary detectors
    using AuxDetList_t = GeometryData_t::AuxDetList_t;


    /// Wires must be found in GDML description within this number of nested
    /// volumes.
    static constexpr std::size_t MaxWireDepthInGDML = 20U;

    /// Value of tolerance for equality comparisons
    static lar::util::RealComparisons<geo::Length_t> coordIs;


    // import iterators
    /**
     * @brief Forward-iterator browsing all cryostat IDs in the detector.
     *
     * Usage example with a while loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::GeometryCore::cryostat_id_iterator
     *   iCryostat = geom->begin_cryostat_id(), cend = geom->end_cryostat_id();
     * while (iCryostat != cend) {
     *   std::cout << "Cryo: " << iCryostat->Cryostat << std::endl;
     *   const geo::CryostatGeo* pCryo = iCryostat.get();
     *   ++iCryostat;
     *   // ...
     * } // while
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The recommended way to iterate is actually to use
     * `GeometryCore::IterateCryostatIDs()` in a range-for loop.
     * It is recommended to save the end iterator rather than calling
     * `GeometryCore::end_cryostat_id()` on every check.
     */
    using cryostat_id_iterator = geo::cryostat_id_iterator;

    /**
     * @brief Forward-iterator browsing all cryostats in the detector.
     *
     * Usage example with a while loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::GeometryCore::cryostat_iterator
     *   iCryostat = geom->begin_cryostat(), cend = geom->end_cryostat();
     * while (iCryostat != cend) {
     *   std::cout << "Cryo: " << iCryostat.ID() << std::endl;
     *   geo::CryostatGeo const& Cryo = *iCryostat;
     *   ++iCryostat;
     *   // ...
     * } // while
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The recommended way to iterate is actually to use
     * `GeometryCore::IterateCryostats()` in a range-for loop.
     * It is recommended to save the end iterator rather than calling
     * `GeometryCore::end_cryostat()` on every check.
     */
    using cryostat_iterator = geo::cryostat_iterator;

    /**
     * @brief Forward-iterator browsing all TPC IDs in the detector.
     *
     * Usage example with a while loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::GeometryCore::TPC_id_iterator iTPC = geom->begin_TPC_id(),
     *   tend = geom->end_TPC_id();
     * while (iTPC != tend) {
     *   std::cout << "TPC: " << *iTPC << std::endl;
     *   // the TPC descriptor object
     *   const geo::TPCGeo* pTPC = iTPC.get();
     *   // the cryostat the TPC is in
     *   geo::CryostatGeo const& Cryo = geom->Cryostat(*iTPC);
     *   ++iTPC;
     *   // ...
     * } // while
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The recommended way to iterate is actually to use
     * `GeometryCore::IterateTPCIDs()` in a range-for loop.
     * It is recommended to save the end iterator rather than calling
     * `GeometryCore::end_TPC_id()` on every check.
     */
    using TPC_id_iterator = geo::TPC_id_iterator;

    /**
     * @brief Forward-iterator browsing all TPCs in the detector.
     *
     * Usage example with a while loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::GeometryCore::TPC_iterator iTPC = geom->begin_TPC(),
     *   tend = geom->end_TPC();
     * while (iTPC != tend) {
     *   std::cout << "TPC: " << iTPC.ID() << std::endl;
     *   // the TPC descriptor object
     *   geo::TPCGeo const& TPC = *iTPC;
     *   ++iTPC;
     *   // ...
     * } // while
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The recommended way to iterate is actually to use
     * `GeometryCore::IterateTPCs()` in a range-for loop.
     * It is recommended to save the end iterator rather than calling
     * `GeometryCore::end_TPC()` on every check.
     */
    using TPC_iterator = geo::TPC_iterator;

    /**
     * @brief Forward-iterator browsing all plane IDs in the detector.
     *
     * Usage example with a while loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::GeometryCore::plane_id_iterator iPlane = geom->begin_plane_id(),
     *   pend = geom->end_plane_id();
     * while (iPlane != pend) {
     *   std::cout << "Plane: " << *iPlane << std::endl;
     *   // the plane descriptor object
     *   const geo::PlaneGeo* pPlane = iPlane.get();
     *   // the TPC the plane is in
     *   geo::TPCGeo const& TPC = geom->TPC(*iPlane);
     *   ++iPlane;
     *   // ...
     * } // while
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The recommended way to iterate is actually to use
     * `GeometryCore::IteratePlaneIDs()` in a range-for loop.
     * It is recommended to save the end iterator rather than calling
     * `GeometryCore::end_plane_id()` on every check.
     */
    using plane_id_iterator = geo::plane_id_iterator;

    /**
     * @brief Forward-iterator browsing all planes in the detector.
     *
     * Usage example with a while loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::GeometryCore::plane_iterator iPlane = geom->begin_plane(),
     *   pend = geom->end_plane();
     * while (iPlane != pend) {
     *   std::cout << "Plane: " << iPlane.ID() << std::endl;
     *   // the plane descriptor object
     *   geo::PlaneGeo const& Plane = *iPlane;
     *   ++iPlane;
     *   // ...
     * } // while
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The recommended way to iterate is actually to use
     * `GeometryCore::IteratePlanes()` in a range-for loop.
     * It is recommended to save the end iterator rather than calling
     * `GeometryCore::end_plane()` on every check.
     */
    using plane_iterator = geo::plane_iterator;

    /**
     * @brief Forward-iterator browsing all wire IDs in the detector.
     *
     * Usage example with a while loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::GeometryCore::wire_id_iterator iWire = geom->begin_wire_id(),
     *   wend = geom->end_wire_id();
     * while (iWire != wend) {
     *   std::cout << "Wire: " << *iWire << std::endl;
     *   // the wire descriptor object
     *   const geo::WireGeo* pWire = iWire.get();
     *   // the TPC the wire is in
     *   geo::TPCGeo const& TPC = geom->TPC(*iWire);
     *   ++iWire;
     *   // ...
     * } // while
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The recommended way to iterate is actually to use
     * `GeometryCore::IterateWireIDs()` in a range-for loop.
     * It is recommended to save the end iterator rather than calling
     * `GeometryCore::end_wire_id()` on every check.
     */
    using wire_id_iterator = geo::wire_id_iterator;

    /**
     * @brief Forward-iterator browsing all wires in the detector.
     *
     * Usage example with a while loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::GeometryCore::wire_iterator iWire = geom->begin_wire(),
     *   wend = geom->end_wire();
     * while (iWire != wend) {
     *   std::cout << "Wire: " << iWire.ID() << std::endl;
     *   // the wire descriptor object
     *   geo::WireGeo const& Wire = *iWire;
     *   ++iWire;
     *   // ...
     * } // while
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The recommended way to iterate is actually to use
     * `GeometryCore::IterateWires()` in a range-for loop.
     * It is recommended to save the end iterator rather than calling
     * `GeometryCore::end_wire()` on every check.
     */
    using wire_iterator = geo::wire_iterator;



    /**
     * @brief Initialize geometry from a given configuration
     * @param pset configuration parameters
     *
     * This constructor does not load any geometry description.
     * The next step is to do exactly that, by GeometryCore::LoadGeometryFile().
     */
    GeometryCore(fhicl::ParameterSet const& pset);

    /// Destructor
    ~GeometryCore();

    // this object is not copiable nor moveable (see also issue #14384);
    // currently, auxiliary detectors are stored as bare pointers,
    // which prevents trivial copy or move.
    GeometryCore(GeometryCore const&) = delete;
    GeometryCore(GeometryCore&&) = delete;
    GeometryCore& operator= (GeometryCore const&) = delete;
    GeometryCore& operator= (GeometryCore&&) = delete;


    /**
     * @brief Returns the tolerance used in looking for positions
     * @return the tolerance value
     *
     * This parameter is used as tolerance ("wiggle") for methods that require
     * it (e.g. `geo::CryostatGeo::FindTPCAtPosition()`).
     * Typically, it's a additional fraction of tolerance: 0 means no tolerance,
     * 0.1 means 10% tolerance.
     *
     * @todo Confirm the definition of wiggle: this one is taken from other doc
     */
    double DefaultWiggle() const { return fPositionWiggle; }

    /**
     * @brief Returns the full directory path to the geometry file source
     * @return the full directory path to the geometry file source
     *
     * This is the full path of the source of the detector geometry GeometryCore
     * relies on.
     */
    std::string ROOTFile() const { return fROOTfile; }

    /**
     * @brief Returns the full directory path to the GDML file source
     * @return the full directory path to the GDML file source
     *
     * This is the full path of the source of the detector geometry handed to
     * the detector simulation (GEANT).
     */
    std::string GDMLFile() const { return fGDMLfile; }



    // BEGIN Detector information
    /// @name Detector information
    /// @{

    //
    // global features
    //
    /// Returns a string with the name of the detector, as configured
    std::string DetectorName() const { return fDetectorName; }


    //
    // position
    //

    /// Returns a pointer to the world volume.
    TGeoVolume const* WorldVolume() const;


    /**
     * @brief Fills the arguments with the boundaries of the world
     * @param xlo (output) pointer to the lower x coordinate
     * @param xlo (output) pointer to the upper x coordinate
     * @param ylo (output) pointer to the lower y coordinate
     * @param ylo (output) pointer to the upper y coordinate
     * @param zlo (output) pointer to the lower z coordinate
     * @param zlo (output) pointer to the upper z coordinate
     * @throw cet::exception (`"GeometryCore"` category) if no world found
     * @see `GetWorldVolumeName()`
     *
     * This method fills the boundaries of the world volume
     * (`GetWorldVolumeName()`).
     *
     * If a pointer is null, its coordinate is skipped.
     *
     * @deprecated Use the version without arguments instead.
     */
    void WorldBox(double* xlo, double* xhi,
                  double* ylo, double* yhi,
                  double* zlo, double* zhi) const;

    /// Returns a box with the extremes of the world volume (from shape axes).
    /// @see `GetWorldVolumeName()`
    geo::BoxBoundedGeo WorldBox() const;

    /**
     * @brief The position of the detector respect to earth surface
     * @return typical y position at surface in units of cm
     *
     * This is the depth (y) of the surface (where earth meets air) for this
     * detector site.
     * The number is expressed in world coordinates and in centimetres,
     * and it represents the y coordinate of earth surface.
     * A negative value means that the origin of coordinates, typically matching
     * the detector centre, is above surface.
     *
     * @todo check that this is actually how it is used
     */
    //
    geo::Length_t SurfaceY() const { return fSurfaceY; }


    //
    // object description and information
    //

    /// Access to the ROOT geometry description manager
    TGeoManager* ROOTGeoManager() const;

    /// Return the name of the world volume (needed by Geant4 simulation)
    const std::string GetWorldVolumeName() const;

    /// Returns the absolute  coordinates of the detector enclosure volume [cm].
    /// @param name name of the volume to be sought (default: `volDetEnclosure`)
    /// @throw cet::exception if the specified volume is not found
    geo::BoxBoundedGeo DetectorEnclosureBox
      (std::string const& name = "volDetEnclosure") const;


    //@{
    /**
     * @brief Returns the name of the deepest volume containing specified point
     * @param point the location to query, in world coordinates
     * @return name of the volume containing the point
     *
     * @todo what happens if none?
     * @todo Unify the coordinates type
     */
    std::string VolumeName(geo::Point_t const& point) const;
    std::string VolumeName(TVector3 const& point) const
      { return VolumeName(geo::vect::toPoint(point)); }
    //@}


    /**
     * @brief Returns all the nodes with volumes with any of the specified names
     * @param vol_names list of names of volumes
     * @return list of nodes found
     *
     * All the nodes in the geometry are checked, and all the ones that contain
     * a volume with a name among the ones specified in vol_names are saved
     * in the collection and returned.
     */
    std::vector<TGeoNode const*> FindAllVolumes
      (std::set<std::string> const& vol_names) const;

    /**
     * @brief Returns paths of all nodes with volumes with the specified names
     * @param vol_names list of names of volumes
     * @return list paths of the found nodes
     *
     * All the nodes in the geometry are checked, and the path of all the ones
     * that contain a volume with a name among the ones specified in vol_names
     * is saved in the collection and returned.
     * A node path is a ordered list of all nodes leading to the final one,
     * starting from thetop level (root) down. The node at the `back()` of the
     * path is the one with name in vol_names.
     * No empty paths are returned.
     */
    std::vector<std::vector<TGeoNode const*>> FindAllVolumePaths
      (std::set<std::string> const& vol_names) const;


    /// Returns the material at the specified position
    TGeoMaterial const* Material(geo::Point_t const& point) const;
    //@{
    /**
     * @brief Name of the deepest material containing the point xyz
     * @return material of the origin by default
     */
    std::string MaterialName(TVector3 const& point) const
      { return MaterialName(geo::vect::toPoint(point)); }
    std::string MaterialName(geo::Point_t const& point) const;
    //@}


    //@{
    /// Returns the total mass [kg] of the specified volume (default: world).
    double TotalMass() const { return TotalMass(GetWorldVolumeName()); }
    double TotalMass(std::string vol) const;
    //@}

    //@{
    /**
     * @brief Returns the column density between two points.
     * @param p1 the first point
     * @param p2 the second point
     * @return the column density [kg / cm&sup2;]
     *
     * The column density is defined as
     * @f$ \int_{\vec{p}_{1}}^{\vec{p}_{2}} \rho(\vec{p}) d\vec{p} @f$
     * where @f$ \rho(\vec{p}) @f$ is the density at point @f$ \vec{p} @f$,
     * which the integral leads from `p1` to `p2` in a straight line.
     *
     * Both points are specified in world coordinates.
     */
    double MassBetweenPoints
      (geo::Point_t const& p1, geo::Point_t const& p2) const;
    double MassBetweenPoints(double *p1, double *p2) const;
    //@}


    /// Prints geometry information with maximum verbosity.
    template <typename Stream>
    void Print(Stream&& out, std::string indent = "  ") const;

    /// @brief Returns a string with complete geometry information.
    /// @see `Print()`
    std::string Info(std::string indent = "  ") const;

    /// @}
    // END Detector information


    /**
     * @brief Returns the ID of the first element of the detector.
     * @tparam GeoID type of the ID to be returned
     * @return ID of the first subelement in the detector
     */
    template <typename GeoID>
    GeoID GetBeginID() const { GeoID id; GetBeginID(id); return id; }

    /**
     * @brief Returns the ID next to the specified one.
     * @tparam GeoID type of the ID to be returned
     * @param id the element ID to be incremented
     * @return ID of the next subelement after `id`
     */
    template <typename GeoID>
    GeoID GetNextID(GeoID const& id) const
      { auto nextID(id); IncrementID(nextID); return nextID; }

    /**
     * @brief Returns the (possibly invalid) ID after the last subelement of
     *        the detector.
     * @tparam GeoID type of the ID to be returned
     * @return ID after the last subelement in the specified geometry element
     */
    template <typename GeoID>
    GeoID GetEndID() const { GeoID id; GetEndID(id); return id; }


    /**
     * @brief Returns the ID of the first subelement of the specified element.
     * @tparam GeoID type of the ID to be returned
     * @tparam ContextID type of the ID of the containing element
     * @param id ID of the containing element
     * @return ID of the first subelement in the specified geometry element
     */
    template <typename GeoID, typename ContextID>
    GeoID GetBeginID(ContextID const& id) const;

    /**
     * @brief Returns the (possibly invalid) ID after the last subelement of
     *        the specified element.
     * @tparam GeoID type of the ID to be returned
     * @tparam ContextID type of the ID of the containing element
     * @param id ID of the containing element
     * @return ID (possibly invalid) after the last subelement in the
     *         specified geometry element
     */
    template <typename GeoID, typename ContextID>
    GeoID GetEndID(ContextID const& id) const;


    /// @name Cryostat access and information
    /// @{

    //
    // group features
    //

    //@{
    /**
     * @brief Returns the number of cryostats in the detector
     *
     * The NElements() and NSiblingElements() methods are overloaded and their
     * return depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Ncryostats() const { return Cryostats().size(); }
    unsigned int NElements() const { return Ncryostats(); }
    unsigned int NSiblingElements(geo::CryostatID const&) const
      { return Ncryostats(); }
    //@}

    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified cryostat
     *
     * The HasElement() method is overloaded and its meaning depends on the type
     * of ID.
     */
    bool HasCryostat(geo::CryostatID const& cryoid) const
      { return cryoid.Cryostat < Ncryostats(); }
    bool HasElement(geo::CryostatID const& cryoid) const
      { return HasCryostat(cryoid); }
    //@}

    //@{
    /**
     * @brief Returns the specified cryostat
     * @param cstat number of cryostat
     * @param cryoid cryostat ID
     * @return a constant reference to the specified cryostat
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     *
     * The GetElement() method is overloaded and its return depends on the type
     * of ID.
     *
     * @todo Make the cryostat number mandatory (as CryostatID)
     */
    CryostatGeo const& Cryostat(geo::CryostatID const& cryoid) const;
    CryostatGeo const& Cryostat(unsigned int const cstat = 0) const
      { return Cryostat(geo::CryostatID(cstat)); }
    CryostatGeo const& GetElement(geo::CryostatID const& cryoid) const
      { return Cryostat(cryoid); }
    //@}

    //@{
    /**
     * @brief Returns the specified cryostat
     * @param cryoid cryostat ID
     * @return a constant pointer to the specified cryostat, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the
     * type of ID.
     */
    CryostatGeo const* CryostatPtr(geo::CryostatID const& cryoid) const
      { return HasCryostat(cryoid)? &(Cryostats()[cryoid.Cryostat]): nullptr; }
    CryostatGeo const* GetElementPtr(geo::CryostatID const& cryoid) const
      { return CryostatPtr(cryoid); }
    //@}

    //@{
    /**
     * @brief Returns the index of the cryostat at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @return the index of the cryostat, or UINT_MAX if no cryostat is there
     *
     * @deprecated Use `PositionToCryostatID()` instead
     */
    geo::CryostatID::CryostatID_t FindCryostatAtPosition
      (geo::Point_t const& worldLoc) const;
    geo::CryostatID::CryostatID_t FindCryostatAtPosition
      (double const worldLoc[3]) const;
    //@}


    /**
     * @brief Returns the cryostat at specified location.
     * @param point the location [cm]
     * @return pointer to the `geo::CryostatGeo` including `point`, or `nullptr`
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    geo::CryostatGeo const* PositionToCryostatPtr
      (geo::Point_t const& point) const;

    /**
     * @brief Returns the ID of the cryostat at specified location.
     * @param point the location [cm]
     * @return ID of the cryostat including `point` (invalid if none)
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    geo::CryostatID PositionToCryostatID(geo::Point_t const& point) const;


    //@{
    /**
     * @brief Returns the cryostat at specified location.
     * @param point the location [cm]
     * @return a constant reference to the `geo::CryostatGeo` containing `point`
     * @throws cet::exception ("Geometry" category) if no cryostat matches
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    CryostatGeo const& PositionToCryostat(geo::Point_t const& point) const;
    CryostatGeo const& PositionToCryostat(double const point[3]) const
      { return PositionToCryostat(geo::vect::makePointFromCoords(point)); }
    //@}

    /**
     * @brief Returns the cryostat at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param cid (output) cryostat ID
     * @return a constant reference to the CryostatGeo object of the cryostat
     * @throws cet::exception ("Geometry" category) if no cryostat matches
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     *
     * @deprecated Use `PositionToCryostat(geo::Point_t const&)` instead.
     */
    CryostatGeo const& PositionToCryostat
      (double const worldLoc[3], geo::CryostatID& cid) const;

    /**
     * @brief Returns the cryostat at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param cstat (output) number of cryostat
     * @return a constant reference to the CryostatGeo object of the cryostat
     * @throws cet::exception ("Geometry" category) if no cryostat matches
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     *
     * @deprecated Use `PositionToCryostat(geo::Point_t const&)` instead.
     */
    CryostatGeo const& PositionToCryostat
      (double const worldLoc[3], unsigned int &cstat) const;

    //
    // iterators
    //

    /// Initializes the specified ID with the ID of the first cryostat
    void GetBeginID(geo::CryostatID& id) const
      { id = geo::CryostatID(0, HasCryostat(geo::CryostatID(0))); }

    /// Initializes the specified ID with the invalid ID after the last cryostat
    void GetEndID(geo::CryostatID& id) const
      { id = geo::CryostatID(Ncryostats(), false); }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(geo::CryostatID& id) const; // inline implementation

    /// Returns an iterator pointing to the first cryostat ID
    cryostat_id_iterator begin_cryostat_id() const
      { return cryostat_id_iterator(this, cryostat_id_iterator::begin_pos); }

    /// Returns an iterator pointing after the last cryostat ID
    cryostat_id_iterator end_cryostat_id() const
      { return cryostat_id_iterator(this, cryostat_id_iterator::end_pos); }

    /// Returns an iterator pointing to the first cryostat
    cryostat_iterator begin_cryostat() const
      { return cryostat_iterator(this, cryostat_iterator::begin_pos); }

    /// Returns an iterator pointing after the last cryostat
    cryostat_iterator end_cryostat() const
      { return cryostat_iterator(this, cryostat_iterator::end_pos); }

    /**
     * @brief Enables ranged-for loops on all cryostat IDs of the detector
     * @returns an object suitable for ranged-for loops on all cryostat IDs
     *
     * Example of usage:
     *
     *     for (geo::CryostatID const& cID: geom->IterateCryostatIDs()) {
     *       geo::CryostatGeo const& Cryo = geom->Cryostat(cID);
     *
     *       // useful code here
     *
     *     } // for all cryostats
     *
     */
    IteratorBox<
      cryostat_id_iterator,
      &GeometryCore::begin_cryostat_id, &GeometryCore::end_cryostat_id
      >
    IterateCryostatIDs() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all cryostats of the detector
     * @returns an object suitable for ranged-for loops on all cryostats
     *
     * Example of usage:
     *
     *     for (geo::CryostatGeo const& Cryo: geom->IterateCryostats()) {
     *
     *       // useful code here
     *
     *     } // for all cryostats
     *
     */
    IteratorBox<
      cryostat_iterator,
      &GeometryCore::begin_cryostat, &GeometryCore::end_cryostat
      >
    IterateCryostats() const { return { this }; }

    //
    // single object features
    //

    //@{
    /// Returns the half width of the cryostat (x direction)
    geo::Length_t CryostatHalfWidth(geo::CryostatID const& cid) const;
    geo::Length_t CryostatHalfWidth(unsigned int cstat = 0) const
      { return CryostatHalfWidth(geo::CryostatID(cstat)); }
    //@}

    //@{
    /// Returns the height of the cryostat (y direction)
    geo::Length_t CryostatHalfHeight(geo::CryostatID const& cid) const;
    geo::Length_t CryostatHalfHeight(unsigned int cstat = 0) const
      { return CryostatHalfHeight(geo::CryostatID(cstat)); }
    //@}

    //@{
    /// Returns the length of the cryostat (z direction)
    geo::Length_t CryostatLength(geo::CryostatID const& cid) const;
    geo::Length_t CryostatLength(unsigned int cstat = 0) const
      { return CryostatLength(geo::CryostatID(cstat)); }
    //@}


    /**
     * @brief Returns the boundaries of the specified cryostat
     * @param boundaries (output) pointer to an area of 6 doubles for boundaries
     * @param cid cryostat ID
     * @throws cet::exception ("GeometryCore" category) if cryostat not present
     * @see CryostatGeo::Boundaries()
     *
     * The boundaries array is filled with:
     * [0] lower x coordinate  [1] upper x coordinate
     * [2] lower y coordinate  [3] upper y coordinate
     * [4] lower z coordinate  [5] upper z coordinate
     *
     * @deprecated Use `CryostatGeo::Boundaries()` (from `Cryostat(cid)`).
     * @todo What happen on invalid cryostat?
     */
    void CryostatBoundaries
      (double* boundaries, geo::CryostatID const& cid) const;

    /**
     * @brief Returns the boundaries of the specified cryostat
     * @param boundaries (output) pointer to an area of 6 doubles for boundaries
     * @param cstat number of cryostat
     * @throws cet::exception ("GeometryCore" category) if cryostat not present
     * @see CryostatGeo::Boundaries()
     *
     * The boundaries array is filled with:
     * [0] lower x coordinate  [1] upper x coordinate
     * [2] lower y coordinate  [3] upper y coordinate
     * [4] lower z coordinate  [5] upper z coordinate
     *
     * @deprecated Use `CryostatBoundaries(double*, geo::CryostatID const&)`
     * or (recommended) `CryostatGeo::Boundaries()` from `Cryostat(cid)` instead
     */
    void CryostatBoundaries
      (double* boundaries, unsigned int cstat = 0) const
      { CryostatBoundaries(boundaries, geo::CryostatID(cstat)); }


    //
    // object description
    //

    //@{
    /**
     * @brief Return the name of LAr TPC volume
     * @param cstat index of the cryostat
     * @return the name of the specified TPC
     *
     * This information is used in the event display.
     *
     * @todo Use a cryostat ID instead
     * @todo What if it does not exist?
     */
    std::string GetCryostatVolumeName(geo::CryostatID const& cid) const;
    std::string GetCryostatVolumeName(unsigned int const cstat = 0) const
      { return GetCryostatVolumeName(geo::CryostatID(cstat)); }
    //@}

    /// @} Cryostat access and information



    /// @name TPC access and information
    /// @{

    //
    // group features
    //

    /**
     * @brief Returns the total number of TPCs in the specified cryostat
     * @param cstat cryostat number
     *
     * @todo Make the cryostat number mandatory (as CryostatID)
     * @todo Change return type to size_t
     * @todo what happens if it does not exist?
     */
    unsigned int NTPC(unsigned int cstat = 0) const
      { return NTPC(geo::CryostatID(cstat)); }

    /// Returns the largest number of TPCs a cryostat in the detector has
    unsigned int MaxTPCs() const;

    /// Returns the total number of TPCs in the detector
    unsigned int TotalNTPC() const;


    /**
     * @brief Returns a container with one entry per TPC.
     * @tparam T type of data in the container
     * @return a container with one default-constructed `T` per TPC
     * @see `geo::TPCDataContainer`
     *
     * The working assumption is that all cryostats have the same number of
     * TPCs. It is always guaranteed that all existing TPCs have an entry in
     * the container, although if the previous working assumption is not
     * satisfied there will be entries in the containers which are not
     * associated to a valid TPC.
     *
     * The interface of the container is detailed in the documentation of the
     * container itself, `geo::TPCDataContainer`. Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto tracksPerTPC
     *   = geom->makeTPCData<std::vector<recob::Track const*>>();
     *
     * for (recob::Track const& track: tracks) {
     *   geo::TPCGeo const* tpc = geom->PositionToTPCptr(track.Start());
     *   if (tpc) tracksPerTPC[tpc->ID()].push_back(&track);
     * } // for
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * where the container will be filled with pointers to all tracks starting
     * from a given TPC (tracks reconstructed as starting outside the TPCs will
     * be not saved in the container).
     */
    template <typename T>
    geo::TPCDataContainer<T> makeTPCData() const
      { return { Ncryostats(), MaxTPCs() }; }

    /**
     * @brief Returns a container with one entry per TPC.
     * @tparam T type of data in the container
     * @param defValue the initial value of all elements in the container
     * @return a container with a value `defValue` per each TPC
     * @see `geo::TPCDataContainer`
     *
     * This function operates as `makeTPCData() const`, except that copies
     * the specified value into all the entries of the container. Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto nTracksPerTPC = geom->makeTPCData(0U);
     *
     * for (recob::Track const& track: tracks) {
     *   geo::TPCGeo const* tpc = geom->PositionToTPCptr(track.Start());
     *   if (tpc) ++(tracksPerTPC[tpc->ID()]);
     * } // for
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    geo::TPCDataContainer<T> makeTPCData(T const& defValue) const
      { return { Ncryostats(), MaxTPCs(), defValue }; }



    //@{
    /**
     * @brief Returns the total number of TPCs in the specified cryostat
     * @param cryoid cryostat number
     * @return number of TPCs in specified cryostat, or 0 if no cryostat found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their
     * return depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int NTPC(geo::CryostatID const& cryoid) const
      {
        CryostatGeo const* pCryo = GetElementPtr(cryoid);
        return pCryo? pCryo->NElements(): 0;
      }
    unsigned int NElements(geo::CryostatID const& cryoid) const
      { return NTPC(cryoid); }
    unsigned int NSiblingElements(geo::TPCID const& tpcid) const
      { return NTPC(tpcid); }
    //@}


    //
    // access
    //
    /// Returns whether we have the specified TPC
    bool HasTPC(geo::TPCID const& tpcid) const
      {
        CryostatGeo const* pCryo = CryostatPtr(tpcid);
        return pCryo? pCryo->HasTPC(tpcid): false;
      }

    /// Returns whether we have the specified TPC
    bool HasElement(geo::TPCID const& tpcid) const { return HasTPC(tpcid); }


    //@{
    /**
     * @brief Returns the specified TPC
     * @param tpcid ID of the tpc
     * @param tpc tpc number within the cryostat
     * @param cstat number of cryostat
     * @return a constant reference to the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     *
     * The GetElement() method is overloaded and its return depends on the type
     * of ID.
     *
     * @todo remove the version with integers
     */
    TPCGeo const& TPC
      (unsigned int const tpc   = 0, unsigned int const cstat = 0) const
      { return TPC(geo::TPCID(cstat, tpc)); }
    TPCGeo const& TPC(geo::TPCID const& tpcid) const
      { return Cryostat(tpcid).TPC(tpcid); }
    TPCGeo const& GetElement(geo::TPCID const& tpcid) const
      { return TPC(tpcid); }
    //@}

    //@{
    /**
     * @brief Returns the specified TPC
     * @param tpcid TPC ID
     * @return a constant pointer to the specified TPC, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the
     * type of ID.
     */
    TPCGeo const* TPCPtr(geo::TPCID const& tpcid) const
      {
        CryostatGeo const* pCryo = CryostatPtr(tpcid);
        return pCryo? pCryo->TPCPtr(tpcid): nullptr;
      } // TPCPtr()
    TPCGeo const* GetElementPtr(geo::TPCID const& tpcid) const
      { return TPCPtr(tpcid); }
    //@}

    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param worldLoc 3D coordinates of the point (world reference frame) [cm]
     * @return the TPC ID, or an invalid one if no TPC is there
     */
    geo::TPCID FindTPCAtPosition(double const worldLoc[3]) const
      { return FindTPCAtPosition(geo::vect::makePointFromCoords(worldLoc)); }

    //@{
    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param worldLoc 3D point (world reference frame, centimeters)
     * @return the TPC ID, or an invalid one if no TPC is there
     */
    geo::TPCID FindTPCAtPosition(geo::Point_t const& point) const;
    geo::TPCID FindTPCAtPosition(TVector3 const& point) const
      { return FindTPCAtPosition(geo::vect::toPoint(point)); }
    //@}

    /**
     * @brief Returns the TPC at specified location.
     * @param point the location [cm]
     * @return the `geo::TPCGeo` including `point`, or `nullptr` if none
     */
    geo::TPCGeo const* PositionToTPCptr(geo::Point_t const& point) const;


    //@{
    /**
     * @brief Returns the TPC at specified location.
     * @param point the location [cm]
     * @return a constant reference to the `geo::TPCGeo` including `point`
     * @throws cet::exception ("Geometry" category) if no TPC matches
     */
    geo::TPCGeo const& PositionToTPC(geo::Point_t const& point) const;
    TPCGeo const& PositionToTPC(double const point[3]) const
      { return PositionToTPC(geo::vect::makePointFromCoords(point)); }
    //@}

    /**
     * @brief Returns the TPC at specified location.
     * @param point the location [cm]
     * @param tpc _(output)_ where to store the number of TPC
     * @param cstat _(output)_ where to store the number of cryostat
     * @return a constant reference to the `geo::TPCGeo` including `point`
     * @throws cet::exception ("Geometry" category) if no TPC matches
     * @deprecated Use `PositionToTPCID()` or `PositionToTPC().ID()`
     */
    TPCGeo const& PositionToTPC
      (double const worldLoc[3], unsigned int &tpc, unsigned int &cstat) const;

    /**
     * @brief Returns the TPC at specified location.
     * @param point the location [cm]
     * @param tpcid _(output)_ where to store the TPC ID
     * @return a constant reference to the `geo::TPCGeo` including `point`
     * @throws cet::exception ("Geometry" category) if no TPC matches
     * @deprecated Use `PositionToTPCID()` or `PositionToTPC().ID()`
     */
    TPCGeo const& PositionToTPC
      (double const worldLoc[3], TPCID& tpcid) const;

    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param point the location [cm]
     * @return ID of the TPC at specified location, invalid if none
     * @see `PositionToTPC()`
     */
    geo::TPCID PositionToTPCID(geo::Point_t const& point) const;

    ///
    /// iterators
    ///

    /// Initializes the specified ID with the ID of the first TPC.
    void GetBeginID(geo::TPCID& id) const
      { GetBeginID(id.asCryostatID()); id.TPC = 0; }

    /// Initializes the specified ID with the invalid ID after the last TPC.
    void GetEndID(geo::TPCID& id) const
      { GetEndID(id.asCryostatID()); id.TPC = 0; }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(geo::TPCID& id) const; // inline implementation

    /// Returns the ID of the first TPC in the specified cryostat.
    geo::TPCID GetBeginTPCID(geo::CryostatID const& id) const
      { return { id, 0 }; }

    /// Returns the (possibly invalid) ID after the last TPC of the specified
    /// cryostat.
    geo::TPCID GetEndTPCID(geo::CryostatID const& id) const
      { return { id.Cryostat + 1, 0 }; }


    /// Returns an iterator pointing to the first TPC ID in the detector.
    TPC_id_iterator begin_TPC_id() const
      { return TPC_id_iterator(this, TPC_id_iterator::begin_pos); }

    /// Returns an iterator pointing after the last TPC ID in the detector.
    TPC_id_iterator end_TPC_id() const
      { return TPC_id_iterator(this, TPC_id_iterator::end_pos); }

    /// Returns an iterator pointing to the first TPC ID in the specified
    /// cryostat.
    TPC_id_iterator begin_TPC_id(geo::CryostatID const& cid) const
      { return TPC_id_iterator(this, GetBeginTPCID(cid)); }

    /// Returns an iterator pointing after the last TPC ID in the specified
    /// cryostat.
    TPC_id_iterator end_TPC_id(geo::CryostatID const& cid) const
      { return TPC_id_iterator(this, GetEndTPCID(cid)); }

    /// Returns an iterator pointing to the first TPC in the detector
    TPC_iterator begin_TPC() const
      { return TPC_iterator(this, TPC_iterator::begin_pos); }

    /// Returns an iterator pointing after the last TPC in the detector
    TPC_iterator end_TPC() const
      { return TPC_iterator(this, TPC_iterator::end_pos); }

    /// Returns an iterator pointing to the first TPC in the detector
    TPC_iterator begin_TPC(geo::CryostatID const& cid) const
      { return TPC_iterator(this, GetBeginTPCID(cid)); }

    /// Returns an iterator pointing after the last TPC in the detector
    TPC_iterator end_TPC(geo::CryostatID const& cid) const
      { return TPC_iterator(this, GetEndTPCID(cid)); }

    /**
     * @brief Enables ranged-for loops on all TPC IDs of the detector
     * @returns an object suitable for ranged-for loops on all TPC IDs
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (geo::TPCID const& tID: geom->IterateTPCIDs()) {
     *   geo::TPCGeo const& TPC = geom->TPC(tID);
     *
     *   // useful code here
     *
     * } // for all TPC
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    IteratorBox<
      TPC_id_iterator,
      &GeometryCore::begin_TPC_id, &GeometryCore::end_TPC_id
      >
    IterateTPCIDs() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all TPC IDs of the specified cryostat.
     * @param cid the ID of the cryostat to loop the TPC IDs of
     * @returns an object suitable for ranged-for loops on TPC IDs
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::CryostatID cid{1}; // cryostat #1 (hope it exists!)
     * for (geo::TPCID const& tID: geom->IterateTPCIDs(cid)) {
     *   geo::TPCGeo const& TPC = geom->TPC(tID);
     *
     *   // useful code here
     *
     * } // for all TPC in cryostat #1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      TPC_id_iterator, geo::CryostatID,
      &GeometryCore::begin_TPC_id, &GeometryCore::end_TPC_id
      >
    IterateTPCIDs(geo::CryostatID const& cid) const { return { this, cid }; }

    /// `IterateTPCIDs()` is not supported on TPC IDs.
    void IterateTPCIDs(geo::TPCID const& pid) const = delete;

    /// `IterateTPCIDs()` is not supported on plane IDs.
    void IterateTPCIDs(geo::PlaneID const& pid) const = delete;

    /// `IterateTPCIDs()` is not supported on wire IDs.
    void IterateTPCIDs(geo::WireID const& pid) const = delete;

    /// `IterateTPCIDs()` is not supported on readout IDs.
    void IterateTPCIDs(readout::TPCsetID const&) const = delete;

    /// `IterateTPCIDs()` is not supported on readout IDs.
    void IterateTPCIDs(readout::ROPID const&) const = delete;

    /**
     * @brief Enables ranged-for loops on all TPCs of the detector.
     * @returns an object suitable for ranged-for loops on all TPCs
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (geo::TPCGeo const& TPC: geom->IterateTPCs()) {
     *
     *   // useful code here
     *
     * } // for TPCs
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    IteratorBox
      <TPC_iterator, &GeometryCore::begin_TPC, &GeometryCore::end_TPC>
    IterateTPCs() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all TPCs of the specified cryostat.
     * @param cid the ID of the cryostat to loop the TPCs of
     * @returns an object suitable for ranged-for loops on TPCs
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::CryostatID cid{1}; // cryostat #1 (hope it exists!)
     * for (geo::TPCGeo const& TPC: geom->IterateTPCs(cid)) {
     *
     *   // useful code here
     *
     * } // for TPCs in cryostat 1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      TPC_iterator, geo::CryostatID,
      &GeometryCore::begin_TPC, &GeometryCore::end_TPC
      >
    IterateTPCs(geo::CryostatID const& cid) const { return { this, cid }; }

    /// `IterateTPCs()` is not supported on TPC IDs.
    void IterateTPCs(geo::TPCID const& pid) const = delete;

    /// `IterateTPCs()` is not supported on plane IDs.
    void IterateTPCs(geo::PlaneID const& pid) const = delete;

    /// `IterateTPCs()` is not supported on wire IDs.
    void IterateTPCs(geo::WireID const& pid) const = delete;

    /// `IterateTPCs()` is not supported on readout IDs.
    void IterateTPCs(readout::TPCsetID const&) const = delete;

    /// `IterateTPCs()` is not supported on readout IDs.
    void IterateTPCs(readout::ROPID const&) const = delete;


    //
    // single object features
    //

    //@{
    /**
     * @brief Returns the half width of the active volume of the specified TPC.
     * @param tpcid ID of the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return the value of the half width of the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @see geo::TPCGeo::ActiveHalfWidth()
     *
     * @todo deprecate this function
     * @todo rename the function
     */
    geo::Length_t DetHalfWidth(geo::TPCID const& tpcid) const;
    geo::Length_t DetHalfWidth
      (unsigned int tpc = 0, unsigned int cstat = 0) const
      { return DetHalfWidth(geo::TPCID(cstat, tpc)); }
    //@}

    //@{
    /**
     * @brief Returns the half height of the active volume of the specified TPC.
     * @param tpcid ID of the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return the value of the half height of the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @see geo::TPCGeo::ActiveHalfHeight()
     *
     * See `geo::TPCGeo::ActiveHalfHeight()` for more details.
     *
     * @todo deprecate this function
     * @todo rename the function
     */
    geo::Length_t DetHalfHeight(geo::TPCID const& tpcid) const;
    geo::Length_t DetHalfHeight
      (unsigned int tpc = 0, unsigned int cstat = 0) const
      { return DetHalfHeight(geo::TPCID(cstat, tpc)); }
    //@}

    //@{
    /**
     * @brief Returns the length of the active volume of the specified TPC.
     * @param tpcid ID of the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return the value of the length of the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @see geo::TPCGeo::ActiveLength()
     *
     * See `geo::TPCGeo::ActiveLength()` for more details.
     *
     * @todo deprecate this function
     * @todo rename the function
     */
    geo::Length_t DetLength(geo::TPCID const& tpcid) const;
    geo::Length_t DetLength(unsigned int tpc = 0, unsigned int cstat = 0) const
      { return DetLength(geo::TPCID(cstat, tpc)); }
    //@}


    //@{
    /**
     * @brief Returns the center of side of the detector facing the beam.
     * @tparam Point _(default: `DefaultPoint_t`)_ return this point type
     * @param tpcid ID of the TPC
     * @return position of center of TPC face toward the beam,
     *         in world coordinates [cm]
     *
     * Effectively, this is the center of the side of TPC active volume
     * which faces the negative _z_ direction, the first that a beam following
     * the positive _z_ direction crosses.
     */
    template <typename Point>
    Point GetTPCFrontFaceCenter(geo::TPCID const& tpcid) const
      { return TPC(tpcid).GetFrontFaceCenter<Point>(); }
    DefaultPoint_t GetTPCFrontFaceCenter(geo::TPCID const& tpcid) const
      { return GetTPCFrontFaceCenter<DefaultPoint_t>(tpcid); }
    //@}
    
    //@{
    /**
     * @brief Returns the center of side of the detector facing the beam.
     * @tparam Point _(default: `DefaultPoint_t`)_ return this point type
     * @param tpc _(default: `0`)_ TPC number within the cryostat `cstat`
     * @param cstat _(default: `0`)_ number of cryostat
     * @return position of center of TPC face toward the beam,
     *         in world coordinates [cm]
     * @see `GetTPCFrontFaceCenter(geo::TPCID const&)`
     *
     * @note Please use `GetTPCFrontFaceCenter(geo::TPCID const&)` instead.
     */
    template <typename Point>
    Point GetTPCFrontFaceCenter
      (unsigned int tpc = 0, unsigned int cstat = 0) const
      { return GetTPCFrontFaceCenter<Point>(geo::TPCID(cstat, tpc)); }
    DefaultPoint_t GetTPCFrontFaceCenter
      (unsigned int tpc = 0, unsigned int cstat = 0) const
      { return GetTPCFrontFaceCenter<DefaultPoint_t>(tpc, cstat); }
    //@}


    //
    // object description
    //

    //@{
    /**
     * @brief Return the name of specified LAr TPC volume
     * @param tpcid ID of the TPC
     * @param tpc index of TPC in the cryostat
     * @param cstat index of the cryostat
     * @return the name of the specified TPC
     *
     * This information is used by Geant4 simulation
     *
     * @todo Use a TPCID instead
     * @todo What if it does not exist?
     */
    std::string GetLArTPCVolumeName(geo::TPCID const& tpcid) const;
    std::string GetLArTPCVolumeName
      (unsigned int const tpc = 0, unsigned int const cstat = 0) const
      { return GetLArTPCVolumeName(geo::TPCID(cstat, tpc)); }
    //@}

    /// @} TPC access and information



    /// @name Plane access and information
    /// @{

    //
    // group features
    //

    /**
     * @brief Returns the total number of wire planes in the specified TPC
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     *
     * @todo Make all the arguments mandatory (as TPCID)
     * @todo Change return type to size_t
     * @todo what happens if TPC does not exist?
     */
    unsigned int Nplanes(unsigned int tpc   = 0, unsigned int cstat = 0) const
      { return Nplanes(geo::TPCID(cstat, tpc)); }

    /// Returns the largest number of planes among all TPCs in this detector
    unsigned int MaxPlanes() const;

    /**
     * @brief Returns a container with one entry per wire plane.
     * @tparam T type of data in the container
     * @return a container with one default-constructed `T` per plane
     * @see `geo::PlaneDataContainer`
     *
     * The working assumption is that all cryostats have the same number of
     * TPCs, and all TPCs have the same number of planes. It is always
     * guaranteed that all existing planes have an entry in the container,
     * although if the previous working assumption is not satisfied there will
     * be entries in the containers which are not associated to a valid plane.
     *
     * The interface of the container is detailed in the documentation of the
     * container itself, `geo::PlaneDataContainer`. Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto hitsPerPlane
     *   = geom->makePlaneData<std::vector<recob::Hit const*>>();
     *
     * for (recob::Hit const& hit: hits) {
     *   if (hit.WireID()) hitsPerPlane[hit.WireID()].push_back(&hit);
     * } // for
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * where the container will be filled with pointers to all hits on the given
     * wire plane (wire IDs are implicitly converted into plane IDs in the index
     * `operator[]` call).
     */
    template <typename T>
    geo::PlaneDataContainer<T> makePlaneData() const
      { return { Ncryostats(), MaxTPCs(), MaxPlanes() }; }

    /**
     * @brief Returns a container with one entry per wire plane.
     * @tparam T type of data in the container
     * @param defValue the initial value of all elements in the container
     * @return a container with one default-constructed `T` per plane
     * @see `geo::PlaneDataContainer`
     *
     * This function operates as `makePlaneData() const`, except that copies
     * the specified value into all the entries of the container. Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto nHitsPerPlane = geom->makePlaneData(0U);
     *
     * for (recob::Hit const& hit: hits) {
     *   if (hit.WireID()) ++(hitsPerPlane[hit.WireID()]);
     * } // for
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    geo::PlaneDataContainer<T> makePlaneData(T const& defValue) const
      { return { Ncryostats(), MaxTPCs(), MaxPlanes(), defValue }; }


    //@{
    /**
     * @brief Returns the total number of planes in the specified TPC
     * @param tpcid TPC ID
     * @return number of planes in specified TPC, or 0 if no TPC found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their
     * return depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nplanes(geo::TPCID const& tpcid) const
      {
        TPCGeo const* pTPC = GetElementPtr(tpcid);
        return pTPC? pTPC->NElements(): 0;
      }
    unsigned int NElements(geo::TPCID const& tpcid) const
      { return Nplanes(tpcid); }
    unsigned int NSiblingElements(geo::PlaneID const& planeid) const
      { return Nplanes(planeid); }
    //@}


    /**
     * @brief Returns the number of views (different wire orientations)
     *
     * Returns the number of different views, or wire orientations, in the
     * detector.
     *
     * The function assumes that all TPCs in all cryostats of a detector have
     * the same number of planes, which should be a safe assumption.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nviews() const;

    /**
     * @brief Returns a list of possible PlaneIDs in the detector
     * @return a constant reference to the set of plane IDs
     *
     * @deprecated use IteratePlaneIDs() instead
     * plane IDs of DUNE FD? probably better to use iterators instead
     */
    [[deprecated("Iterate through geo::GeometryCore::IteratePlaneIDs() instead")]]
    std::set<PlaneID> const& PlaneIDs() const;


    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified plane
     *
     * The HasElement() method is overloaded and its meaning depends on the type
     * of ID.
     *
     */
    bool HasPlane(geo::PlaneID const& planeid) const
      {
        geo::TPCGeo const* pTPC = TPCPtr(planeid);
        return pTPC? pTPC->HasPlane(planeid): false;
      }
    bool HasElement(geo::PlaneID const& planeid) const
      { return HasPlane(planeid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param planeid ID of the plane
     * @param p plane number within the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return a constant reference to the specified plane
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @throw cet::exception (`PlaneOutOfRange` category) if no such plane
     *
     * The GetElement() method is overloaded and its return depends on the type
     * of ID.
     *
     * @todo remove the version with integers
     */
    PlaneGeo const& Plane
      (unsigned int const p, unsigned int const tpc   = 0, unsigned int const cstat = 0)
      const
      { return Plane(geo::PlaneID(cstat, tpc, p)); }
    PlaneGeo const& Plane(geo::PlaneID const& planeid) const
      { return TPC(planeid).Plane(planeid); }
    PlaneGeo const& GetElement(geo::PlaneID const& planeid) const
      { return Plane(planeid); }
    //@}

    //@{
    /**
     * @brief Returns the specified plane
     * @param planeid plane ID
     * @return a constant pointer to the specified plane, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the
     * type of ID.
     */
    PlaneGeo const* PlanePtr(geo::PlaneID const& planeid) const
      {
        geo::TPCGeo const* pTPC = TPCPtr(planeid);
        return pTPC? pTPC->PlanePtr(planeid): nullptr;
      } // PlanePtr()
    PlaneGeo const* GetElementPtr(geo::PlaneID const& planeid) const
      { return PlanePtr(planeid); }
    //@}

    //
    // iterators
    //

    /// Initializes the specified ID with the ID of the first plane.
    void GetBeginID(geo::PlaneID& id) const
      { GetBeginID(id.asTPCID()); id.Plane = 0; }

    /// Initializes the specified ID with the invalid ID after the last plane.
    void GetEndID(geo::PlaneID& id) const
      { GetEndID(id.asTPCID()); id.Plane = 0; }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(geo::PlaneID& id) const; // inline implementation

    /// Returns the ID of the first plane of the specified cryostat.
    geo::PlaneID GetBeginPlaneID(geo::CryostatID const& id) const
      { return { GetBeginTPCID(id), 0 }; }

    /// Returns the (possibly invalid) ID after the last plane of the specified
    /// cryostat.
    geo::PlaneID GetEndPlaneID(geo::CryostatID const& id) const
      { return { GetEndTPCID(id), 0 }; }

    /// Returns the ID of the first plane of the specified TPC.
    geo::PlaneID GetBeginPlaneID(geo::TPCID const& id) const
      { return { id, 0 }; }

    /// Returns the (possibly invalid) ID after the last plane of the specified
    /// TPC.
    geo::PlaneID GetEndPlaneID(geo::TPCID const& id) const
      { return { GetNextID(id), 0 }; }

    /// Returns an iterator pointing to the first plane ID in the detector
    plane_id_iterator begin_plane_id() const
      { return plane_id_iterator(this, plane_id_iterator::begin_pos); }

    /// Returns an iterator pointing after the last plane ID in the detector
    plane_id_iterator end_plane_id() const
      { return plane_id_iterator(this, plane_id_iterator::end_pos); }

    /// Returns an iterator pointing to the first plane ID in the specified
    /// cryostat.
    plane_id_iterator begin_plane_id(geo::CryostatID const& ID) const
      { return plane_id_iterator(this, GetBeginPlaneID(ID)); }

    /// Returns an iterator pointing after the last plane ID in the specified
    /// cryostat.
    plane_id_iterator end_plane_id(geo::CryostatID const& ID) const
      { return plane_id_iterator(this, GetEndPlaneID(ID)); }

    /// Returns an iterator pointing to the first plane ID in the specified
    /// TPC.
    plane_id_iterator begin_plane_id(geo::TPCID const& ID) const
      { return plane_id_iterator(this, GetBeginPlaneID(ID)); }

    /// Returns an iterator pointing after the last plane ID in the specified
    /// TPC.
    plane_id_iterator end_plane_id(geo::TPCID const& ID) const
      { return plane_id_iterator(this, GetEndPlaneID(ID)); }

    /// Returns an iterator pointing to the first plane in the detector
    plane_iterator begin_plane() const
      { return plane_iterator(this, plane_iterator::begin_pos); }

    /// Returns an iterator pointing after the last plane in the detector
    plane_iterator end_plane() const
      { return plane_iterator(this, plane_iterator::end_pos); }

    /// Returns an iterator pointing to the first plane in the specified
    /// cryostat.
    plane_iterator begin_plane(geo::CryostatID const& ID) const
      { return plane_iterator(this, GetBeginPlaneID(ID)); }

    /// Returns an iterator pointing after the last plane in the specified
    /// cryostat.
    plane_iterator end_plane(geo::CryostatID const& ID) const
      { return plane_iterator(this, GetEndPlaneID(ID)); }

    /// Returns an iterator pointing to the first plane in the specified TPC.
    plane_iterator begin_plane(geo::TPCID const& ID) const
      { return plane_iterator(this, GetBeginPlaneID(ID)); }

    /// Returns an iterator pointing after the last plane in the specified TPC.
    plane_iterator end_plane(geo::TPCID const& ID) const
      { return plane_iterator(this, GetEndPlaneID(ID)); }

    /**
     * @brief Enables ranged-for loops on all plane IDs of the detector.
     * @returns an object suitable for ranged-for loops on all plane IDs
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (geo::PlaneID const& pID: geom->IteratePlaneIDs()) {
     *   geo::PlaneGeo const& Plane = geom->Plane(pID);
     *
     *   // useful code here
     *
     * } // for all plane IDs
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    IteratorBox<
      plane_id_iterator,
      &GeometryCore::begin_plane_id, &GeometryCore::end_plane_id
      >
    IteratePlaneIDs() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all plane IDs of the specified
     *        cryostat.
     * @param cid the ID of the cryostat to loop the plane IDs of
     * @returns an object suitable for ranged-for loops on plane IDs
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::CryostatID cid{1}; // cryostat #1 (hope it exists!)
     * for (geo::PlaneID const& pID: geom->IteratePlaneIDs(cid)) {
     *   geo::PlaneGeo const& plane = geom->Plane(pID);
     *
     *   // useful code here
     *
     * } // for all planes in cryostat #1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      plane_id_iterator, geo::CryostatID,
      &GeometryCore::begin_plane_id, &GeometryCore::end_plane_id
      >
    IteratePlaneIDs(geo::CryostatID const& cid) const { return { this, cid }; }

    /**
     * @brief Enables ranged-for loops on all plane IDs of the specified TPC.
     * @param tid the ID of the TPC to loop the plane IDs of
     * @returns an object suitable for ranged-for loops on plane IDs
     *
     * If the TPC ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::TPCID tid{ 0, 1 }; // C:0 T:1 (hope it exists!)
     * for (geo::PlaneID const& pID: geom->IteratePlaneIDs(tid)) {
     *   geo::PlaneGeo const& plane = geom->Plane(pID);
     *
     *   // useful code here
     *
     * } // for all planes in C:0 T:1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      plane_id_iterator, geo::TPCID,
      &GeometryCore::begin_plane_id, &GeometryCore::end_plane_id
      >
    IteratePlaneIDs(geo::TPCID const& tid) const { return { this, tid }; }


    /// `IteratePlaneIDs()` is not supported on plane IDs.
    void IteratePlaneIDs(geo::PlaneID const& pid) const = delete;

    /// `IteratePlaneIDs()` is not supported on wire IDs.
    void IteratePlaneIDs(geo::WireID const& pid) const = delete;

    /// `IteratePlaneIDs()` is not supported on readout IDs.
    void IteratePlaneIDs(readout::TPCsetID const&) const = delete;

    /// `IteratePlaneIDs()` is not supported on readout IDs.
    void IteratePlaneIDs(readout::ROPID const&) const = delete;


    /**
     * @brief Enables ranged-for loops on all planes of the detector.
     * @returns an object suitable for ranged-for loops on all planes
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (geo::PlaneGeo const& Plane: geom->IteratePlanes()) {
     *
     *   // useful code here
     *
     * } // for all planes
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    IteratorBox<
      plane_iterator,
      &GeometryCore::begin_plane, &GeometryCore::end_plane
      >
    IteratePlanes() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all planes of the specified cryostat.
     * @param cid the ID of the cryostat to loop the planes of
     * @returns an object suitable for ranged-for loops on planes
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::CryostatID cid{1}; // cryostat #1 (hope it exists!)
     * for (geo::PlaneGeo const& plane: geom->IteratePlanes(cid)) {
     *
     *   // useful code here
     *
     * } // for planes in cryostat 1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      plane_iterator, geo::CryostatID,
      &GeometryCore::begin_plane, &GeometryCore::end_plane
      >
    IteratePlanes(geo::CryostatID const& cid) const { return { this, cid }; }

    /**
     * @brief Enables ranged-for loops on all planes of the specified TPC.
     * @param tid the ID of the TPC to loop the planes of
     * @returns an object suitable for ranged-for loops on planes
     *
     * If the TPC ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::TPCID tid{ 0, 1 }; // C:0 T:1 (hope it exists!)
     * for (geo::PlaneGeo const& plane: geom->IteratePlanes(tid)) {
     *
     *   // useful code here
     *
     * } // for planes in C:0 T:1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      plane_iterator, geo::TPCID,
      &GeometryCore::begin_plane, &GeometryCore::end_plane
      >
    IteratePlanes(geo::TPCID const& tid) const { return { this, tid }; }

    /// `IteratePlanes()` is not supported on plane IDs.
    void IteratePlanes(geo::PlaneID const& pid) const = delete;

    /// `IteratePlanes()` is not supported on wire IDs.
    void IteratePlanes(geo::WireID const& pid) const = delete;

    /// `IteratePlanes()` is not supported on readout IDs.
    void IteratePlanes(readout::TPCsetID const&) const = delete;

    /// `IteratePlanes()` is not supported on readout IDs.
    void IteratePlanes(readout::ROPID const&) const = delete;


    //
    // single object features
    //

    //@{
    /**
     * @brief Returns the distance between two planes
     * @param p1 index of the first plane
     * @param p2 index of the second plane
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     * @return distance between the planes
     *
     * @todo add a version with plane IDs
     * @todo deprecate this function
     * @todo add a default version for a given TPCID
     * @todo add a version with two plane indices for a given TPCID
     * @todo return the absolute value of the distance (makes the order unimportant)
     * @todo document what will happen (in the future methods) with planes on different TPCs
     */
    geo::Length_t PlanePitch(
      geo::TPCID const& tpcid,
      geo::PlaneID::PlaneID_t p1 = 0, geo::PlaneID::PlaneID_t p2 = 1
      )
      const;
    geo::Length_t PlanePitch(geo::PlaneID const& pid1, geo::PlaneID const& pid2) const;
    geo::Length_t PlanePitch(unsigned int p1 = 0,
                             unsigned int p2 = 1,
                             unsigned int tpc = 0,
                             unsigned int cstat = 0) const;
    //@}

    /**
     * @brief Returns the view (wire orientation) on the channels of specified TPC plane
     * @param plane TPC plane ID
     * @return the type of signal on the specified plane, or geo::kUnknown
     */
    View_t View(geo::PlaneID const& pid) const;

    /**
     * @brief Returns the type of signal on the channels of specified TPC plane
     * @param plane TPC plane ID
     * @return the type of signal on the specified plane, or geo::kMysteryType
     *
     * Assumes that all the channels on the plane have the same signal type.
     *
     * @todo verify that kMysteryType is returned on invalid plane
     */
    SigType_t SignalType(geo::PlaneID const& pid) const;


    /// @} Plane access and information


    /// @name Wire access and information
    /// @{

    //
    // group features
    //

    /**
     * @brief Returns the total number of wires in the specified plane
     * @param p plane number within the TPC
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     *
     * @todo Make all the arguments mandatory (as PlaneID)
     * @todo Change return type to size_t
     * @todo what happens if it does not exist?
     */
    unsigned int Nwires
      (unsigned int p, unsigned int tpc   = 0, unsigned int cstat = 0) const
      { return Nwires(geo::PlaneID(cstat, tpc, p)); }

    //@{
    /**
     * @brief Returns the total number of wires in the specified plane
     * @param planeid plane ID
     * @return number of wires in specified plane, or 0 if no plane found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their
     * return depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nwires(geo::PlaneID const& planeid) const
      {
        PlaneGeo const* pPlane = GetElementPtr(planeid);
        return pPlane? pPlane->NElements(): 0;
      }
    unsigned int NElements(geo::PlaneID const& planeid) const
      { return Nwires(planeid); }
    unsigned int NSiblingElements(geo::WireID const& wireid) const
      { return Nwires(wireid); }

    /// Returns the largest number of wires among all planes in this detector
    unsigned int MaxWires() const;

    //@}


    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified wire
     *
     * The HasElement() method is overloaded and its meaning depends on the type
     * of ID.
     */
    bool HasWire(geo::WireID const& wireid) const
      {
        geo::PlaneGeo const* pPlane = PlanePtr(wireid);
        return pPlane? pPlane->HasWire(wireid): false;
      }
    bool HasElement(geo::WireID const& wireid) const { return HasWire(wireid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param wireid wire ID
     * @return a constant pointer to the specified wire, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the
     * type of ID.
     */
    WireGeo const* WirePtr(geo::WireID const& wireid) const
      {
        geo::PlaneGeo const* pPlane = PlanePtr(wireid);
        return pPlane? pPlane->WirePtr(wireid): nullptr;
      } // WirePtr()
    WireGeo const* GetElementPtr(geo::WireID const& wireid) const
      { return WirePtr(wireid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param wireid ID of the wire
     * @return a constant reference to the specified wire
     * @throw cet::exception if not found
     *
     * The GetElement() method is overloaded and its return depends on the type
     * of ID.
     */
    WireGeo const& Wire(geo::WireID const& wireid) const
      { return Plane(wireid).Wire(wireid); }
    WireGeo const& WireIDToWireGeo(geo::WireID const& wireid) const
      { return Wire(wireid); }
    WireGeo const& GetElement(geo::WireID const& wireid) const
      { return Wire(wireid); }
    //@}

    //
    // iterators
    //

    /// Initializes the specified ID with the ID of the first wire.
    void GetBeginID(geo::WireID& id) const
      { GetBeginID(id.asPlaneID()); id.Wire = 0; }

    /// Initializes the specified ID with the invalid ID after the last wire.
    void GetEndID(geo::WireID& id) const
      { GetEndID(id.asPlaneID()); id.Wire = 0; }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(geo::WireID& id) const; // inline implementation

    /// Returns the ID of the first wire in the specified cryostat.
    geo::WireID GetBeginWireID(geo::CryostatID const& id) const
      { return { GetBeginPlaneID(id), 0 }; }

    /// Returns the (possibly invalid) ID after the last wire in the specified
    /// cryostat.
    geo::WireID GetEndWireID(geo::CryostatID const& id) const
      { return { GetEndPlaneID(id), 0 }; }

    /// Returns the ID of the first wire of the specified TPC.
    geo::WireID GetBeginWireID(geo::TPCID const& id) const
      { return { geo::PlaneID(id, 0), 0 }; }

    /// Returns the (possibly invalid) ID after the last wire of the specified
    /// TPC.
    geo::WireID GetEndWireID(geo::TPCID const& id) const
      { return { geo::PlaneID(GetNextID(id), 0), 0 }; }

    /// Returns the ID of the first wire of the specified wire plane.
    geo::WireID GetBeginWireID(geo::PlaneID const& id) const
      { return { id, 0 }; }

    /// Returns the (possibly invalid) ID after the last wire of the specified
    /// wire plane.
    geo::WireID GetEndWireID(geo::PlaneID const& id) const
      { return { GetNextID(id), 0 }; }

    /// Returns an iterator pointing to the first wire ID in the detector.
    wire_id_iterator begin_wire_id() const
      { return wire_id_iterator(this, wire_id_iterator::begin_pos); }

    /// Returns an iterator pointing after the last wire ID in the detector.
    wire_id_iterator end_wire_id() const
      { return wire_id_iterator(this, wire_id_iterator::end_pos); }

    /// Returns an iterator pointing to the first wire ID in specified cryostat.
    wire_id_iterator begin_wire_id(geo::CryostatID const& id) const
      { return wire_id_iterator(this, GetBeginWireID(id)); }

    /// Returns an iterator pointing after the last wire ID in specified
    /// cryostat.
    wire_id_iterator end_wire_id(geo::CryostatID const& id) const
      { return wire_id_iterator(this, GetEndWireID(id)); }

    /// Returns an iterator pointing to the first wire ID in specified TPC.
    wire_id_iterator begin_wire_id(geo::TPCID const& id) const
      { return wire_id_iterator(this, GetBeginWireID(id)); }

    /// Returns an iterator pointing after the last wire ID in specified TPC.
    wire_id_iterator end_wire_id(geo::TPCID const& id) const
      { return wire_id_iterator(this, GetEndWireID(id)); }

    /// Returns an iterator pointing to the first wire ID in specified plane.
    wire_id_iterator begin_wire_id(geo::PlaneID const& id) const
      { return wire_id_iterator(this, GetBeginWireID(id)); }

    /// Returns an iterator pointing after the last wire ID in specified plane.
    wire_id_iterator end_wire_id(geo::PlaneID const& id) const
      { return wire_id_iterator(this, GetEndWireID(id)); }

    /// Returns an iterator pointing to the first wire in the detector
    wire_iterator begin_wire() const
      { return wire_iterator(this, wire_iterator::begin_pos); }

    /// Returns an iterator pointing after the last wire in the detector
    wire_iterator end_wire() const
      { return wire_iterator(this, wire_iterator::end_pos); }

    /// Returns an iterator pointing to the first wire in specified cryostat.
    wire_iterator begin_wire(geo::CryostatID const& id) const
      { return wire_iterator(begin_wire_id(id)); }

    /// Returns an iterator pointing after the last wire in specified cryostat.
    wire_iterator end_wire(geo::CryostatID const& id) const
      { return wire_iterator(end_wire_id(id)); }

    /// Returns an iterator pointing to the first wire in specified TPC.
    wire_iterator begin_wire(geo::TPCID const& id) const
      { return wire_iterator(begin_wire_id(id)); }

    /// Returns an iterator pointing after the last wire in specified TPC.
    wire_iterator end_wire(geo::TPCID const& id) const
      { return wire_iterator(end_wire_id(id)); }

    /// Returns an iterator pointing to the first wire in specified plane.
    wire_iterator begin_wire(geo::PlaneID const& id) const
      { return wire_iterator(begin_wire_id(id)); }

    /// Returns an iterator pointing after the last wire in specified plane.
    wire_iterator end_wire(geo::PlaneID const& id) const
      { return wire_iterator(end_wire_id(id)); }

    /**
     * @brief Enables ranged-for loops on all wire IDs of the detector.
     * @returns an object suitable for ranged-for loops on all wire IDs
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (geo::WireID const& wID: geom->IterateWireIDs()) {
     *   geo::WireGeo const& Wire = geom->Wire(wID);
     *
     *   // useful code here
     *
     * } // for all wires
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     */
    IteratorBox<
      wire_id_iterator,
      &GeometryCore::begin_wire_id, &GeometryCore::end_wire_id
      >
    IterateWireIDs() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all wire IDs of specified cryostat.
     * @param cid the ID of the cryostat to loop the wires of
     * @returns an object suitable for ranged-for loops on cryostat wire IDs
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::CryostatID cid{1}; // cryostat #1 (hope it exists!)
     * for (geo::WireID const& wID: geom->IterateWireIDs(cid)) {
     *   geo::WireGeo const& Wire = geom->Wire(wID);
     *
     *   // useful code here
     *
     * } // for all wires
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      wire_id_iterator,
      geo::CryostatID,
      &GeometryCore::begin_wire_id, &GeometryCore::end_wire_id
      >
    IterateWireIDs(geo::CryostatID const& cid) const { return { this, cid }; }

    /**
     * @brief Enables ranged-for loops on all wire IDs of specified TPC.
     * @param tid the ID of the TPC to loop the wires of
     * @returns an object suitable for ranged-for loops on TPC wire IDs
     *
     * If the TPC ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::TPCID tid{0, 1}; // C:0 T:1 (hope it exists!)
     * for (geo::WireID const& wID: geom->IterateWireIDs(tid)) {
     *   geo::WireGeo const& Wire = geom->Wire(wID);
     *
     *   // useful code here
     *
     * } // for all wires in C:0 T:1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      wire_id_iterator,
      geo::TPCID,
      &GeometryCore::begin_wire_id, &GeometryCore::end_wire_id
      >
    IterateWireIDs(geo::TPCID const& tid) const { return { this, tid }; }

    /**
     * @brief Enables ranged-for loops on all wire IDs of specified wire plane.
     * @param pid the ID of the wire plane to loop the wires of
     * @returns an object suitable for ranged-for loops on plane wire IDs
     *
     * If the wire plane ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::PlaneID pid{0, 0, 1}; // C:0 T:0 P:1
     * for (geo::WireID const& wID: geom->IterateWireIDs(pid)) {
     *   geo::WireGeo const& Wire = geom->Wire(wID);
     *
     *   // useful code here
     *
     * } // for all wires in C:0 T:0 P:1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      wire_id_iterator,
      geo::PlaneID,
      &GeometryCore::begin_wire_id, &GeometryCore::end_wire_id
      >
    IterateWireIDs(geo::PlaneID const& pid) const { return { this, pid }; }

    /// `IterateWireIDs()` is not supported on wire IDs.
    void IterateWireIDs(geo::WireID const& pid) const = delete;

    /// `IterateWireIDs()` is not supported on readout IDs.
    void IterateWireIDs(readout::TPCsetID const&) const = delete;

    /// `IterateWireIDs()` is not supported on readout IDs.
    void IterateWireIDs(readout::ROPID const&) const = delete;


    /**
     * @brief Enables ranged-for loops on all wires of the detector.
     * @returns an object suitable for ranged-for loops on all wires
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (geo::WireGeo const& Wire: geom->IterateWires()) {
     *
     *   // useful code here
     *
     * } // for all wires
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    IteratorBox<
      wire_iterator,
      &GeometryCore::begin_wire, &GeometryCore::end_wire
      >
    IterateWires() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all wires of specified cryostat.
     * @param cid the ID of the cryostat to loop the wires of
     * @returns an object suitable for ranged-for loops on cryostat wires
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::CryostatID cid{1}; // cryostat #1 (hope it exists!)
     * for (geo::WireID const& Wire: geom->IterateWires(cid)) {
     *
     *   // useful code here
     *
     * } // for all wires
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      wire_iterator,
      geo::CryostatID,
      &GeometryCore::begin_wire, &GeometryCore::end_wire
      >
    IterateWires(geo::CryostatID const& cid) const { return { this, cid }; }

    /**
     * @brief Enables ranged-for loops on all wires of specified TPC.
     * @param tid the ID of the TPC to loop the wires of
     * @returns an object suitable for ranged-for loops on TPC wires
     *
     * If the TPC ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::TPCID tid{0, 1}; // C:0 T:1 (hope it exists!)
     * for (geo::WireID const& Wire: geom->IterateWires(tid)) {
     *
     *   // useful code here
     *
     * } // for all wires in C:0 T:1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      wire_iterator,
      geo::TPCID,
      &GeometryCore::begin_wire, &GeometryCore::end_wire
      >
    IterateWires(geo::TPCID const& tid) const { return { this, tid }; }

    /**
     * @brief Enables ranged-for loops on all wires of specified wire plane.
     * @param pid the ID of the wire plane to loop the wires of
     * @returns an object suitable for ranged-for loops on plane wires
     *
     * If the wire plane ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::PlaneID pid{0, 1}; // C:0 T:0 P:1
     * for (geo::WireID const& Wire: geom->IterateWires(pid)) {
     *
     *   // useful code here
     *
     * } // for all wires in C:0 T:0 T:1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      wire_iterator,
      geo::PlaneID,
      &GeometryCore::begin_wire, &GeometryCore::end_wire
      >
    IterateWires(geo::PlaneID const& tid) const { return { this, tid }; }

    /// `IterateWires()` is not supported on wire IDs.
    void IterateWires(geo::WireID const& pid) const = delete;

    /// `IterateWires()` is not supported on readout IDs.
    void IterateWires(readout::TPCsetID const&) const = delete;

    /// `IterateWires()` is not supported on readout IDs.
    void IterateWires(readout::ROPID const&) const = delete;

    //
    // single object features
    //

    //@{
    /**
     * @brief Returns the distance between two consecutive wires.
     * @param p plane number within the TPC
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     * @return the distance between the two wires
     *
     * @note The current geometry assumptions imply that wire pitch is constant
     *       between all wires on the same wire plane. This is an assumption
     *       non-trivial to remove.
     *
     * @todo add a version with wire IDs
     * @todo deprecate this function
     * @todo document what will happen (in the future methods) with wires on different planes
     *
     */
    geo::Length_t WirePitch(geo::PlaneID const& planeid) const;
    geo::Length_t WirePitch(unsigned int plane = 0,
                            unsigned int tpc = 0,
                            unsigned int cstat = 0) const
      { return WirePitch(geo::PlaneID(cstat, tpc, plane)); }
    //@}

    /**
     * @brief Returns the distance between two wires in the specified view
     * @param w1 index of the first wire
     * @param w2 index of the second wire
     * @param p plane number within the TPC
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     * @return the distance between the two wires
     *
     * This method assumes that all the wires on all the planes on the specified
     * view of all TPCs have the same pitch.
     */
    geo::Length_t WirePitch(geo::View_t view) const;


    //@{
    /**
     * @brief Returns the angle of the wires in the specified view from vertical
     * @param view the view
     * @param TPC the index of the TPC in the specified cryostat
     * @param Cryo the cryostat
     * @param tpcid ID of the TPC
     * @return the angle [radians]
     * @throw cet::exception ("GeometryCore" category) if no such view
     *
     * The angle is defined as in WireGeo::ThetaZ().
     *
     * This method assumes all wires in the view have the same angle (it queries
     * for the first).
     *
     * @deprecated This does not feel APA-ready
     */
    double WireAngleToVertical(geo::View_t view, geo::TPCID const& tpcid) const;
    double WireAngleToVertical(geo::View_t view, int TPC=0, int Cryo=0) const
      { return WireAngleToVertical(view, geo::TPCID(Cryo, TPC)); }
    //@}

    /// @} Wire access and information



    /**
     * @name Wire geometry queries
     *
     * Please note the differences between functions:
     * ChannelsIntersect(), WireIDsIntersect() and IntersectionPoint()
     * all calculate wires intersection using the same equation.
     * ChannelsIntersect() and WireIdsIntersect() will return true
     * if the two wires cross, return false if they don't.
     * IntersectionPoint() does not check if the two wires cross.
     */
    /// @{

    //
    // simple geometry queries
    //

    /**
     * @brief Fills two arrays with the coordinates of the wire end points
     * @param wireid ID of the wire
     * @param xyzStart (output) an array with the start coordinate
     * @param xyzEnd (output) an array with the end coordinate
     * @throws cet::exception wire not present
     *
     * The starting point is the wire end with lower z coordinate.
     *
     * @deprecated use the wire ID interface instead (but note that it does not
     *             sort the ends)
     */
    void WireEndPoints
      (geo::WireID const& wireid, double *xyzStart, double *xyzEnd) const;

    /**
     * @brief Fills two arrays with the coordinates of the wire end points
     * @param cstat cryostat number
     * @param tpc tpc number within the cryostat
     * @param plane plane number within the TPC
     * @param wire wire number within the plane
     * @param xyzStart (output) an array with the start coordinate
     * @param xyzEnd (output) an array with the end coordinate
     * @throws cet::exception wire not present
     *
     * The starting point is the wire end with lower z coordinate.
     *
     * @deprecated use the wire ID interface instead (but note that it does not
     *             sort the ends)
     */
    void WireEndPoints(
      unsigned int cstat, unsigned int tpc, unsigned int plane, unsigned int wire,
      double *xyzStart, double *xyzEnd
      ) const
      { WireEndPoints(geo::WireID(cstat, tpc, plane, wire), xyzStart, xyzEnd); }

    //@{
    /**
     * @brief Returns a segment whose ends are the wire end points
     * @param wireid ID of the wire
     * @return a segment whose ends are the wire end points
     * @throws cet::exception wire not present
     *
     * The start and end are assigned as returned from the geo::WireGeo object.
     * The rules for this assignment are documented in that class.
     *
     * @deprecated use the wire ID interface instead (but note that it does not
     *             sort the ends)
     */
    template <typename Point>
    Segment<Point> WireEndPoints(geo::WireID const& wireID) const;
    Segment<DefaultPoint_t> WireEndPoints(geo::WireID const& wireID) const
      { return WireEndPoints<DefaultPoint_t>(wireID); }

    //@}

    //
    // closest wire
    //

    /**
     * @brief Returns the ID of wire closest to position in the specified TPC.
     * @param point the point to be tested [cm]
     * @param planeid ID of the plane
     * @return the ID of the wire, or an invalid wire ID
     * @see `geo::PlaneGeo::ClosestWireID()`
     * @bug Instead of returning an invalid wire ID, an exception is thrown!
     *
     * If the nearest wire is not closer than half a wire pitch, the result is
     * marked invalid. The returned (invalid) ID will contain the non-existing
     * wire that would be the nearest, if it existed.
     *
     * If the wire ID is invalid and the existing closest wire is desired,
     * a possible solution is (when the BUG will be solved):
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::WireID wireID = geom->NearestWireID(point, planeID);
     * if (!wireID) wireID = geom->Plane(planeID).ClosestWireID(wireID);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Note however that this will execute plane lookup twice, and a more
     * efficient approach would be to ask the plane everything directly:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::PlaneGeo const& plane = geom->Plane(planeID);
     * geo::WireID wireID = plane.NearestWireID(point);
     * if (!wireID) wireID = plane.ClosestWireID(wireID);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Until the BUG is fixed, the actual working code is:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::PlaneGeo const& plane = geom->Plane(planeID);
     * geo::WireID wireID;
     * try {
     *   wireID = plane.NearestWireID(point);
     * }
     * catch (geo::InvalidWireError const& e) {
     *   if (!e.hasSuggestedWire()) throw;
     *   wireID = plane.ClosestWireID(e.suggestedWireID());
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    geo::WireID NearestWireID
      (geo::Point_t const& point, geo::PlaneID const& planeid) const;

    //@{
    /**
     * @brief Returns the ID of wire closest to position in the specified TPC.
     * @param point the point to be tested [cm]
     * @param planeid ID of the plane
     * @param PlaneNo plane number within the TPC
     * @param TPCNo tpc number within the cryostat
     * @param cstat cryostat number
     * @return the ID of the wire, or an invalid wire ID
     * @bug Instead of returning an invalid wire ID, an exception is thrown!
     *
     * The different versions allow different way to provide the position.
     *
     * @deprecated Use the version with a `geo::Point_t` and `PlaneID` arguments
     * @todo remove the integers version
     */
    geo::WireID NearestWireID
      (const double point[3], geo::PlaneID const& planeid) const;
    geo::WireID  NearestWireID
      (std::vector<double> const& point, geo::PlaneID const& planeid)  const;
    geo::WireID   NearestWireID
      (const TVector3& point, geo::PlaneID const& planeid) const
      { return NearestWireID(geo::vect::toPoint(point), planeid); }
    geo::WireID   NearestWireID(const double point[3],
                                      unsigned int const PlaneNo,
                                      unsigned int const TPCNo = 0,
                                      unsigned int const cstat = 0)  const
     { return NearestWireID(point, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    geo::WireID   NearestWireID(std::vector<double> const& point,
                                      unsigned int const PlaneNo,
                                      unsigned int const TPCNo = 0,
                                      unsigned int const cstat = 0)  const
     { return NearestWireID(point, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    geo::WireID   NearestWireID(const TVector3& point,
                                      unsigned int const PlaneNo,
                                      unsigned int const TPCNo = 0,
                                      unsigned int const cstat = 0)  const
     { return NearestWireID(point, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    geo::WireID   NearestWireID(geo::Point_t const& point,
                                      unsigned int const PlaneNo,
                                      unsigned int const TPCNo = 0,
                                      unsigned int const cstat = 0)  const
     { return NearestWireID(point, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    //@}

    /**
     * @brief Returns the index of wire closest to position in the specified TPC
     * @param point the point to be tested [cm]
     * @param planeid ID of the plane
     * @return the index of the wire, or `geo::WireID::InvalidID` on failure
     * @bug Actually, on failure an exception `geo::InvalidWireError` is thrown
     *
     * @deprecated Use NearestWireID() instead.
     */
    geo::WireID::WireID_t NearestWire
      (geo::Point_t const& point, geo::PlaneID const& planeid) const;

    //@{
    /**
     * @brief Returns the index of wire closest to position in the specified TPC
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param planeid ID of the plane
     * @param PlaneNo plane number within the TPC
     * @param TPCNo tpc number within the cryostat
     * @param cstat cryostat number
     * @return the index of the wire
     *
     * The different versions allow different way to provide the position.
     *
     * @deprecated Use NearestWireID() instead.
     * @todo remove the integers version
     * @todo what happens when no wire is found?
     */
    unsigned int       NearestWire
      (const double worldLoc[3], geo::PlaneID const& planeid)  const;
    unsigned int       NearestWire
      (std::vector<double> const& worldLoc, geo::PlaneID const& planeid) const;
    unsigned int       NearestWire
      (const TVector3& worldLoc, geo::PlaneID const& planeid)  const
      { return NearestWire(geo::vect::toPoint(worldLoc), planeid); }
    unsigned int       NearestWire(const double worldLoc[3],
                                   unsigned int const PlaneNo,
                                   unsigned int const TPCNo = 0,
                                   unsigned int const cstat = 0)  const
      { return NearestWire(worldLoc, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    unsigned int       NearestWire(std::vector<double> const& worldLoc,
                                   unsigned int const PlaneNo,
                                   unsigned int const TPCNo = 0,
                                   unsigned int const cstat = 0) const
      { return NearestWire(worldLoc, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    unsigned int       NearestWire(const TVector3& worldLoc,
                                   unsigned int const PlaneNo,
                                   unsigned int const TPCNo = 0,
                                   unsigned int const cstat = 0)  const
      { return NearestWire(worldLoc, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    unsigned int       NearestWire(geo::Point_t const& worldLoc,
                                   unsigned int const PlaneNo,
                                   unsigned int const TPCNo = 0,
                                   unsigned int const cstat = 0)  const
      { return NearestWire(worldLoc, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    //@}


    /**
     * @brief Returns the index of the nearest wire to the specified position
     * @param YPos y coordinate on the wire plane
     * @param ZPos z coordinate on the wire plane
     * @param planeid ID of the plane
     * @return an index interpolation between the two nearest wires
     * @deprecated Use
     *             `WireCoordinate(geo::Point_t const&, geo::PlaneID const&)`
     *             instead
     *
     * Respect to `NearestWireID()`, this method returns a real number,
     * representing a continuous coordinate in the wire axis, with the round
     * values corresponding to the actual wires.
     *
     * @todo Unify (y, z) coordinate
     */
    geo::Length_t WireCoordinate
      (double YPos, double ZPos, geo::PlaneID const& planeid) const;

    /**
     * @brief Returns the index of the nearest wire to the specified position
     * @param YPos y coordinate on the wire plane
     * @param ZPos z coordinate on the wire plane
     * @param PlaneNo number of plane
     * @param TPCNo number of TPC
     * @param cstat number of cryostat
     * @return an index interpolation between the two nearest wires
     * @see ChannelMapAlg::WireCoordinate()
     *
     * @deprecated Use the version with plane ID instead
     */
    geo::Length_t WireCoordinate(double YPos, double ZPos,
                                 unsigned int PlaneNo,
                                 unsigned int TPCNo,
                                 unsigned int cstat) const
      { return WireCoordinate(YPos, ZPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }

    //@{
    /**
     * @brief Returns the index of the nearest wire to the specified position
     * @param pos world coordinates of the position (it will be projected)
     * @param planeid ID of the plane
     * @return an index interpolation between the two nearest wires
     * @see ChannelMapAlg::WireCoordinate()
     *
     * Respect to NearestWireID(), this method returns a real number,
     * representing a continuous coordinate in the wire axis, with the round
     * values corresponding to the actual wires.
     */
    geo::Length_t WireCoordinate
      (geo::Point_t const& pos, geo::PlaneID const& planeid) const;
    geo::Length_t WireCoordinate
      (TVector3 const& pos, geo::PlaneID const& planeid) const
      { return WireCoordinate(geo::vect::toPoint(pos), planeid); }
    //@}

    //
    // wire intersections
    //

    // The following functions are utilized to determine if two wires
    // in the TPC intersect or not, and if they do then
    // determine the coordinates of the intersection.

    /**
     * @brief Computes the intersection between two lines on a plane
     * @param A_start_x x coordinate of one point of the first segment
     * @param A_start_y y coordinate of one point of the first segment
     * @param A_end_x x coordinate of another point of the first segment
     * @param A_end_y y coordinate of another point of the first segment
     * @param B_start_x x coordinate of one point of the second segment
     * @param B_start_y y coordinate of one point of the second segment
     * @param B_end_x x coordinate of another point of the second segment
     * @param B_end_y y coordinate of another point of the second segment
     * @param x _(output)_ variable to store the x coordinate of intersection
     * @param y _(output)_ variable to store the y coordinate of intersection
     * @return whether intersection exists
     *
     * The order of the ends is not relevant.
     * The return value is `false` if the two segments are parallel.
     * In that case, `x` and `y` variables are not changed.
     * Otherwise, they hold the intersection coordinate, even if the
     * intersection point is beyond one or both the segments.
     */
    bool IntersectLines(
      double A_start_x, double A_start_y, double A_end_x, double A_end_y,
      double B_start_x, double B_start_y, double B_end_x, double B_end_y,
      double& x, double& y
      ) const;

    /**
     * @brief Computes the intersection between two segments on a plane
     * @param A_start_x x coordinate of the start of the first segment
     * @param A_start_y y coordinate of the start of the first segment
     * @param A_end_x x coordinate of the end of the first segment
     * @param A_end_y y coordinate of the end of the first segment
     * @param B_start_x x coordinate of the start of the second segment
     * @param B_start_y y coordinate of the start of the second segment
     * @param B_end_x x coordinate of the end of the second segment
     * @param B_end_y y coordinate of the end of the second segment
     * @param x _(output)_ variable to store the x coordinate of intersection
     * @param y _(output)_ variable to store the y coordinate of intersection
     * @return whether intersection exists and is on both segments
     *
     * The order of the ends is not relevant.
     * The return value is `false` if the two segments are parallel, or if their
     * intersection point is not on _both_ the segments.
     * If the segments are parallel, x and y variables are not changed.
     * Otherwise, they hold the intersection coordinate, even if the
     * intersection point is beyond one or both the segments.
     */
    bool IntersectSegments(
      double A_start_x, double A_start_y, double A_end_x, double A_end_y,
      double B_start_x, double B_start_y, double B_end_x, double B_end_y,
      double& x, double& y
      ) const;

    //@{
    /**
     * @brief Computes the intersection between two wires.
     * @param wid1 ID of the first wire
     * @param wid2 ID of the other wire
     * @param intersection (output) the intersection point (global coordinates)
     * @return whether an intersection was found inside the TPC the wires belong
     *
     * The "intersection" refers to the projection of the wires into the same
     * wire plane. The coordinate along the drift direction is arbitrarily set
     * to the one of the first wire.
     * Wires are assumed to have at most one intersection.
     * If wires are parallel, `intersection` will have all components set to
     * infinity (`std::numeric_limits<>::infinity()`) and `false` is returned.
     * If the intersection is outside the TPC, `false` is also returned, but the
     * `intersection` point will contain that intersection.
     *
     * To test that the result is not infinity (nor NaN), use
     * `geo::vect::isfinite(intersection)` etc.
     *
     */
    bool WireIDsIntersect(
      WireID const& wid1, WireID const& wid2,
      geo::Point_t& intersection
      ) const;
    bool WireIDsIntersect
      (WireID const& wid1, WireID const& wid2, TVector3& intersection) const;
    //@}

    /**
     * @brief Computes the intersection between two wires.
     * @param wid1 ID of the first wire
     * @param wid2 ID of the other wire
     * @param widIntersect (output) the coordinate of the intersection point
     * @return whether an intersection was found within the TPC
     *
     * The "intersection" refers to the projection of the wires into the same
     * @f$ x = 0 @f$ plane.
     * Wires are assumed to have at most one intersection.
     * If wires are parallel, `widIntersect` will have the two components set to
     * infinity (`std::numeric_limits<>::infinity()`) and the TPC number set to
     * invalid (`geo::TPCID::InvalidID`). Also, `false` is returned.
     * If the intersection is outside the TPC, `false` is also returned, but the
     * `widIntersect` will contain the coordinates of that intersection. The TPC
     * number is still set to invalid, although the intersection _might_ belong
     * to a valid TPC somewhere else.
     *
     *
     * @deprecated This method uses arbitrary assumptions and should not be
     *             used. Use the interface returning a full vector instead.
     */
    bool WireIDsIntersect
      (WireID const& wid1, WireID const& wid2, WireIDIntersection& widIntersect)
      const;

    /**
     * @brief Returns the intersection point of two wires
     * @param wid1 ID of the first wire
     * @param wid2 ID of the other wire
     * @param y (output) y coordinate of the intersection point
     * @param z (output) z coordinate of the intersection point
     * @return whether an intersection was found within the TPC
     * @see WireIDsIntersect()
     *
     * The behaviour of this method reflects the one of `WireIDsIntersect()`,
     * which supersedes this one.
     *
     * To test if the result is infinity, use e.g. `std::isfinite(y)`.
     *
     * @deprecated This method uses arbitrary assumptions and should not be
     *             used. Use `WireIDsIntersect()` returning a vector, instead.
     */
    bool IntersectionPoint(geo::WireID const& wid1,
                           geo::WireID const& wid2,
                           double &y,
                           double &z) const;

    /**
     * @brief Returns the intersection point of two wires
     * @param wire1 wire index of the first wire
     * @param wire2 wire index of the other wire
     * @param plane1 plane index of the first wire
     * @param plane2 plane index of the other wire
     * @param cstat cryostat number
     * @param tpc tpc number within the cryostat where the planes belong
     * @param y (output) y coordinate of the intersection point
     * @param z (output) z coordinate of the intersection point
     * @return whether an intersection was found
     *
     * No check is performed, not any information provided, about the validity
     * of the result.
     *
     * @deprecated This method uses arbitrary assumptions and should not be
     *             used. Use `WireIDsIntersect()` returning a vector, instead.
     */
    bool IntersectionPoint(unsigned int wire1,
                           unsigned int wire2,
                           unsigned int plane1,
                           unsigned int plane2,
                           unsigned int cstat,
                           unsigned int tpc,
                           double &y,
                           double &z) const
      {
        return IntersectionPoint(
          geo::WireID(cstat, tpc, plane1, wire1),
          geo::WireID(cstat, tpc, plane2, wire2),
          y, z
          );
      }


    /**
     * @brief Returns the plane that is not in the specified arguments
     * @param pid1 a plane
     * @param pid2 another plane
     * @return the ID to the third plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     * @throws cet::exception (category: "GeometryCore") if pid1 and pid2 match
     *
     * This function requires a geometry with exactly three planes.
     * If the two input planes are not on the same TPC, the result is undefined.
     */
    geo::PlaneID ThirdPlane
      (geo::PlaneID const& pid1, geo::PlaneID const& pid2) const;


    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param pid2 ID of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @param output_plane ID of the plane on which to calculate the slope
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if input planes match
     *
     * Given a slope as projected in two planes, returns the slope as projected
     * in the specified output plane.
     * The slopes are defined in uniform units; they should be computed as
     * distance ratios (or tangent of a geometrical angle; the formula is still
     * valid using dt/dw directly in case of equal wire pitch in all planes
     * and uniform drift velocity.
     */
    double ThirdPlaneSlope(geo::PlaneID const& pid1, double slope1,
                           geo::PlaneID const& pid2, double slope2,
                           geo::PlaneID const& output_plane) const;

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param pid2 ID of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a slope as projected in two planes, returns the slope as projected
     * in the third plane.
     * This function is a shortcut assuming exactly three wire planes in the
     * TPC, in which case the output plane is chosen as the one that is neither
     * of the input planes.
     */
    double ThirdPlaneSlope(geo::PlaneID const& pid1, double slope1,
                           geo::PlaneID const& pid2, double slope2) const;

    //@{
    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param plane1 index of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param plane2 index of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @param tpcid TPC where the two planes belong
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a slope as projected in two planes, returns the slope as projected
     * in the third plane.
     */
    double ThirdPlaneSlope(geo::PlaneID::PlaneID_t plane1, double slope1,
                           geo::PlaneID::PlaneID_t plane2, double slope2,
                           geo::TPCID const& tpcid) const
      {
        return ThirdPlaneSlope(
          geo::PlaneID(tpcid, plane1), slope1,
          geo::PlaneID(tpcid, plane2), slope2
          );
      }
    double ThirdPlaneSlope(unsigned int plane1, double slope1,
                           unsigned int plane2, double slope2,
                           unsigned int tpc, unsigned int cstat) const
      {
        return ThirdPlaneSlope
          (plane1, slope1, plane2, slope2, geo::TPCID(cstat, tpc));
      }
    //@}


    /**
     * @brief Returns dT/dW on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first dT/dW
     * @param dTdW1 dT/dW as seen on the first plane
     * @param pid2 ID of the plane of the second dT/dW
     * @param dTdW2 dT/dW  as seen on the second plane
     * @param output_plane ID of the plane on which to calculate the slope
     * @return dT/dW on the third plane, or -999. if dT/dW would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a dT/dW as projected in two planes, returns the dT/dW as projected
     * in the third plane.
     * The dT/dW are defined in time ticks/wide number units.
     */
    double ThirdPlane_dTdW(geo::PlaneID const& pid1, double slope1,
                           geo::PlaneID const& pid2, double slope2,
                           geo::PlaneID const& output_plane) const;

    /**
     * @brief Returns dT/dW on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first dT/dW
     * @param dTdW1 dT/dW as seen on the first plane
     * @param pid2 ID of the plane of the second dT/dW
     * @param dTdW2 dT/dW  as seen on the second plane
     * @return dT/dW on the third plane, or -999. if dT/dW would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a dT/dW as projected in two planes, returns the dT/dW as projected
     * in the third plane.
     * This function is a shortcut assuming exactly three wire planes in the
     * TPC, in which case the output plane is chosen as the one that is neither
     * of the input planes.
     */
    double ThirdPlane_dTdW(geo::PlaneID const& pid1, double slope1,
                           geo::PlaneID const& pid2, double slope2) const;


    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param angle1 angle or the wires on the first plane
     * @param slope1 slope as observed on the first plane
     * @param angle2 angle or the wires on the second plane
     * @param slope2 slope as observed on the second plane
     * @param angle_target angle or the wires on the target plane
     * @return the slope as measure on the third plane, or 999 if infinity
     *
     * This function will return a small slope if both input slopes are small.
     */
    static double ComputeThirdPlaneSlope(
      double angle1, double slope1,
      double angle2, double slope2,
      double angle_target
      );

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param angle1 angle or the wires on the first plane
     * @param pitch1 wire pitch on the first plane
     * @param dTdW1 slope in dt/dw units as observed on the first plane
     * @param angle2 angle or the wires on the second plane
     * @param pitch2 wire pitch on the second plane
     * @param dTdW2 slope in dt/dw units as observed on the second plane
     * @param angle_target angle or the wires on the target plane
     * @param pitch_target wire pitch on the target plane
     * @return dt/dw slope as measured on the third plane, or 999 if infinity
     *
     * The input slope must be specified in dt/dw non-homogeneous coordinates.
     *
     * This function will return a small slope if both input slopes are small.
     */
    static double ComputeThirdPlane_dTdW(
      double angle1, double pitch1, double dTdW1,
      double angle2, double pitch2, double dTdW2,
      double angle_target, double pitch_target
      );

    /// @} Wire geometry queries



    /**
     * @name Optical detector geometry access and information
     * @anchor GeometryCoreOpDetGeometry
     * @see @ref GeometryCoreOpDetChannel "optical detector channel information"
     *
     * There are a number of ways to identify an optical detector or channel:
     *
     * * geometric:
     *     * cryostat (e.g. `geo::CryostatID`) and relative optical detector
     *       number within it
     *     * unique optical detector number
     * * readout:
     *     * optical detector channel
     *     * "hardware" channel
     *
     * And they all should be better documented!
     */
    /// @{

    //
    // group features
    //

    //@{
    /// Number of OpDets in the whole detector
    unsigned int NOpDets() const;
    //@}


    //
    // access
    //
    
    //@{
    /**
     * @brief Returns the optical detector with the specified ID.
     * @param oid ID of the optical detector
     * @return the optical detector `oid`
     * @throw cet::exception (category `"GeometryCore"`) if `oid` ID is invalid
     * @throw cet::exception (category `"GeometryCore"`) if `oid` does not
     *        represent any optical detector
     */
    geo::OpDetGeo const& OpDet(geo::OpDetID const& oid) const;
    //@}
    
    
    //@{
    /**
     * @brief Returns the `geo::OpDetGeo` object for the given channel number.
     * @param OpChannel optical detector unique channel number
     * @see GeometryCoreOpDetGeometry "optical detector identification"
     */
    OpDetGeo const& OpDetGeoFromOpChannel(unsigned int OpChannel) const;
    //@}

    //@{
    /**
     * @brief Returns the `geo::OpDetGeo` object for the given detector number.
     * @param OpDet optical detector unique number
     * @see GeometryCoreOpDetGeometry "optical detector identification"
     */
    OpDetGeo const& OpDetGeoFromOpDet(unsigned int OpDet) const;
    //@}


    //@{
    /**
     * @brief Find the nearest OpChannel to some point
     * @param xyz point to be queried, in world coordinates
     * @return the nearest OpChannel to the point,
     *         or `std::numeric_limits<unsigned int>::max()` if invalid point
     *
     * @deprecated This method does not tell in which cryostat the detector is;
     *             use `geo::CryostatGeo::GetClosestOpDet()` instead
     *             (find the cryostat with `PositionToCryostatPtr()`).
     *
     */
    unsigned int GetClosestOpDet(geo::Point_t const& point) const;
    unsigned int GetClosestOpDet(double const* point) const;
    //@}
    
    
    //@{
    /**
     * @brief Returns the optical detector in the cryostat closest to `point`.
     * @param oid ID of the optical detector
     * @return the TPC closest to the optical detector `oid`
     * @see `geo::CryostatGeo::GetTPCclosestToOpDet()`
     * @throw cet::exception (category `"GeometryCore"`) if `oid` invalid
     * 
     * The TPC closest to the center of the specified optical detector is
     * returned.
     * 
     * See `geo::CryostatGeo::GetTPCclosestToOpDet()` for details.
     */
    geo::TPCGeo const& GetTPCclosestToOpDet(geo::OpDetID const& oid) const;
    //@}


    //
    // object description
    //

    //@{
    /**
     * @brief Returns gdml string which gives sensitive opdet name
     * @param c ID of the cryostat the detector is in
     *
     * This name is defined in the geometry (GDML) description.
          *
     * @todo Change to use CryostatID
     */
    std::string OpDetGeoName(unsigned int c = 0) const;
    //@}

    /// @} Optical detector access and information



    /// @name Auxiliary detectors access and information
    /// @{

    /// @todo use a AutDetID_t instead of unsigned int?

    //
    // group features
    //

    /**
     * @brief Returns the number of auxiliary detectors
     *
     * This method returns the total number of scintillator paddles
     * (Auxiliary Detectors aka AuxDet) outside of the cryostat
     *
     * @todo Change return type to size_t
     */
    unsigned int NAuxDets() const { return AuxDets().size(); }

    /**
     * @brief Returns the number of sensitive components of auxiliary detector
     * @param aid ID of the auxiliary detector
     * @return number of sensitive components in the auxiliary detector aid
     * @thrws cet::exception (category "Geometry") if aid does not exist
     */
    unsigned int NAuxDetSensitive(size_t const& aid) const;

    //
    // access
    //

    /**
     * @brief Returns the specified auxiliary detector
     * @param ad the auxiliary detector index
     * @return a constant reference to the specified auxiliary detector
     *
     * @todo what happens if it does not exist?
     * @todo remove the default parameter?
     */
    AuxDetGeo const& AuxDet(unsigned int const ad = 0) const;

    /**
     * @brief Returns the index of the auxiliary detector at specified location.
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @return the index of the detector, or
     *        `std::numeric_limits<unsigned int>::max()` if no detector is there
     *
     * @bug Actually, an exception is thrown.
     * @deprecated Use the version with `geo::Point_t`.
     */
    unsigned int FindAuxDetAtPosition(double const worldLoc[3]) const;

    /**
     * @brief Returns the index of the auxiliary detector at specified location.
     * @param point location to be tested
     * @return the index of the detector, or
     *        `std::numeric_limits<unsigned int>::max()` if no detector is there
     *
     * @bug Actually, an exception is thrown.
     */
    unsigned int FindAuxDetAtPosition(geo::Point_t const& point) const;

    /**
     * @brief Fills the indices of the sensitive auxiliary detector at location
     * @param point location to be tested
     * @param adg _(output)_ auxiliary detector index
     * @param sv _(output)_ sensitive volume index
     */
    void  FindAuxDetSensitiveAtPosition(geo::Point_t const& point,
                                        std::size_t       & adg,
                                        std::size_t       & sv) const;

    /**
     * @brief Fills the indices of the sensitive auxiliary detector at location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param adg (output) auxiliary detector index
     * @param sv (output) sensitive volume index
     * @deprecated Use the version with `geo::Point_t`.
     */
    void  FindAuxDetSensitiveAtPosition(double const worldLoc[3],
                                        size_t     & adg,
                                        size_t     & sv) const;

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param point location to be tested
     * @param ad _(output)_ the auxiliary detector index
     * @return constant reference to AuxDetGeo object of the auxiliary detector
     *
     * @todo what happens if it does not exist?
     */
    AuxDetGeo const& PositionToAuxDet
      (geo::Point_t const& point, unsigned int& ad) const;

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param ad (output) the auxiliary detector index
     * @return constant reference to AuxDetGeo object of the auxiliary detector
     *
     * @deprecated Use the version with `geo::Point_t`.
     * @todo what happens if it does not exist?
     */
    AuxDetGeo const& PositionToAuxDet
      (double const worldLoc[3], unsigned int& ad) const;

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param point location to be tested
     * @param ad _(output)_ the auxiliary detector index
     * @param sv _(output)_ the auxiliary detector sensitive volume index
     * @return reference to AuxDetSensitiveGeo object of the auxiliary detector
     *
     * @todo what happens if it does not exist?
     */
    const AuxDetSensitiveGeo& PositionToAuxDetSensitive
      (geo::Point_t const& point, size_t& ad, size_t& sv) const;

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param ad (output) the auxiliary detector index
     * @param sv (output) the auxiliary detector sensitive volume index
     * @return reference to AuxDetSensitiveGeo object of the auxiliary detector
     *
     * @todo what happens if it does not exist?
     * @deprecated Use the version with `geo::Point_t`.
     */
    const AuxDetSensitiveGeo& PositionToAuxDetSensitive(double const worldLoc[3],
                                                        size_t     & ad,
                                                        size_t     & sv) const;

    const AuxDetGeo&         ChannelToAuxDet(std::string const& auxDetName,
                                             uint32_t    const& channel) const; // return the AuxDetGeo for the given detector
                                                                                // name and channel

    const AuxDetSensitiveGeo& ChannelToAuxDetSensitive(std::string const& auxDetName,
                                                       uint32_t    const& channel) const; // return the AuxDetSensitiveGeo for the given

    /// @} Auxiliary detectors access and information



    /// @name TPC readout channels and views
    /// @{

    //
    // group features
    //

    /// Returns the number of TPC readout channels in the detector
    unsigned int Nchannels() const;

    /// @brief Returns the number of channels in the specified ROP
    /// @return number of channels in the specified ROP, 0 if non-existent
    unsigned int Nchannels(readout::ROPID const& ropid) const;

    /// @brief Returns an std::vector<ChannelID_t> in all TPCs in a TPCSet
    std::vector<raw::ChannelID_t> ChannelsInTPCs() const;
    //
    /**
     * @brief Returns a list of possible views in the detector.
     * @return the set of views
     */
    std::set<geo::View_t> const& Views() const { return allViews; }


    //
    // access
    //

    /**
     * @brief Returns whether the specified channel exists and is valid
     * @param channel the ID of the channel
     * @return whether the specified channel exists
     *
     * A channel is defined as existing and valid if its ID is not invalid and
     * if the channel is physical.
     */
    bool HasChannel(raw::ChannelID_t channel) const;

    //@{
    /**
     * @brief Returns the ID of the TPC channel connected to the specified wire
     * @param plane the number of plane
     * @param wire the number of wire
     * @param tpc the number of TPC
     * @param cryostat the number of cryostat
     * @param wireid the ID of the wire
     * @return the ID of the channel, or raw::InvalidChannelID if invalid wire
     *
     * @todo Verify the raw::InvalidChannelID part
     * @todo remove the integers version
     */
    raw::ChannelID_t  PlaneWireToChannel(WireID const& wireid) const;
    raw::ChannelID_t  PlaneWireToChannel(unsigned int const plane,
                                         unsigned int const wire,
                                         unsigned int const tpc = 0,
                                         unsigned int const cstat = 0) const
      { return PlaneWireToChannel(geo::WireID(cstat, tpc, plane, wire)); }
    //@}

    //
    // single object features
    //

    /**
     * @brief Returns the type of signal on the specified TPC channel
     * @param channel TPC channel ID
     * @return the type of signal on the specified channel, or geo::kMysteryType
     *
     * @todo verify that kMysteryType is returned on invalid channel
     */
    SigType_t SignalType(raw::ChannelID_t const channel) const;


    /**
     * @brief Returns the view (wire orientation) on the specified TPC channel
     * @param channel TPC channel ID
     * @return the type of signal on the specified channel, or geo::kUnknown
     *
     * The view of the readout plane `channel` belongs to is returned, as in
     * `View(readout::ROPID const&) const`.
     */
    View_t View(raw::ChannelID_t const channel) const;


    /**
     * @brief Returns a list of wires connected to the specified TPC channel
     * @param channel TPC channel ID
     * @return vector containing the ID of all the connected wires
     * @throws cet::exception (category: "Geometry") if non-existent channel
     */
    std::vector<geo::WireID> ChannelToWire
      (raw::ChannelID_t const channel) const;


    /// Returns the ID of the ROP the channel belongs to
    /// @throws cet::exception (category: "Geometry") if non-existent channel
    readout::ROPID ChannelToROP(raw::ChannelID_t channel) const;


    //
    // geometry queries
    //

    /**
     * @brief Returns the ID of the channel nearest to the specified position
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param planeid ID of the wire plane the channel must belong to
     * @return the ID of the channel, or `raw::InvalidChannelID` if invalid wire
     * @bug on invalid wire, a `geo::InvalidWireError` exception is thrown
     *
     */
    raw::ChannelID_t  NearestChannel
      (geo::Point_t const& worldLoc, geo::PlaneID const& planeid) const;

    //@{
    /**
     * @brief Returns the ID of the channel nearest to the specified position
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param PlaneNo the number of plane
     * @param TPCNo the number of TPC
     * @param cstat the number of cryostat
     * @return the ID of the channel, or raw::InvalidChannelID if invalid wire
     * @bug on invalid wire, a `geo::InvalidWireError` exception is thrown
     *
     * The different versions allow different way to provide the position.
     *
     * @todo remove the integers version
     */
    raw::ChannelID_t  NearestChannel
      (const double worldLoc[3], geo::PlaneID const& planeid) const;
    raw::ChannelID_t  NearestChannel
      (std::vector<double> const& worldLoc, geo::PlaneID const& planeid) const;
    raw::ChannelID_t  NearestChannel
      (const TVector3& worldLoc, geo::PlaneID const& planeid) const
      { return NearestChannel(geo::vect::toPoint(worldLoc), planeid); }
    raw::ChannelID_t  NearestChannel(const double worldLoc[3],
                                       unsigned int const PlaneNo,
                                       unsigned int const TPCNo = 0,
                                       unsigned int const cstat = 0) const
      { return NearestChannel(worldLoc, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    raw::ChannelID_t  NearestChannel(std::vector<double> const& worldLoc,
                                       unsigned int const PlaneNo,
                                       unsigned int const TPCNo = 0,
                                       unsigned int const cstat = 0) const
      { return NearestChannel(worldLoc, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    raw::ChannelID_t  NearestChannel(const TVector3& worldLoc,
                                       unsigned int const PlaneNo,
                                       unsigned int const TPCNo = 0,
                                       unsigned int const cstat = 0) const
      { return NearestChannel(worldLoc, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    raw::ChannelID_t  NearestChannel(geo::Point_t const& worldLoc,
                                       unsigned int const PlaneNo,
                                       unsigned int const TPCNo = 0,
                                       unsigned int const cstat = 0) const
      { return NearestChannel(worldLoc, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
    //@}

    /**
     * @brief Returns an intersection point of two channels
     * @param c1 one channel ID
     * @param c2 the other channel ID
     * @param y (output) y coordinate of the intersection
     * @param z (output) z coordinate of the intersection
     * @return whether a intersection point was found
     *
     * @todo what happens for channels from different TPCs?
     * @todo what happens for channels with multiple intersection points?
     *
     * @deprecated This is clearly not APA-aware
     */
    bool ChannelsIntersect
      (raw::ChannelID_t c1, raw::ChannelID_t c2, double &y, double &z) const;

    /// @} TPC readout channels



    /// @name TPC set information
    /// @{

    //
    // group features
    //

    //@{
    /**
     * @brief Returns the total number of TPC sets in the specified cryostat
     * @param cryoid cryostat ID
     * @return number of TPC sets in the cryostat, or 0 if no cryostat found
     *
     * The NSiblingElements() method is overloaded and its
     * return depends on the type of ID.
     */
    unsigned int NTPCsets(readout::CryostatID const& cryoid) const;
    unsigned int NSiblingElements(readout::TPCsetID const& tpcsetid) const
      { return NTPCsets(tpcsetid); }
    //@}

    /// Returns the largest number of TPC sets any cryostat in the detector has
    unsigned int MaxTPCsets() const;


    /**
     * @brief Returns a container with one entry per TPC set.
     * @tparam T type of data in the container
     * @return a container with one default-constructed `T` per TPC set
     * @see `readout::TPCsetDataContainer`
     *
     * The working assumption is that all cryostats have the same number of
     * TPC sets. It is always guaranteed that all existing TPC sets have an
     * entry in the container, although if the previous working assumption is
     * not satisfied there will be entries in the containers which are not
     * associated to a valid TPC set.
     *
     * The interface of the container is detailed in the documentation of the
     * container itself, `readout::TPCsetDataContainer`. Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto tracksPerTPCset
     *   = geom->makeTPCsetData<std::vector<recob::Track const*>>();
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    readout::TPCsetDataContainer<T> makeTPCsetData() const
      { return { Ncryostats(), MaxTPCsets() }; }

    /**
     * @brief Returns a container with one entry per TPC set.
     * @tparam T type of data in the container
     * @param defValue the initial value of all elements in the container
     * @return a container with a value `defValue` per each TPC set
     * @see `readout::TPCsetDataContainer`
     *
     * This function operates as `makeTPCsetData() const`, except that it copies
     * the specified value into all the entries of the container. Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto nTracksPerTPCset = geom->makeTPCsetData(0U);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    readout::TPCsetDataContainer<T> makeTPCsetData(T const& defValue) const
      { return { Ncryostats(), MaxTPCsets(), defValue }; }


    //
    // access
    //
    /// Returns whether we have the specified TPC set
    /// @return whether the TPC set is valid and exists
    bool HasTPCset(readout::TPCsetID const& tpcsetid) const;

    /// Returns whether we have the specified TPC set
    bool HasElement(readout::TPCsetID const& tpcsetid) const
      { return HasTPCset(tpcsetid); }


    /**
     * @brief Returns the ID of the TPC set at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @return the TPC set ID, or an invalid one if no TPC set is there
     */
    readout::TPCsetID FindTPCsetAtPosition(double const worldLoc[3]) const;

    //
    // mapping
    //
    /// Returns the ID of the TPC set tpcid belongs to
    readout::TPCsetID TPCtoTPCset(geo::TPCID const& tpcid) const;

    /**
     * @brief Returns a list of ID of TPCs belonging to the specified TPC set
     * @param tpcsetid ID of the TPC set to convert into TPC IDs
     * @return the list of TPCs, empty if TPC set is invalid
     *
     * Note that the check is performed on the validity of the TPC set ID, that
     * does not necessarily imply that the TPC set specified by the ID actually
     * exists. Check the existence of the TPC set first (HasTPCset()).
     * Behaviour on valid, non-existent TPC set IDs is undefined.
     */
    std::vector<geo::TPCID> TPCsetToTPCs
      (readout::TPCsetID const& tpcsetid) const;


    ///
    /// iterators
    ///

    /// Initializes the specified ID with the ID of the first TPC set
    void GetBeginID(readout::TPCsetID& id) const
      { GetBeginID(id.asCryostatID()); id.TPCset = 0; }

    /// Initializes the specified ID with the invalid ID after the last TPC set
    void GetEndID(readout::TPCsetID& id) const
      { GetEndID(id.asCryostatID()); id.TPCset = 0; }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(readout::TPCsetID& id) const; // inline implementation

    /// Returns the ID of the first TPC set in the specified cryostat.
    readout::TPCsetID GetBeginTPCsetID(geo::CryostatID const& id) const
      { return { id, 0 }; }

    /// Returns the (possibly invalid) ID after the last TPC set of the
    /// specified cryostat.
    readout::TPCsetID GetEndTPCsetID(geo::CryostatID const& id) const
      { return { id.Cryostat + 1, 0 }; }


    /// Returns an iterator pointing to the first TPC set ID in the detector
    TPCset_id_iterator begin_TPCset_id() const
      { return TPCset_id_iterator(this, TPCset_id_iterator::begin_pos); }

    /// Returns an iterator pointing after the last TPC set ID in the detector
    TPCset_id_iterator end_TPCset_id() const
      { return TPCset_id_iterator(this, TPCset_id_iterator::end_pos); }

    /// Returns an iterator pointing to the first TPC set ID in the specified
    /// cryostat.
    TPCset_id_iterator begin_TPCset_id(geo::CryostatID const& cid) const
      { return TPCset_id_iterator(this, GetBeginTPCsetID(cid)); }

    /// Returns an iterator pointing after the last TPC set ID in the specified
    /// cryostat.
    TPCset_id_iterator end_TPCset_id(geo::CryostatID const& cid) const
      { return TPCset_id_iterator(this, GetEndTPCsetID(cid)); }

    /**
     * @brief Enables ranged-for loops on all TPC set IDs of the detector
     * @returns an object suitable for ranged-for loops on all TPC set IDs
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (readout::TPCsetID const& sID: geom->IterateTPCsetIDs()) {
     *
     *   // useful code here
     *
     * } // for all TPC sets
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     */
    IteratorBox<
      TPCset_id_iterator,
      &GeometryCore::begin_TPCset_id, &GeometryCore::end_TPCset_id
      >
    IterateTPCsetIDs() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all TPC set IDs of the specified
     *        cryostat.
     * @param cid the ID of the cryostat to loop the TPC set IDs of
     * @returns an object suitable for ranged-for loops on TPC set IDs
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::CryostatID cid{1}; // cryostat #1 (hope it exists!)
     * for (readout::TPCsetID const& tID: geom->IterateTPCsetIDs(cid)) {
     *
     *   // useful code here
     *
     * } // for all TPC sets in cryostat #1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      TPCset_id_iterator, geo::CryostatID,
      &GeometryCore::begin_TPCset_id, &GeometryCore::end_TPCset_id
      >
    IterateTPCsetIDs(geo::CryostatID const& cid) const { return { this, cid }; }

    /// @} TPC set information


    /// @name Readout plane information
    /// @{

    //
    // group features
    //

    //@{
    /**
     * @brief Returns the total number of ROP in the specified TPC set
     * @param tpcsetid TPC set ID
     * @return number of readout planes in the TPC set, or 0 if no TPC set found
     *
     * Note that this methods explicitly check the existence of the TPC set.
     *
     * The NSiblingElements() method is overloaded and its
     * return depends on the type of ID.
     */
    unsigned int NROPs(readout::TPCsetID const& tpcsetid) const;
    unsigned int NSiblingElements(readout::ROPID const& ropid) const
      { return NROPs(ropid); }
    //@}

    /// Returns the largest number of ROPs a TPC set in the detector has
    unsigned int MaxROPs() const;


    /**
     * @brief Returns a container with one entry per readout plane.
     * @tparam T type of data in the container
     * @return a container with one default-constructed `T` per readout plane
     * @see `readout::ROPDataContainer`
     *
     * The working assumption is that all cryostats have the same number of
     * TPC sets, and all TPC sets have the same number of readout planes.
     * It is always guaranteed that all existing readout planes have an entry
     * in the container, although if the previous working assumption is not
     * satisfied there will be entries in the container which are not
     * associated to a valid readout plane.
     *
     * The interface of the container is detailed in the documentation of the
     * container itself, `readout::ROPDataContainer`. Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto hitsPerROP
     *   = geom->makeROPdata<std::vector<recob::Hit const*>>();
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    readout::ROPDataContainer<T> makeROPdata() const
      { return { Ncryostats(), MaxTPCsets(), MaxROPs() }; }

    /**
     * @brief Returns a container with one entry per readout plane.
     * @tparam T type of data in the container
     * @param defValue the initial value of all elements in the container
     * @return a container with one default-constructed `T` per readout plane
     * @see `readout::ROPDataContainer`
     *
     * This function operates as `makeROPdata() const`, except that copies
     * the specified value into all the entries of the container. Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto nHitsPerROP = geom->makeROPdata(0U);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    readout::ROPDataContainer<T> makeROPdata(T const& defValue) const
      { return { Ncryostats(), MaxTPCsets(), MaxROPs(), defValue }; }


    //
    // access
    //
    /// Returns whether we have the specified readout plane
    /// @return whether the readout plane is valid and exists
    bool HasROP(readout::ROPID const& ropid) const;

    /// Returns whether we have the specified readout plane
    /// @return whether the readout plane is valid and exists
    bool HasElement(readout::ROPID const& ropid) const { return HasROP(ropid); }


    //
    // mapping
    //
    /**
     * @brief Returns the ID of the ROP planeid belongs to
     * @param planeid ID of the wire plane
     * @return the ID of the ROP planeid belongs to
     *
     * If planeid is an invalid ID, an invalid ROP ID is returned.
     * If planeid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent wire plane, the result is undefined.
     * Use HasPlaneID() to check if the wire plane actually exists.
     */
    readout::ROPID WirePlaneToROP(geo::PlaneID const& planeid) const;

    /**
     * @brief Returns a list of ID of planes belonging to the specified ROP
     * @param ropid ID of the readout plane
     * @return list of ID of wire planes belonging to the specified ROP
     *
     * If ropid is an invalid ID, an empty list is returned.
     * If ropid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent readout plane, the result is undefined.
     * Use HasROP() to check if the readout plane actually exists.
     */
    std::vector<geo::PlaneID> ROPtoWirePlanes
      (readout::ROPID const& ropid) const;

    /**
     * @brief Returns a list of ID of TPCs the specified ROP spans
     * @param ropid ID of the readout plane
     * @return the list of TPC IDs, empty if readout plane ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour on non-existing readout planes is undefined.
     */
    std::vector<geo::TPCID> ROPtoTPCs(readout::ROPID const& ropid) const;


    /**
     * @brief Returns the ID of the first channel in the specified readout plane
     * @param ropid ID of the readout plane
     * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour for non-existing readout planes is undefined.
     */
    raw::ChannelID_t FirstChannelInROP(readout::ROPID const& ropid) const;

    ///
    /// iterators
    ///

    /// Initializes the specified ID with the ID of the first readout plane.
    void GetBeginID(readout::ROPID& id) const
      { GetBeginID(id.asTPCsetID()); id.ROP = 0; }

    /// Initializes the specified ID with the invalid ID after the last ROP.
    void GetEndID(readout::ROPID& id) const
      { GetEndID(id.asTPCsetID()); id.ROP = 0; }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(readout::ROPID& id) const; // inline implementation

    /// Returns the ID of the first readout plane of the specified cryostat.
    readout::ROPID GetBeginROPID(geo::CryostatID const& id) const
      { return { GetBeginTPCsetID(id), 0 }; }

    /// Returns the (possibly invalid) ID after the last readout plane of the
    /// specified cryostat.
    readout::ROPID GetEndROPID(geo::CryostatID const& id) const
      { return { GetEndTPCsetID(id), 0 }; }

    /// Returns the ID of the first readout plane of the specified TPC set.
    readout::ROPID GetBeginROPID(readout::TPCsetID const& id) const
      { return { id, 0 }; }

    /// Returns the (possibly invalid) ID after the last readout plane of the
    /// specified TPC set.
    readout::ROPID GetEndROPID(readout::TPCsetID const& id) const
      { return { GetNextID(id), 0 }; }

    /// Returns an iterator pointing to the first ROP ID in the detector.
    ROP_id_iterator begin_ROP_id() const
      { return ROP_id_iterator(this, ROP_id_iterator::begin_pos); }

    /// Returns an iterator pointing after the last ROP ID in the detector.
    ROP_id_iterator end_ROP_id() const
      { return ROP_id_iterator(this, ROP_id_iterator::end_pos); }

    /// Returns an iterator pointing to the first readout plane ID in the
    /// specified cryostat.
    ROP_id_iterator begin_ROP_id(geo::CryostatID const& ID) const
      { return ROP_id_iterator(this, GetBeginROPID(ID)); }

    /// Returns an iterator pointing after the last readout plane ID in the
    /// specified cryostat.
    ROP_id_iterator end_ROP_id(geo::CryostatID const& ID) const
      { return ROP_id_iterator(this, GetEndROPID(ID)); }

    /// Returns an iterator pointing to the first readout plane ID in the
    /// specified TPC set.
    ROP_id_iterator begin_ROP_id(readout::TPCsetID const& ID) const
      { return ROP_id_iterator(this, GetBeginROPID(ID)); }

    /// Returns an iterator pointing after the last readout plane ID in the
    /// specified TPC set.
    ROP_id_iterator end_ROP_id(readout::TPCsetID const& ID) const
      { return ROP_id_iterator(this, GetEndROPID(ID)); }

    /**
     * @brief Enables ranged-for loops on all readout plane IDs of the detector.
     * @returns an object suitable for ranged-for loops on all ROP IDs
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (readout::ROPID const& rID: geom->IterateROPIDs()) {
     *
     *   // useful code here
     *
     * } // for all ROPs
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     */
    IteratorBox<
      ROP_id_iterator,
      &GeometryCore::begin_ROP_id, &GeometryCore::end_ROP_id
      >
    IterateROPIDs() const { return { this }; }

    /**
     * @brief Enables ranged-for loops on all readout plane IDs of the specified
     *        cryostat.
     * @param cid the ID of the cryostat to loop the readout plane IDs of
     * @returns an object suitable for ranged-for loops on readout plane IDs
     *
     * If the cryostat ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::CryostatID cid{1}; // cryostat #1 (hope it exists!)
     * for (readout::ROPID const& rID: geom->IterateROPIDs(cid)) {
     *
     *   // useful code here
     *
     * } // for all readout planes in cryostat #1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      ROP_id_iterator, geo::CryostatID,
      &GeometryCore::begin_ROP_id, &GeometryCore::end_ROP_id
      >
    IterateROPIDs(geo::CryostatID const& cid) const { return { this, cid }; }

    /**
     * @brief Enables ranged-for loops on all readout plane IDs of the specified
     *        TPC set.
     * @param sid the ID of the TPC set to loop the readout plane IDs of
     * @returns an object suitable for ranged-for loops on readout plane IDs
     *
     * If the TPC set ID is invalid, the effect is undefined.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * readout::TPCsetID sid{ 0, 1 }; // C:0 S:1 (hope it exists!)
     * for (readout::ROPID const& rID: geom->IterateROPIDs(sid)) {
     *
     *   // useful code here
     *
     * } // for all readout planes in C:0 S:1
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    LocalIteratorBox<
      ROP_id_iterator, readout::TPCsetID,
      &GeometryCore::begin_ROP_id, &GeometryCore::end_ROP_id
      >
    IterateROPIDs(readout::TPCsetID const& sid) const { return { this, sid }; }


    /**
     * @brief Returns the view of the channels in the specified readout plane
     * @param ropid readout plane ID
     * @return the type of signal on the specified ROP
     *
     * Returns the view (wire orientation) on the channels of specified readout
     * plane.
     * If ropid is an invalid ID, geo::kUnknown is returned.
     * If ropid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent readout plane, the result is undefined.
     * Use HasROP() to check if the readout plane actually exists.
     */
    geo::View_t View(readout::ROPID const& ropid) const;

    /**
     * @brief Returns the type of signal of channels in specified readout plane
     * @param ropid readout plane ID
     * @return the type of signal on the specified ROP
     *
     * Assumes that all the channels on the readout plane have the same signal
     * type.
     * If ropid is an invalid ID, geo::kMysteryType is returned.
     * If ropid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent readout plane, the result is undefined.
     * Use HasROP() to check if the readout plane actually exists.
     */
    geo::SigType_t SignalType(readout::ROPID const& ropid) const;


    /// @} Readout plane information



    /**
     * @name Optical readout channels
     * @anchor GeometryCoreOpDetChannel
     * @see @ref GeometryCoreOpDetGeometry "optical detector geometry information"
     */
    /// @todo add explanation of the different IDs
    /// @{

    //
    // group features
    //

    /// Number of electronics channels for all the optical detectors
    unsigned int NOpChannels() const;

    /// Largest optical channel number
    unsigned int MaxOpChannel() const;

    // Number of hardware channels for a given optical detector
    unsigned int NOpHardwareChannels(int opDet) const;


    //
    // access
    //

    /// Is this a valid OpChannel number?
    bool IsValidOpChannel(int opChannel) const;

    /// Convert detector number and hardware channel to unique channel
    unsigned int OpChannel(int detNum, int hardwareChannel) const;

    /// Convert unique channel to detector number
    unsigned int OpDetFromOpChannel(int opChannel) const;

    /// Convert unique channel to hardware channel
    unsigned int HardwareChannelFromOpChannel(int opChannel) const;

    /// Get unique opdet number from cryo and internal count
    unsigned int OpDetFromCryo(unsigned int o, unsigned int c) const;

    /// @} Optical readout channels


    //
    // unsorted methods
    //

    /**
     * @brief Returns whether a value is within the specified range
     * @param value the value to be tested
     * @param min the lower boundary
     * @param max the upper boundary
     * @return whether the value is within range
     *
     * If min is larger than max, they are swapped.
     * A tolerance of 10^-6 (absolute) is used.
     *
     * @todo Use wiggle instead of 10^-6
     * @todo resort source code for a bit of speed up
     */
    bool ValueInRange(double value, double min, double max) const;


    /// @name Geometry initialization
    /// @{

    /**
     * @brief Loads the geometry information from the specified files
     * @param gdmlfile path to file to be used for Geant4 simulation
     * @param rootfile path to file for internal geometry representation
     * @param builder algorithm to be used for the interpretation of geometry
     * @param bForceReload reload even if there is already a valid geometry
     * @see ApplyChannelMap()
     *
     * Both paths must directly resolve to an available file, as no search
     * is performed for them.
     *
     * The gdmlfile parameter does not have to necessarily be in GDML format,
     * as long as it's something supported by Geant4. This file is not used by
     * the geometry, but its path is provided on request by the simulation
     * modules (see LArSoft `LArG4` module).
     * The rootfile also does not need to be a ROOT file, but just anything
     * that TGeoManager::Import() supports. This file is parsed immediately
     * and the internal geometry representation is built out of it.
     *
     * @note After calling this method, the detector geometry information can
     * be considered complete, but the geometry service provider is not fully
     * initialized yet, since it's still necessary to provide or update the
     * channel mapping.
     */
    void LoadGeometryFile(
      std::string gdmlfile, std::string rootfile,
      geo::GeometryBuilder& builder,
      bool bForceReload = false
      );

    /**
     * @brief Loads the geometry information from the specified files
     * @param gdmlfile path to file to be used for Geant4 simulation
     * @param rootfile path to file for internal geometry representation
     * @param bForceReload reload even if there is already a valid geometry
     * @see ApplyChannelMap()
     *
     * This legacy version of `LoadGeometryFile()` uses a standard
     * `geo::GeometryBuilder` implementation.
     * Do not rely on it if you can avoid it.
     */
    void LoadGeometryFile
      (std::string gdmlfile, std::string rootfile, bool bForceReload = false);

    /**
     * @brief Initializes the geometry to work with this channel map
     * @param pChannelMap a pointer to the channel mapping algorithm to be used
     * @see LoadGeometryFile()
     *
     * The specified channel mapping is used with this geometry.
     * These modifications typically involve some resorting of the objects.
     *
     * This method needs to be called after LoadGeometryFile() to complete the
     * geometry initialization.
     */
    void ApplyChannelMap(std::unique_ptr<geo::ChannelMapAlg> pChannelMap);
    /// @}


  protected:
    /// Sets the detector name
    void SetDetectorName(std::string const& new_name) { fDetectorName = new_name; }

    /// Returns the object handling the channel map
    geo::ChannelMapAlg const* ChannelMap() const
      { return fChannelMapAlg.get(); }

    //@{
    /// Return the internal cryostat list
    CryostatList_t&       Cryostats()       { return fGeoData.cryostats; }
    CryostatList_t const& Cryostats() const { return fGeoData.cryostats; }
    //@}

    //@{
    /// Return the interfal auxiliary detectors list
    AuxDetList_t&       AuxDets()       { return fGeoData.auxDets; }
    AuxDetList_t const& AuxDets() const { return fGeoData.auxDets; }
    //@}

  private:

    GeometryData_t fGeoData;        ///< The detector description data

    double         fSurfaceY;       ///< The point where air meets earth for this detector.
    std::string    fDetectorName;   ///< Name of the detector.
    std::string    fGDMLfile;       ///< path to geometry file used for Geant4 simulation
    std::string    fROOTfile;       ///< path to geometry file for geometry in GeometryCore
    double         fMinWireZDist;   ///< Minimum distance in Z from a point in which
                                    ///< to look for the closest wire
    double         fPositionWiggle; ///< accounting for rounding errors when testing positions

    /// Configuration for the geometry builder
    /// (needed since builder is created after construction).
    fhicl::ParameterSet fBuilderParameters;
    std::unique_ptr<const geo::ChannelMapAlg> fChannelMapAlg;
                                    ///< Object containing the channel to wire mapping

    // cached values
    std::set<geo::View_t> allViews; ///< All views in the detector.


    std::vector<TGeoNode const*> FindDetectorEnclosure
      (std::string const& name = "volDetEnclosure") const;

    bool FindFirstVolume
      (std::string const& name, std::vector<const TGeoNode*>& path) const;

    /// Parses ROOT geometry nodes and builds LArSoft geometry representation.
    /// @param builder the algorithm to be used
    void BuildGeometry(geo::GeometryBuilder& builder);

    /// Wire ID check for WireIDsIntersect methods
    bool WireIDIntersectionCheck
      (const geo::WireID& wid1, const geo::WireID& wid2) const;

    /// Returns whether x and y are within both specified ranges (A and B).
    static bool PointWithinSegments(
      double A_start_x, double A_start_y, double A_end_x, double A_end_y,
      double B_start_x, double B_start_y, double B_end_x, double B_end_y,
      double x, double y
      );

    /// Runs the sorting of geometry with the sorter provided by channel mapping
    void SortGeometry(geo::GeoObjectSorter const& sorter);

    /// Performs all the updates needed after sorting
    void UpdateAfterSorting();

    /// Deletes the detector geometry structures
    void ClearGeometry();

    /// Throws an exception ("GeometryCore" category) unless pid1 and pid2
    /// are on different planes of the same TPC (ID validity is not checked)
    static void CheckIndependentPlanesOnSameTPC
      (geo::PlaneID const& pid1, geo::PlaneID const& pid2, const char* caller);

  }; // class GeometryCore



  /** **************************************************************************
   * @brief Iterator to navigate through all the nodes
   *
   * Note that this is not a fully standard forward iterator in that it lacks
   * of the postfix operator. The reason is that it's too expensive and it
   * should be avoided.
   * Also I did not bother declaring the standard type definitions
   * (that's just laziness).
   *
   * An example of iteration:
   *
   *     TGeoNode const* pCurrentNode;
   *
   *     ROOTGeoNodeForwardIterator iNode(geom->ROOTGeoManager()->GetTopNode());
   *     while ((pCurrentNode = *iNode)) {
   *       // do something with pCurrentNode
   *       ++iNode;
   *     } // while
   *
   * These iterators are one use only, and they can't be reset after a loop
   * is completed.
   */
  class ROOTGeoNodeForwardIterator {
      public:
    /// Constructor: start from this node
    ROOTGeoNodeForwardIterator(TGeoNode const* start_node)
      { init(start_node); }

    /// Returns the pointer to the current node, or nullptr if none
    TGeoNode const* operator* () const
      { return current_path.empty()? nullptr: current_path.back().self; }

    /// Points to the next node, or to nullptr if there are no more
    ROOTGeoNodeForwardIterator& operator++ ();

    /// Returns the full path of the current node
    std::vector<TGeoNode const*> get_path() const;

      protected:
    using Node_t = TGeoNode const*;
    struct NodeInfo_t {
      Node_t self; int sibling;
      NodeInfo_t(Node_t new_self, int new_sibling)
        : self(new_self), sibling(new_sibling) {}
    }; // NodeInfo_t

    /// which node, which sibling?
    std::vector<NodeInfo_t> current_path;

    void reach_deepest_descendant();

    void init(TGeoNode const* start_node);

  }; // class ROOTGeoNodeForwardIterator

  /// @}
  // END Geometry group --------------------------------------------------------

} // namespace geo



//******************************************************************************
//*** inline implementation
//***
inline bool geo::GeometryCore::IncrementID(geo::CryostatID& id) const {
  ++id.Cryostat;
  if (id) id.isValid = HasCryostat(id); // if invalid already, it stays so
  return bool(id);
} // geo::GeometryCore::IncrementID(geo::CryostatID)

inline bool geo::GeometryCore::IncrementID(geo::TPCID& id) const {
  unsigned int const nTPCsInCryo = NTPC(id);
  if (++id.TPC < nTPCsInCryo) return bool(id); // if was invalid, it stays so
  // no more TPCs in this cryostat
  id.TPC = 0;
  return IncrementID(id.asCryostatID()); // also sets validity
} // geo::GeometryCore::IncrementID(geo::TPCID)

inline bool geo::GeometryCore::IncrementID(geo::PlaneID& id) const {
  // this implementation is non-optimal, in that the cryostat lookup is
  // performed both here and, potentially, in IncrementID(TPCID)
  unsigned int const nPlanesInTPC = Nplanes(id);
  if (++id.Plane < nPlanesInTPC) return bool(id); // if was invalid, stays so
  // no more planes in this TPCs
  id.Plane = 0;
  return IncrementID(id.asTPCID()); // also sets validity
} // geo::GeometryCore::IncrementID(geo::PlaneID)

inline bool geo::GeometryCore::IncrementID(geo::WireID& id) const {
  // this implementation is non-optimal, in that the TPC lookup is
  // performed both here and, potentially, in IncrementID(PlaneID)
  unsigned int const nWiresInPlane = Nwires(id);
  if (++id.Wire < nWiresInPlane) return bool(id); // if was invalid, stays so
  // no more wires in this plane
  id.Wire = 0;
  return IncrementID(id.asPlaneID()); // also sets validity
} // geo::GeometryCore::IncrementID(geo::WireID)

inline bool geo::GeometryCore::IncrementID(readout::TPCsetID& id) const {
  unsigned int const nTPCsetsInCryo = NTPCsets(id);
  if (++id.TPCset < nTPCsetsInCryo)
    return bool(id); // if was invalid, it stays so
  // no more TPC sets in this cryostat
  id.TPCset = 0;
  return IncrementID(id.asCryostatID()); // also sets validity
} // geo::GeometryCore::IncrementID(readout::TPCsetID)

inline bool geo::GeometryCore::IncrementID(readout::ROPID& id) const {
  // this implementation is non-optimal, in that the cryostat lookup is
  // performed both here and, potentially, in IncrementID(TPCsetID)
  unsigned int const nROPinTPC = NROPs(id);
  if (++id.ROP < nROPinTPC) return bool(id); // if was invalid, stays so
  // no more readout planes in this TPC set
  id.ROP = 0;
  return IncrementID(id.asTPCsetID()); // also sets validity
} // geo::GeometryCore::IncrementID(readout::ROPID)



//******************************************************************************
//***  template implementation
//***
//------------------------------------------------------------------------------
template <typename Point>
geo::GeometryCore::Segment<Point> geo::GeometryCore::WireEndPoints
  (geo::WireID const& wireid) const
{
  geo::WireGeo const& wire = Wire(wireid);
  return { wire.GetStart<Point>(), wire.GetEnd<Point>() };
} // geo::GeometryCore::WireEndPoints(WireID)


//------------------------------------------------------------------------------
template <typename Stream>
void geo::GeometryCore::Print
  (Stream&& out, std::string indent /* = "  " */) const
{

  out << "Detector " << DetectorName() << " has "
    << Ncryostats() << " cryostats and "
    << NAuxDets() << " auxiliary detectors:";

  auto const& detEnclosureBox = DetectorEnclosureBox();
  out << "\n" << indent << "Detector enclosure: "
    << detEnclosureBox.Min() << " -- " << detEnclosureBox.Max()
    << " cm => ( " << detEnclosureBox.SizeX() << " x "
    << detEnclosureBox.SizeY() << " x "
    << detEnclosureBox.SizeZ() << " ) cm^3"
    ;

  for (geo::CryostatGeo const& cryostat: IterateCryostats()) {
    out << "\n" << indent;
    cryostat.PrintCryostatInfo
      (std::forward<Stream>(out), indent + "  ", cryostat.MaxVerbosity);

    const unsigned int nTPCs = cryostat.NTPC();
    for(unsigned int t = 0;  t < nTPCs; ++t) {
      const geo::TPCGeo& tpc = cryostat.TPC(t);

      out << "\n" << indent << "  ";
      tpc.PrintTPCInfo
        (std::forward<Stream>(out), indent + "    ", tpc.MaxVerbosity);

      const unsigned int nPlanes = tpc.Nplanes();
      for(unsigned int p = 0; p < nPlanes; ++p) {
        const geo::PlaneGeo& plane = tpc.Plane(p);
        const unsigned int nWires = plane.Nwires();

        out << "\n" << indent << "    ";
        plane.PrintPlaneInfo
          (std::forward<Stream>(out), indent + "      ", plane.MaxVerbosity);
        geo::SigType_t const sigType = SignalType(plane.ID());
        out << "\n" << indent << "      "
          << "signal type: " << SignalTypeName(sigType)
          << " (" << static_cast<int>(sigType) << ")";

        for(unsigned int w = 0;  w < nWires; ++w) {
          const geo::WireGeo& wire = plane.Wire(w);
          geo::WireID wireID(plane.ID(), w);

          // the wire should be aligned on z axis, half on each side of 0,
          // in its local frame
          out << "\n" << indent << "      " << wireID << " ";
          wire.PrintWireInfo
            (std::forward<Stream>(out), indent + "      ", wire.MaxVerbosity);
        } // for wire
      } // for plane
    } // for TPC

    unsigned int const nOpDets = cryostat.NOpDet();
    for (unsigned int iOpDet = 0; iOpDet < nOpDets; ++iOpDet) {
      geo::OpDetGeo const& opDet = cryostat.OpDet(iOpDet);
      out << "\n" << indent << "  [OpDet #" << iOpDet << "] "
        << "(closest: " << cryostat.GetTPCclosestToOpDet(opDet.ID()).ID()
        << ") ";
      
      opDet.PrintOpDetInfo
        (std::forward<Stream>(out), indent + "  ", opDet.MaxVerbosity);
    } // for
  } // for cryostat

  unsigned int const nAuxDets = NAuxDets();
  for (unsigned int iDet = 0; iDet < nAuxDets; ++iDet) {
    geo::AuxDetGeo const& auxDet = AuxDet(iDet);

    out << "\n" << indent << "[#" << iDet << "] ";
    auxDet.PrintAuxDetInfo
      (std::forward<Stream>(out), indent + "  ", auxDet.MaxVerbosity);

    unsigned int const nSensitive = auxDet.NSensitiveVolume();
    switch (nSensitive) {
      case 0: break;
      case 1: {
        geo::AuxDetSensitiveGeo const& auxDetS = auxDet.SensitiveVolume(0U);
        out << "\n" << indent << "  ";
        auxDetS.PrintAuxDetInfo
          (std::forward<Stream>(out), indent + "    ", auxDetS.MaxVerbosity);
        break;
      }
      default:
        for (unsigned int iSens = 0; iSens < nSensitive; ++iSens) {
          out << "\n" << indent << "[#" << iSens << "] ";
          geo::AuxDetSensitiveGeo const& auxDetS
            = auxDet.SensitiveVolume(iSens);
          auxDetS.PrintAuxDetInfo
            (std::forward<Stream>(out), indent + "    ", auxDetS.MaxVerbosity);
        } // for
        break;
    } // if sensitive detectors

  } // for auxiliary detector

  out << '\n';

} // geo::GeometryCore::Print()


//------------------------------------------------------------------------------
// template member function specializations
namespace geo {

  template <>
  inline geo::TPCID GeometryCore::GetBeginID<geo::TPCID, geo::CryostatID>
    (geo::CryostatID const& id) const
    { return GetBeginTPCID(id); }

  template <>
  inline geo::TPCID GeometryCore::GetEndID<geo::TPCID, geo::CryostatID>
    (geo::CryostatID const& id) const
    { return GetEndTPCID(id); }

  template <>
  inline geo::PlaneID GeometryCore::GetBeginID<geo::PlaneID, geo::CryostatID>
    (geo::CryostatID const& id) const
    { return GetBeginPlaneID(id); }

  template <>
  inline geo::PlaneID GeometryCore::GetEndID<geo::PlaneID, geo::CryostatID>
    (geo::CryostatID const& id) const
    { return GetEndPlaneID(id); }

} // namespace geo

//******************************************************************************
//
// geo::details::cryostat_id_iterator_base<>
//
template <typename GEOID>
inline geo::details::cryostat_id_iterator_base<GEOID>::operator bool() const
  { return geometry() && geometry()->HasElement(localID()); }

template <typename GEOID>
inline auto geo::details::cryostat_id_iterator_base<GEOID>::get() const
  -> ElementPtr_t
  { return geometry()->GetElementPtr(localID()); }

template <typename GEOID>
inline void geo::details::cryostat_id_iterator_base<GEOID>::set_local_limits()
  { limit = geometry()->NSiblingElements(localID()); }

template <typename GEOID>
inline void geo::details::cryostat_id_iterator_base<GEOID>::set_begin()
  { geometry()->GetBeginID(ID()); }

template <typename GEOID>
inline void geo::details::cryostat_id_iterator_base<GEOID>::set_end()
  { geometry()->GetEndID(ID()); }

template <typename GEOID>
void geo::details::cryostat_id_iterator_base<GEOID>::next() {
  if (at_end()) return;
  if (++local_index() < limit) return;
  localID().isValid = false;
} // geo::cryostat_id_iterator_base<GEOID>::next()


//
// geo::details::TPC_id_iterator_base<>
//
template <typename GEOID>
inline geo::details::TPC_id_iterator_base<GEOID>::operator bool() const {
  return upper_iterator::geometry()
    && upper_iterator::geometry()->HasElement(localID());
} // geo::details::TPC_id_iterator_base<>::operator bool()


template <typename GEOID>
inline
auto geo::details::TPC_id_iterator_base<GEOID>::get() const -> ElementPtr_t
  { return upper_iterator::geometry()->GetElementPtr(localID()); }

template <typename GEOID>
inline void geo::details::TPC_id_iterator_base<GEOID>::set_local_limits() {
  // limit is how many sibling TPCs there are
  limit = upper_iterator::geometry()->NSiblingElements(localID());
} // geo::details::TPC_id_iterator_base<GEOID>::set_local_limits()

template <typename GEOID>
inline void geo::details::TPC_id_iterator_base<GEOID>::next() {
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
  set_local_limits();
} // geo::details::TPC_id_iterator_base<GEOID>::next()


//
// geo::details::plane_id_iterator_base<>
//
template <typename GEOID>
inline geo::details::plane_id_iterator_base<GEOID>::operator bool() const {
  return upper_iterator::geometry()
    && upper_iterator::geometry()->HasElement(localID());
} // geo::details::plane_id_iterator_base<>::operator bool()


template <typename GEOID>
inline auto geo::details::plane_id_iterator_base<GEOID>::get() const
  -> ElementPtr_t
  { return upper_iterator::geometry()->GetElementPtr(localID()); }

template <typename GEOID>
inline void geo::details::plane_id_iterator_base<GEOID>::set_local_limits() {
  // limit is how many sibling planes there are
  limit = upper_iterator::geometry()->NSiblingElements(localID());
} // geo::details::plane_id_iterator_base<GEOID>::set_local_limits()

template <typename GEOID>
inline void geo::details::plane_id_iterator_base<GEOID>::next() {
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
  set_local_limits();
} // geo::details::plane_id_iterator_base<GEOID>::next()


//
// geo::details::wire_id_iterator_base<>
//
template <typename GEOID>
inline geo::details::wire_id_iterator_base<GEOID>::operator bool() const {
  return upper_iterator::geometry()
    && upper_iterator::geometry()->HasElement(localID());
} // geo::details::wire_id_iterator_base<>::operator bool()

template <typename GEOID>
inline auto geo::details::wire_id_iterator_base<GEOID>::get() const
  -> ElementPtr_t
  { return upper_iterator::geometry()->GetElementPtr(localID()); }

template <typename GEOID>
inline void geo::details::wire_id_iterator_base<GEOID>::set_local_limits() {
  // limit is how many sibling wires there are
  limit = upper_iterator::geometry()->NSiblingElements(localID());
} // geo::details::wire_id_iterator_base<>::set_local_limits()

template <typename GEOID>
inline void geo::details::wire_id_iterator_base<GEOID>::next() {
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
  set_local_limits();
} // geo::details::wire_id_iterator_base<>::next()


//
// comparison operators between ID iterators and element iterators
//
template <typename GEOIDITER>
bool geo::details::operator==
  (geometry_element_iterator<GEOIDITER> const& iter, GEOIDITER const& id_iter)
{
  return iter.id_iterator() == id_iter;
} // operator==(iterator_t, id_iterator_t)

template <typename GEOIDITER>
bool geo::details::operator!=
  (geometry_element_iterator<GEOIDITER> const& iter, GEOIDITER const& id_iter)
{
  return iter.id_iterator() != id_iter;
} // operator!=(iterator_t, id_iterator_t)


//
// geo::details::TPCset_id_iterator_base<>
//
template <typename GEOID>
inline geo::details::TPCset_id_iterator_base<GEOID>::operator bool() const {
  return upper_iterator::geometry()
    && upper_iterator::geometry()->HasElement(localID());
} // geo::details::TPCset_id_iterator_base<>::operator bool()


template <typename GEOID>
inline void geo::details::TPCset_id_iterator_base<GEOID>::set_local_limits() {
  // limit is how many sibling TPCs there are
  limit = upper_iterator::geometry()->NSiblingElements(localID());
} // geo::details::TPCset_id_iterator_base<GEOID>::set_local_limits()

template <typename GEOID>
inline void geo::details::TPCset_id_iterator_base<GEOID>::next() {
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
  set_local_limits();
} // geo::details::TPCset_id_iterator_base<GEOID>::next()


//
// geo::details::ROP_id_iterator_base<>
//
template <typename GEOID>
inline geo::details::ROP_id_iterator_base<GEOID>::operator bool() const {
  return upper_iterator::geometry()
    && upper_iterator::geometry()->HasElement(localID());
} // geo::details::ROP_id_iterator_base<>::operator bool()


template <typename GEOID>
inline void geo::details::ROP_id_iterator_base<GEOID>::set_local_limits() {
  // limit is how many sibling planes there are
  limit = upper_iterator::geometry()->NSiblingElements(localID());
} // geo::details::ROP_id_iterator_base<GEOID>::set_local_limits()

template <typename GEOID>
inline void geo::details::ROP_id_iterator_base<GEOID>::next() {
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
  set_local_limits();
} // geo::details::ROP_id_iterator_base<GEOID>::next()



//******************************************************************************

#endif // LARCOREALG_GEOMETRY_GEOMETRYCORE_H
