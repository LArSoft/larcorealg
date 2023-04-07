/**
 * @file   Partitions.h
 * @brief  Classes describing partition of an area with associated data.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 13, 2017
 *
 * This is used by `geo::DriftPartitions` to describe the partition of a plane
 * across TPCs.
 */

#ifndef LARCOREALG_GEOMETRY_PARTITIONS_H
#define LARCOREALG_GEOMETRY_PARTITIONS_H

// LArSoft libraries
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::demangle(), ...
#include "larcorealg/Geometry/SimpleGeo.h"   // lar::util::simple_geo::Rectangle

// C/C++ standard libraries
#include <algorithm> // std::upper_bound()
#include <cassert>
#include <cstdlib>  // std::size_t
#include <iterator> // std::next(), ...
#include <memory>   // std::unique_ptr<>
#include <sstream>
#include <string>
#include <type_traits> // std::declval(), std::true_type...
#include <utility>     // std::move(), std::forward()
#include <vector>

/// Partition-related utilities.
namespace geo::part {

  //-------------------------------------------------------------------------
  /// A basic interface for objects owning an area.
  class AreaOwner {
  public:
    /// Type of area covered by the partition.
    using Area_t = lar::util::simple_geo::Rectangle<double>;

    /// Type of pointer to Area_t data member of type Range_t.
    using AreaRangeMember_t = Area_t::Range_t(Area_t::*);

    /// Constructor: sets the covered area and no subpartitions.
    AreaOwner(Area_t const& area) : myArea(area) {}

    /// Returns whether the specified point is covered by this object.
    bool contains(double w, double d) const { return area().contains(w, d); }

    /// Returns the covered area.
    Area_t const& area() const { return myArea; }

    /// Output the owned area into an output stream.
    template <typename Stream>
    void dumpArea(Stream&& out) const
    {
      std::forward<Stream>(out) << area();
    }

  private:
    Area_t myArea; ///< Covered area.

  }; // class AreaOwner

  namespace details {

    /// Ordering class to sort partition by specified range (lower boundary).
    template <AreaOwner::AreaRangeMember_t Range>
    struct PartitionSorterByAreaRangeLower;

  } // namespace details

  //*************************************************************************
  //***  Fundamental classes
  //-------------------------------------------------------------------------

  /**
     * @brief Class providing custom dump for data contained in the partition.
     * @tparam Data the type of data in the partition
     *
     * This data type is expected to do its printout in a constructor with
     * signature:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * template <typename Stream>
     * PartitionDataDescriber(
     *   Stream&& out, Data const* data,
     *   std::string indent = "", std::string firstIndent = ""
     *   );
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * This class provides a terse, non-specific dump of any datum.
     * Specializations are expected for specific data types.
     */
  template <typename Data>
  struct PartitionDataDescriber {
    /// Constructor; see `describePartitionData()` for argument description.
    template <typename Stream>
    PartitionDataDescriber(Stream&& out,
                           Data const* data,
                           std::string indent = "",
                           std::string firstIndent = "");
  }; // struct PartitionDataDescriber<>

  /**
     * @brief Describes a data object for `Partition::describe()` method.
     * @tparam Stream type of output stream
     * @tparam Data type of data to be described
     * @param out output stream
     * @param data pointer to the data to be described
     * @param indent indentation string applied to all lines except the first
     * @param firstIndent special indentation string for the first line
     *
     * The description of the data pointed by `data` is printed into the output stream
     * `out`. The first line of the output is indented with `firstIndent` string (by
     * default: no indent), while all the others are indented with the string `indent` (by
     * default, also no indent).  The last output line is _not_ ended by a end-of-line
     * character.
     *
     * This function relies on `PartitionDataDescriber` template class, that should be
     * specialized to customize the printout of known data types.
     */
  template <typename Stream, typename Data>
  void describePartitionData(Stream&& out,
                             Data const* data,
                             std::string indent = "",
                             std::string firstIndent = "");

  //--------------------------------------------------------------------------
  /**
     * @brief Non-template definitions and data for `Partition` class hierarchy.
     *
     * The partition base class provides a common non-templated ground for all `Partition`
     * hierarchies.  The class defines an area that the partition cover, as a rectangle.
     * The dimensions of this rectangle, called "width" and "depth", don't have to match
     * any axis from any 3D coordinate system.
     */
  class PartitionBase : public AreaOwner {
  public:
    // imported types
    using Area_t = AreaOwner::Area_t;
    using AreaRangeMember_t = AreaOwner::AreaRangeMember_t;

    /// Constructor: sets the covered area and no subpartitions.
    PartitionBase(Area_t const& area) : AreaOwner(area) {}

  protected:
    /// Returns a description of the partition area.
    std::string describeArea(std::string indent, std::string firstIndent) const;

  }; // class PartitionBase

  //--------------------------------------------------------------------------
  /**
     * @brief Base element of a partitioned structure.
     * @tparam Data type of data contained in the partition
     *
     * An area partition is represented by a hierarchy of partition objects, each one
     * containing a non-owning pointer to the data pertaining its area.
     *
     * The partition classes all derive from this `Partition` class, which provides access
     * to the data, the area definition and a few utility functions:
     *
     * * mainly for debugging purposes, the hierarchy can be dumped into a stream by the
     *   `describe()` method
     * * method `walk()` applies a specified operation to all the partitions in the
     *   hierarchy, in an undefined order
     * * method `atPoint()` returns the data of the partition covering the specified
     *   location
     *
     * The hierarchy is such that a partition object (`Partition`) knows or can find of
     * all the subpartitions it contains, but it does not have any knowledge of the
     * partition containing it (if any).
     *
     * The partition classes do not provide algorithms to establish their relations.
     */
  template <typename Data>
  class Partition : public PartitionBase {
  public:
    using Data_t = Data;                  ///< Type of data stored in the partition.
    using Partition_t = Partition<Data>;  ///< This type.
    using Area_t = PartitionBase::Area_t; ///< Type of area.

    /// Type of list of subpartitions. It needs to preserve polymorphism.
    using Subpartitions_t = std::vector<std::unique_ptr<Partition_t const>>;

    /// Constructor: sets the covered area and no subpartitions.
    Partition(Area_t const& area) : PartitionBase(area) {}

    /// Destructor (default, virtual).
    virtual ~Partition() = default;

    /// Returns the datum directly stored (nullptr if none).
    virtual Data_t* data() const { return nullptr; }

    /**
      * @brief Returns the (sub)partition including the specified coordinates.
      * @param w width coordinate of the point
      * @param d depth coordinate of the point
      * @return a pointer to the partition, or `nullptr` if none
      *
      * This method returns a pointer to the datum associated to the area the specified
      * point (`w`, `d`) belongs to.  The area is searched for among all the included
      * partitions.
      */
    virtual Data_t* atPoint(double w, double d) const = 0;

    /// Returns a description of the partition.
    std::string describe(std::string indent, std::string firstIndent) const
    {
      return doDescribe(indent, firstIndent);
    }

    /// Returns a description of the partition.
    std::string describe(std::string indent = "") const { return describe(indent, indent); }

    /**
       * @brief Applies `pred` to all partitions.
       * @tparam Pred a predicate type
       * @param pred the predicate to be applied
       *
       * The predicate `pred` is applied to this partition first, and then to all
       * subpartitions in no specified order.
       *
       * The predicate is any object behaving like a unary function of signature:
       * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
       * void predicate(Partition<Data> const& part);
       * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       * If the predicate returns a value, that value is ignored.  The predicate is
       * forwarded while "walking" through the partitions.
       */
    template <typename Pred>
    void walk(Pred&& pred) const
    {
      walk(this, pred);
    }

    /// Returns the number of subparts in the partition (0 if simple element).
    std::size_t nParts() const { return parts().size(); }

  protected:
    static Subpartitions_t const NoSubparts; ///< Subpartitions (if any).

    /// Returns a list of all subpartitions.
    virtual Subpartitions_t const& parts() const { return NoSubparts; }

    /// Returns a description of the partition.
    virtual std::string doDescribe(std::string indent, std::string firstIndent) const
    {
      return PartitionBase::describeArea(indent, firstIndent);
    }

    /// Applies `pred` to start partition first, and then to all subpartitions.
    template <typename Pred>
    static void walk(Partition_t const* start, Pred&& pred);

  }; // class Partition<>

  //**************************************************************************
  //***  Partition class hierarchy
  //--------------------------------------------------------------------------
  /// Partition also containing data directly.
  template <typename Data>
  class PartitionWithData : public Partition<Data> {

  public:
    using Base_t = Partition<Data>;              ///< Base class.
    using Partition_t = Partition<Data>;         ///< Base type of the partition.
    using Data_t = typename Partition_t::Data_t; ///< Type of contained data.
    using Area_t = typename Partition_t::Area_t; ///< Type of covered area.

    /// Constructor: sets the covered area and the contained datum
    PartitionWithData(Area_t const& area, Data_t* myData) : Base_t(area), myData(myData) {}

    /// Returns the datum directly stored (nullptr if none).
    virtual Data_t* data() const override { return myData; }

    /// Returns stored datum only if point is covered, `nullptr` otherwise.
    virtual Data_t* atPoint(double w, double d) const override
    {
      return Base_t::contains(w, d) ? myData : nullptr;
    }

  private:
    Data_t* myData; ///< The contained datum.

    /// Returns a description of the partition.
    virtual std::string doDescribe(std::string indent, std::string firstIndent) const override;

  }; // class PartitionWithData

  //--------------------------------------------------------------------------
  /// Unpartitioned element ("leaf") of a partitioned area.
  template <typename Data>
  class PartitionElement : public PartitionWithData<Data> {
  public:
    using Base_t = PartitionWithData<Data>; ///< Base class.

    // Import constructors
    using Base_t::Base_t;
  };

  //--------------------------------------------------------------------------
  /// Partition divided in subpartitions (abstract).
  template <typename Data>
  class PartitionContainer : public PartitionWithData<Data> {
  public:
    using Base_t = PartitionWithData<Data>;
    using Partition_t = Partition<Data>; ///< Base type of the partition.

    // inherited types
    using Data_t = typename Partition_t::Data_t;
    using Area_t = typename Partition_t::Area_t;
    using Subpartitions_t = typename Partition_t::Subpartitions_t;

    /// Returns stored datum only if point is covered, `nullptr` otherwise.
    virtual Data_t* atPoint(double w, double d) const override;

  protected:
    Subpartitions_t myParts; ///< List of subpartitions.

    /// Returns the number of contained subpartitions.
    std::size_t size() const { return parts().size(); }

    /// Returns a list of the subpartitions owned.
    virtual Subpartitions_t const& parts() const override { return myParts; }

    /**
      * @brief Constructor: sets the partition.
      * @param area overall area covered
      * @param subpartitions list of subpartitions (pointers)
      * @param defData datum to be returned for points not covered by
      *        subpartitions
      *
      * The subpartitions will be moved from the argument.
      *
      * It is required and assumed that the subpartitions do not overlap and that the
      * points covered by them are a subset of `area`.  Neither of theses requirements is
      * checked.
      */
    PartitionContainer(Area_t const& area,
                       Subpartitions_t&& subpartitions,
                       Data_t* defData = nullptr)
      : Base_t(area, defData), myParts(std::move(subpartitions))
    {}

    /// Returns the only partition which could contain the specified width.
    virtual Partition_t const* findPart(double w, double d) const = 0;

    /// Describes this and each of the subpartitions.
    virtual std::string doDescribe(std::string indent, std::string firstIndent) const override;

    /// Introduction to the description of the subpartitions.
    virtual std::string describeIntro() const;

  }; // class PartitionContainer<>

  //--------------------------------------------------------------------------
  /**
     * @brief Partition of area sorted across a dimension.
     * @tparam Data type of data contained in the partition
     * @tparam Sorter type of functor providing comparison of partitions
     *
     * The sorter is a functor containing comparison functions. It must be compatible with
     * both `std::sort()` and `std::lower_bound()` functions.  The former requirement
     * implies that the sorter can compare two constant pointers to partitions. The latter
     * also implies that the sorter can compare a constant pointer to partition to a "key"
     * (a real number), and vice versa. The meaning of this comparison is not prescribed;
     * existing implementations interpret that value as a width or depth coordinate.
     *
     * A copy of the sorter is kept in this partition.
     */
  template <typename Data, typename Sorter>
  class SortedPartition : public PartitionContainer<Data> {
  public:
    using Base_t = PartitionContainer<Data>;
    using Partition_t = Partition<Data>; ///< Base type of the partition.
    using Sorter_t = Sorter;             ///< Type of sorter being used.

    // inherited data types
    using Data_t = typename Partition_t::Data_t;
    using Area_t = typename Partition_t::Area_t;
    using Subpartitions_t = typename Base_t::Subpartitions_t;

    /**
      * @brief Constructor: sets the partition.
      * @param area overall area covered
      * @param subpartitions list of subpartitions (pointers)
      * @param defData datum to be returned for points not covered by
      *        subpartitions
      * @param sorter instance of the sorter to be used
      *
      * The subpartitions will be moved from the argument and will be sorted using the
      * comparison contained in the `sorter`.  Note that this will invalidate existing
      * pointers to the sub-partitions.
      *
      * It is required and assumed that the subpartitions do not overlap and that the
      * points covered by them are a subset of `area`.  Neither of theses requirements is
      * checked.
      */
    SortedPartition(Area_t const& area,
                    Subpartitions_t&& subpartitions,
                    Data_t* defData = nullptr,
                    Sorter_t sorter = {})
      : Base_t(area, std::move(subpartitions), defData), sorter(sorter)
    {
      initParts();
    }

  protected:
    Sorter_t sorter; ///< Object used for sorting and binary search.

    /// Returns the only partition which could contain the specified `key`.
    Partition_t const* findPartWithKey(double key) const;

    /// Performs initialization on the specified subpartition list.
    void initParts();

  }; // class SortedPartition<>

  //--------------------------------------------------------------------------
  /// Partition of area along a area range dimension (width or depth).
  template <typename Data, PartitionBase::AreaRangeMember_t Range>
  class PartitionSortedByRange
    : public SortedPartition<Data, details::PartitionSorterByAreaRangeLower<Range>> {
  public:
    using Base_t = SortedPartition<Data, details::PartitionSorterByAreaRangeLower<Range>>;
    using Base_t::Base_t; // import inherited constructors
  };

  //--------------------------------------------------------------------------
  /// Partition of area along the depth dimension.
  template <typename Data>
  class DepthPartition : public PartitionSortedByRange<Data, &PartitionBase::Area_t::depth> {
  public:
    using Base_t = PartitionSortedByRange<Data, &PartitionBase::Area_t::depth>;

    using Base_t::Base_t; // import inherited constructors

    /// Returns the only partition which could contain the specified depth.
    virtual typename Base_t::Partition_t const* findPart(double /* w */, double d) const override
    {
      return Base_t::findPartWithKey(d);
    }

  private:
    virtual std::string describeIntro() const override;

  }; // class DepthPartition

  //--------------------------------------------------------------------------
  /// Partition of area along the width dimension.
  template <typename Data>
  class WidthPartition : public PartitionSortedByRange<Data, &PartitionBase::Area_t::width> {
  public:
    /// Base class.
    using Base_t = PartitionSortedByRange<Data, &PartitionBase::Area_t::width>;

    using Base_t::Base_t; // inherited constructors

    /// Returns the only partition which could contain the specified depth.
    virtual typename Base_t::Partition_t const* findPart(double w, double /* d */) const override
    {
      return Base_t::findPartWithKey(w);
    }

  private:
    virtual std::string describeIntro() const override;

  }; // class WidthPartition

  //--------------------------------------------------------------------------
  /// A container of partitions organised in a width/depth rectangular grid.
  template <typename Data>
  class GridPartition : public PartitionContainer<Data> {
  public:
    using Partition_t = Partition<Data>;     ///< Base type of the partition.
    using Base_t = PartitionContainer<Data>; ///< Base class.

    // Inherited types
    using Data_t = typename Partition_t::Data_t; ///< Type of contained data.
    using Area_t = typename Partition_t::Area_t; ///< Type of covered area.
    using Subpartitions_t = typename Base_t::Subpartitions_t;

    /**
      * @brief Creates a partition with a grid of subpartitions.
      * @param area total area covered by this partition
      * @param subpartitions all subpartitions, row by row
      * @param nDepthPartitions number of partitions on depth direction ("rows")
      * @param nWidthPartitions number of partitions on width direction ("columns")
      * @param defData partition data for areas not covered by subpartitions
      *
      * The content of the collection of subpartitions is stolen.  The subpartitions in
      * that collection are expected to be organized by row: (0;0), (0,1), (0,2)... where
      * the first index spans `nDepthPartitions` values and the second one spans
      * `nWidthPartitions` values.
      *
      */
    GridPartition(Area_t const& area,
                  Subpartitions_t&& subpartitions,
                  unsigned int nDepthPartitions,
                  unsigned int nWidthPartitions,
                  Data_t* defData = nullptr);

    /// Constructor: autodetects `nWidthPartitions` from number of subpartitions.
    GridPartition(Area_t const& area,
                  Subpartitions_t&& subpartitions,
                  unsigned int nDepthPartitions,
                  Data_t* defData = nullptr);

  private:
    std::vector<double> widthSeps; ///< Separators for width dimension.
    std::vector<double> depthSeps; ///< Separators for depth dimension.

    /// Number of partitions on width direction.
    std::size_t nWidthParts() const { return widthSeps.size(); }
    ///< Number of partitions on depth direction.
    std::size_t nDepthParts() const { return depthSeps.size(); }

    auto part(std::size_t iDepth, std::size_t iWidth) -> decltype(auto)
    {
      return Base_t::parts()[iDepth * nWidthParts() + iWidth];
    }
    auto part(std::size_t iDepth, std::size_t iWidth) const -> decltype(auto)
    {
      return Base_t::parts()[iDepth * nWidthParts() + iWidth];
    }

    /// Returns the only partition which could contain the specified depth.
    virtual Partition_t const* findPart(double w, double d) const override;

    /// Computes and returns width separation levels proper for `widthSeps`.
    std::vector<double> computeWidthSeps(unsigned int nD, unsigned int nW) const;
    /// Computes and returns width separation levels proper for `depthSeps`.
    std::vector<double> computeDepthSeps(unsigned int nD, unsigned int nW) const;

    /// Prints the information about the partition grid.
    virtual std::string doDescribe(std::string indent, std::string firstIndent) const override;

    template <PartitionBase::AreaRangeMember_t Range, typename BeginIter, typename EndIter>
    static std::vector<double> detectSeparators(BeginIter b,
                                                EndIter e,
                                                std::size_t const nGroups,
                                                std::size_t const startDelta,
                                                std::size_t const stride);

  }; // class GridPartition<>

  //--------------------------------------------------------------------------

} // namespace geo::part

//******************************************************************************
//***  inline implementation
//******************************************************************************
inline std::string geo::part::PartitionBase::describeArea(std::string,
                                                          std::string firstIndent) const
{
  std::ostringstream sstr;
  sstr << firstIndent << "partition covers ";
  dumpArea(sstr);
  return sstr.str();
}

//******************************************************************************
//***  template implementation
//******************************************************************************

namespace geo::part::details {

  //------------------------------------------------------------------------
  //---  some metaprogramming utilities
  //------------------------------------------------------------------------
  /// Trait type evaluating true if `T` is derived from `PartitionBase`.
  template <typename T, typename = void>
  struct is_partition_type : public std::false_type {};

  template <typename Part>
  struct is_partition_type<
    Part,
    std::enable_if_t<std::is_base_of<PartitionBase, std::decay_t<Part>>::value>>
    : public std::true_type {};

  /// Constant true if `T` is derived from `PartitionBase`.
  template <typename T>
  constexpr bool is_partition_type_v = is_partition_type<T>();

  //------------------------------------------------------------------------
  /// Trait type evaluating true if `T` is pointer to some `PartitionBase`.
  template <typename, typename = void>
  struct is_partition_ptr : public std::false_type {};

  template <typename PartPtr>
  struct is_partition_ptr<PartPtr,
                          std::enable_if_t<is_partition_type_v<decltype(*std::declval<PartPtr>())>>>
    : public std::true_type {};

  /// Constant true if `T` is pointer to some `PartitionBase`.
  template <typename T>
  constexpr bool is_partition_ptr_v = is_partition_ptr<T>();

  //------------------------------------------------------------------------
  /// Trait type evaluating true if `T` is iterator to some `PartitionBase`.
  template <typename, typename = void>
  struct is_partition_ptr_iterator : public std::false_type {};

  template <typename Iter>
  struct is_partition_ptr_iterator<
    Iter,
    std::enable_if_t<is_partition_ptr_v<std::decay_t<typename Iter::value_type>>>>
    : public std::true_type {};

  //------------------------------------------------------------------------
  //---  Sorting objects
  //------------------------------------------------------------------------
  /// Class extracting the lower bound of the specified range of an area.
  template <AreaOwner::AreaRangeMember_t Range>
  struct RangeLowerBoundExtractor {
    static constexpr auto range = Range;
    using Area_t = AreaOwner::Area_t;

    double operator()(double lower) const { return lower; }
    double operator()(Area_t::Range_t const& r) const { return (*this)(r.lower); }
    double operator()(Area_t const& area) const { return (*this)(area.*range); }
    double operator()(AreaOwner const& area) const { return (*this)(area.area()); }
    double operator()(AreaOwner const* ptr) const { return (*this)(*ptr); }
  };

  //------------------------------------------------------------------------
  /// Class extracting the lower bound of the specified range of a partition area.
  template <PartitionBase::AreaRangeMember_t Range>
  struct PartitionRangeLowerBoundExtractor : public RangeLowerBoundExtractor<Range> {
    using Base_t = RangeLowerBoundExtractor<Range>;
    using Area_t = typename Base_t::Area_t;

    using Partition_t = PartitionBase; ///< Base type of the partition.

    using Base_t::operator(); // import inherited versions

    auto operator()(Partition_t const& part) { return Base_t::operator()(part.area()); }
    template <typename PartPtr,
              typename = std::enable_if_t<details::is_partition_ptr<PartPtr>::value>>
    auto operator()(PartPtr const& part)
    {
      return operator()(*part);
    }
  };

  //------------------------------------------------------------------------
  template <PartitionBase::AreaRangeMember_t Range>
  struct PartitionSorterByAreaRangeLower {
    using Sorter_t = PartitionSorterByAreaRangeLower<Range>;

    using KeyExtractor_t = PartitionRangeLowerBoundExtractor<Range>;

    static constexpr auto range = Range;

    template <typename T>
    static auto key(T const& obj)
    {
      return KeyExtractor_t()(obj);
    }

    /// Type of sorting key. In short: `double`.
    using Key_t = decltype(Sorter_t::key(std::declval<PartitionBase>()));

    static Key_t key(Key_t k) { return k; } // shortcut
    static bool sortKey(Key_t a, Key_t b) { return a < b; }

    template <typename A, typename B>
    bool operator()(A const& a, B const& b) const
    {
      return sortKey(key(a), key(b));
    }
  };

  //------------------------------------------------------------------------

} // namespace geo::part::details

//------------------------------------------------------------------------------
//---  partition data description
//---

template <typename Data>
template <typename Stream>
geo::part::PartitionDataDescriber<Data>::PartitionDataDescriber(Stream&& out,
                                                                Data const* data,
                                                                std::string indent /* = "" */,
                                                                std::string firstIndent /* = "" */
)
{
  out << firstIndent;
  std::string typeName = lar::debug::demangle<Data>();
  if (data) { out << typeName << "[" << ((void*)data) << "]"; }
  else {
    out << "no '" << typeName << "' data";
  }
}

template <typename Stream, typename Data>
void geo::part::describePartitionData(Stream&& out,
                                      Data const* data,
                                      std::string indent /* = "" */,
                                      std::string firstIndent /* = "" */
)
{
  PartitionDataDescriber<Data>(std::forward<Stream>(out), data, indent, firstIndent);
}

//------------------------------------------------------------------------------
//--- geo::part::Partition
//---
template <typename Data>
typename geo::part::Partition<Data>::Subpartitions_t const geo::part::Partition<Data>::NoSubparts;

//------------------------------------------------------------------------------
template <typename Data>
template <typename Pred>
void geo::part::Partition<Data>::walk(Partition_t const* start, Pred&& pred)
{
  if (!start) return;
  pred(*start);

  // recursive implementation
  for (auto const& subPart : start->parts())
    subPart->walk(std::forward<Pred>(pred));
}

//------------------------------------------------------------------------------
//---  geo::part::PartitionWithData
//---
template <typename Data>
std::string geo::part::PartitionWithData<Data>::doDescribe(std::string indent,
                                                           std::string firstIndent) const
{
  std::string msg = Base_t::doDescribe(indent, firstIndent);
  if (data()) {
    std::ostringstream sstr;
    sstr << ": ";
    describePartitionData(sstr, data(), indent);
    msg += sstr.str();
  }
  else {
    msg += " (no data)";
  }
  return msg;
}

//------------------------------------------------------------------------------
//--- geo::part::PartitionContainer
//---
//------------------------------------------------------------------------------
template <typename Data>
auto geo::part::PartitionContainer<Data>::atPoint(double w, double d) const -> Data_t*
{
  if (!Base_t::contains(w, d)) return nullptr; // not our point at all
  // it's ours; see if it belongs to a subpart
  auto part = findPart(w, d);
  return part ? part->atPoint(w, d) : Base_t::data();
}

//------------------------------------------------------------------------------
template <typename Data>
std::string geo::part::PartitionContainer<Data>::doDescribe(std::string indent,
                                                            std::string firstIndent) const
{
  std::string msg = firstIndent + describeIntro();
  if (Base_t::data()) {
    std::ostringstream sstr;
    sstr << ", and ";
    describePartitionData(sstr, Base_t::data(), indent, "");
    msg += sstr.str();
  }

  for (auto const& part : parts()) {
    msg += "\n" + indent + " * ";
    msg += part->describe(indent + "  ", "");
  }

  return msg;
}

//------------------------------------------------------------------------------
template <typename Data>
std::string geo::part::PartitionContainer<Data>::describeIntro() const
{
  return std::to_string(parts().size()) + " subpartitions";
}

//------------------------------------------------------------------------------
//---  geo::part::SortedPartition
//---
template <typename Data, typename Sorter>
auto geo::part::SortedPartition<Data, Sorter>::findPartWithKey(double key) const
  -> Partition_t const*
{
  auto pbegin = Base_t::parts().cbegin();
  auto iPart = std::upper_bound(pbegin, Base_t::parts().cend(), key, sorter);
  return (iPart == pbegin) ? nullptr : (--iPart)->get();
}

//------------------------------------------------------------------------------
template <typename Data, typename Sorter>
void geo::part::SortedPartition<Data, Sorter>::initParts()
{
  /*
   * Initialization tasks:
   * - ensure that the parts are sorted by increasing depth
    */
  std::sort(Base_t::myParts.begin(), Base_t::myParts.end(), sorter);
}

//------------------------------------------------------------------------------
//---  geo::part::DepthPartition
//---
template <typename Data>
std::string geo::part::DepthPartition<Data>::describeIntro() const
{
  std::ostringstream sstr;
  sstr << Base_t::size() << " partitions along depth covering " << Base_t::area();
  return sstr.str();
}

//------------------------------------------------------------------------------
//---  geo::part::WidthPartition
//---
template <typename Data>
std::string geo::part::WidthPartition<Data>::describeIntro() const
{
  std::ostringstream sstr;
  sstr << Base_t::size() << " partitions along width covering " << Base_t::area();
  return sstr.str();
}

//------------------------------------------------------------------------------
//---  geo::part::GridPartition
//---
template <typename Data>
geo::part::GridPartition<Data>::GridPartition(Area_t const& area,
                                              Subpartitions_t&& subpartitions,
                                              unsigned int nDepthPartitions,
                                              unsigned int nWidthPartitions,
                                              Data_t* defData /* = nullptr */
                                              )
  : Base_t(area, std::move(subpartitions), defData)
  , widthSeps(computeWidthSeps(nDepthPartitions, nWidthPartitions))
  , depthSeps(computeDepthSeps(nDepthPartitions, nWidthPartitions))
{
  assert(nWidthPartitions * nDepthPartitions == Base_t::size());
}

template <typename Data>
geo::part::GridPartition<Data>::GridPartition(Area_t const& area,
                                              Subpartitions_t&& subpartitions,
                                              unsigned int nDepthPartitions,
                                              Data_t* defData /* = nullptr */
                                              )
  : GridPartition(area,
                  std::move(subpartitions),
                  nDepthPartitions,
                  (nDepthPartitions ? subpartitions.size() / nDepthPartitions : 0),
                  defData)
{}

//------------------------------------------------------------------------------
template <typename Data>
std::string geo::part::GridPartition<Data>::doDescribe(std::string indent,
                                                       std::string firstIndent) const
{
  std::ostringstream sstr;
  sstr << firstIndent << Base_t::describeIntro() << " in a (WxD) = " << nWidthParts() << " x "
       << nDepthParts() << " grid";
  if (Base_t::data()) {
    sstr << ", and ";
    describePartitionData(sstr, Base_t::data(), indent, "");
  }
  for (std::size_t iDepth = 0; iDepth < nDepthParts(); ++iDepth) {
    for (std::size_t iWidth = 0; iWidth < nWidthParts(); ++iWidth) {
      sstr << "\n"
           << indent << " [" << iDepth << "][" << iWidth << "] "
           << part(iDepth, iWidth)->describe(indent + "  ", "");
    } // for width
  }   // for depth
  return sstr.str();
}

//------------------------------------------------------------------------------
template <typename Data>
auto geo::part::GridPartition<Data>::findPart(double w, double d) const -> Partition_t const*
{
  auto const iWidth = std::upper_bound(widthSeps.cbegin(), widthSeps.cend(), w);
  if (iWidth == widthSeps.cbegin()) return nullptr;
  auto const iDepth = std::upper_bound(depthSeps.cbegin(), depthSeps.cend(), d);
  if (iDepth == depthSeps.cbegin()) return nullptr;
  return part(std::distance(depthSeps.cbegin(), iDepth) - 1U,
              std::distance(widthSeps.cbegin(), iWidth) - 1U)
    .get();
}

//------------------------------------------------------------------------------
template <typename Data>
std::vector<double> geo::part::GridPartition<Data>::computeWidthSeps(unsigned int,
                                                                     unsigned int nW) const
{
  return detectSeparators<&Area_t::width>(
    Base_t::parts().cbegin(), Base_t::parts().cend(), nW, 1U, nW);
}

template <typename Data>
std::vector<double> geo::part::GridPartition<Data>::computeDepthSeps(unsigned int nD,
                                                                     unsigned int nW) const
{
  return detectSeparators<&Area_t::depth>(
    Base_t::parts().cbegin(), Base_t::parts().cend(), nD, nW, 1U);
}

template <typename Data>
template <geo::part::AreaOwner::AreaRangeMember_t Range, typename BeginIter, typename EndIter>
std::vector<double> geo::part::GridPartition<Data>::detectSeparators(BeginIter b,
                                                                     EndIter e,
                                                                     std::size_t const nGroups,
                                                                     std::size_t const startDelta,
                                                                     std::size_t const stride)
{
  /*
   * The iterators are better be random access.
   * The range [b,e[ is considered to be a 2D table stored row after row.
   * This function can operate on rows or columns, given the proper arguments.
   *
   * The full range is split in nGroups "groups" (e.g. rows or columns).
   * Each group g starts at the element g x startDelta.
   * All members of that group are stride elements far from each other.
   *
   * Be the table of size (d x w).
   * To process data row by row:
   * - the group is a row
   * - the number of groups is the number of rows: nGroups = d
   * - each group will have (d x w)/nGroups = w elements
   * - the start offset is the number of group times the size of it:
   *   startDelta = w
   * - the elements are contiguous: stride = 1
   *
   * To process data column by column:
   * - the group is a column
   * - the number of groups is the number of columns: nGroups = w
   * - each group will have (d x w)/nGroups = d elements
   * - the start offset matches the number of group: startDelta = 1
   * - the elements are separated by a full row: stride = w
   *
   */
  // the separator is on the lower bound of selected range of the partition area

  static_assert(details::is_partition_ptr_iterator<BeginIter>(),
                "Begin iterator does not point to a pointer to partition type");

  details::PartitionRangeLowerBoundExtractor<Range> lowerBound;

  std::size_t const nParts = std::distance(b, e);
  std::size_t const nPartsInGroup = nParts / nGroups;

  auto const part = [b](std::size_t index) { return std::next(b, index)->get(); };

  std::vector<double> seps(nGroups);
  for (size_t g = 0; g < nGroups; ++g) {

    double& sep = seps[g];

    // indices of an element in the previous and next group, respectively
    std::size_t index = g * startDelta;
    sep = lowerBound(part(index));

    std::size_t const iend = index + nPartsInGroup * stride;
    while ((index += stride) < iend) {
      double const l = lowerBound(part(index));
      if (sep > l) sep = l;
    } // while (element)

  } // while (groups)
  return seps;
}

//******************************************************************************

#endif // LARCOREALG_GEOMETRY_PARTITIONS_H
