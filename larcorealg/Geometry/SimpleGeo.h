/**
 * @file   SimpleGeo.h
 * @brief  Some simple functions to represent geometry entities
 * @ingroup Geometry
 *
 * This library is simple and header-only.
 *
 */

#ifndef LARCOREALG_GEOMETRY_SIMPLEGEO_H
#define LARCOREALG_GEOMETRY_SIMPLEGEO_H

// C/C++ standard library
#include <algorithm> // std::min(), std::max()
#include <array>
#include <stdexcept> // std::runtime_error

namespace lar {
  namespace util {

    // --- BEGIN Geometry group ------------------------------------------------
    /// @ingroup Geometry
    /// @{

    /**
     * @brief Simple class definitions for geometry concepts.
     *
     * This namespace provides data structures for "points" and extended
     * structures.
     *
     * Points are either 2D (`Point2D`) or 3D (`Point3D`).
     *
     * The extended structures are offered in two flavours:
     * * point-based (area, volume): definitions of the boundaries are
     *   provided as two points of proper dimensionality: `Point2D` for `Area`,
     *   `pOint3D` for `Volume`
     * * dimension-based (area only so far): definitions of the boundaries are
     *   provided as one range (`Range`) per dimension
     *
     *
     */
    namespace simple_geo {

      /// @{
      /// @name Dimensionless objects (points)

      /// 2D point (x, y) (by default, with double precision)
      template <typename Data = double>
      struct Point2D {
        using Data_t = Data;
        Data_t x = 0.;
        Data_t y = 0.;

        Point2D() = default;
        Point2D(Data_t x, Data_t y) : x(x), y(y) {}
        Point2D(Data_t const* p) : x(p[0]), y(p[1]) {}

      }; // struct Point2D<>

      template <typename T>
      bool operator==(Point2D<T> const& a, Point2D<T> const& b)
      {
        return (a.x == b.x) && (a.y == b.y);
      }
      template <typename T>
      bool operator!=(Point2D<T> const& a, Point2D<T> const& b)
      {
        return (a.x != b.x) || (a.y != b.y);
      }
      template <typename T>
      Point2D<T> operator+(Point2D<T> const& a, Point2D<T> const& b)
      {
        return {a.x + b.x, a.y + b.y};
      }
      template <typename T>
      Point2D<T> operator*(Point2D<T> const& p, typename Point2D<T>::Data_t f)
      {
        return {p.x * f, p.y * f};
      }
      template <typename T>
      Point2D<T> operator/(Point2D<T> const& p, typename Point2D<T>::Data_t f)
      {
        return {p.x / f, p.y / f};
      }
      template <typename Stream, typename T>
      Stream& operator<<(Stream&& out, Point2D<T> const& p)
      {
        out << "( " << p.x << " ; " << p.y << " )";
        return out;
      }

      /// 3D point (x, y, z) (by default, with double precision)
      template <typename Data = double>
      struct Point3D {
        using Data_t = Data;
        Data_t x = 0.;
        Data_t y = 0.;
        Data_t z = 0.;

        Point3D() = default;
        Point3D(Data_t x, Data_t y, Data_t z) : x(x), y(y), z(z) {}
        Point3D(Data_t const* p) : x(p[0]), y(p[1]), z(p[2]) {}
      }; // struct Point3D<>

      template <typename T>
      bool operator==(Point3D<T> const& a, Point3D<T> const& b)
      {
        return (a.x == b.x) && (a.y == b.y) && (a.z == b.z);
      }
      template <typename T>
      bool operator!=(Point3D<T> const& a, Point3D<T> const& b)
      {
        return (a.x != b.x) || (a.y != b.y) || (a.z != b.z);
      }
      template <typename T>
      Point3D<T> operator+(Point3D<T> const& a, Point3D<T> const& b)
      {
        return {a.x + b.x, a.y + b.y, a.z + b.z};
      }
      template <typename T>
      Point3D<T> operator*(Point3D<T> const& p, typename Point2D<T>::Data_t f)
      {
        return {p.x * f, p.y * f, p.z * f};
      }
      template <typename T>
      Point3D<T> operator/(Point3D<T> const& p, typename Point2D<T>::Data_t f)
      {
        return {p.x / f, p.y / f, p.z / f};
      }
      template <typename Stream, typename T>
      Stream& operator<<(Stream&& out, Point3D<T> const& p)
      {
        out << "( " << p.x << " ; " << p.y << " ; " << p.z << " )";
        return out;
      }

      /// @}

      /// @{
      /// @name Point-bounded volumes

      /// Area/volume delimited by points: base/1D implementation
      template <typename Point>
      class AreaBase {
      public:
        using Point_t = Point;
        using Area_t = AreaBase<Point>;

        /// Exception thrown when result of intersection is null
        class NullIntersection : public std::runtime_error {
          using std::runtime_error::runtime_error;
        };

        AreaBase() = default;

        AreaBase(Point_t const& a, Point_t const& b) { set_sorted(min.x, max.x, a.x, b.x); }

        Point_t const& Min() const { return min; }
        Point_t const& Max() const { return max; }
        Point_t Center() const { return (min + max) * 0.5; }

        auto DeltaX() const { return Max().x - Min().x; }
        bool isEmptyX() const { return (DeltaX() == 0); }
        bool isNullX() const { return Max().x < Min().x; }
        bool isNull() const { return isNullX(); }
        unsigned int nonEmptyDims() const { return isEmptyX() ? 0 : 1; }
        bool isEmpty() const { return nonEmptyDims() == 0; }
        bool isLine() const { return nonEmptyDims() == 1; }

        auto thinnestSize() const { return DeltaX(); }

        /// Returns the index of the thinnest side (0 is x)
        unsigned int thinnestSide() const { return 0; }

        void IncludePoint(Point_t const& point) { set_min_max(min.x, max.x, point.x); }

        void Include(Area_t const& area)
        {
          IncludePoint(area.min);
          IncludePoint(area.max);
        }

        void Intersect(Area_t const& area)
        {
          set_max(min.x, area.min.x);
          set_min(max.x, area.max.x);
          if (isNullX()) throw NullIntersection("null x dimension");
        }

        bool operator==(Area_t const& as) const { return (min == as.min) && (max == as.max); }
        bool operator!=(Area_t const& as) const { return (min != as.min) || (max != as.max); }

      protected:
        using Data_t = typename Point_t::Data_t;

        Point_t min, max;

        static void set_min(Data_t& var, Data_t val)
        {
          if (val < var) var = val;
        }

        static void set_max(Data_t& var, Data_t val)
        {
          if (val > var) var = val;
        }

        static void set_min_max(Data_t& min_var, Data_t& max_var, Data_t val)
        {
          set_min(min_var, val);
          set_max(max_var, val);
        }

        static void set_sorted(Data_t& min_var, Data_t& max_var, Data_t a, Data_t b)
        {
          if (a > b) {
            min_var = b;
            max_var = a;
          }
          else {
            min_var = a;
            max_var = b;
          }
        }

      }; // class AreaBase<>

      /// Area delimited by two points
      template <typename Point = Point2D<double>>
      class Area : public AreaBase<Point> {
        using Base_t = AreaBase<Point>;

      public:
        using Point_t = typename Base_t::Point_t;
        using Area_t = Area<Point_t>;

        Area() = default;

        Area(Point_t const& a, Point_t const& b) : Base_t(a, b)
        {
          Base_t::set_sorted(Base_t::min.y, Base_t::max.y, a.y, b.y);
        }

        auto DeltaY() const { return Base_t::Max().y - Base_t::Min().y; }
        bool isEmptyY() const { return DeltaY() == 0; }
        unsigned int nonEmptyDims() const { return Base_t::nonEmptyDims() + (isEmptyY() ? 0 : 1); }
        bool isNullY() const { return Base_t::Max().y < Base_t::Min().y; }
        bool isNull() const { return Base_t::isNull() || isNullY(); }
        bool isEmpty() const { return nonEmptyDims() == 0; }
        bool isLine() const { return nonEmptyDims() == 1; }
        bool isPlane() const { return nonEmptyDims() == 2; }

        auto thinnestSize() const { return std::min(DeltaY(), Base_t::thinnestSize()); }

        /// Returns the index of the thinnest side (0 is x, 1 is y)
        unsigned int thinnestSide() const
        {
          return (DeltaY() < Base_t::thinnestSize()) ? 1 : Base_t::thinnestSide();
        }

        void IncludePoint(Point_t const& point)
        {
          Base_t::IncludePoint(point);
          Base_t::set_min_max(Base_t::min.y, Base_t::max.y, point.y);
        }

        void Include(Area_t const& area)
        {
          IncludePoint(area.min);
          IncludePoint(area.max);
        }

        void Intersect(Area_t const& area)
        {
          Base_t::Intersect(area);
          Base_t::set_max(Base_t::min.y, area.min.y);
          Base_t::set_min(Base_t::max.y, area.max.y);
          if (isNullY()) throw typename Base_t::NullIntersection("null y dimension");
        }

      }; // class Area<>

      template <typename Stream, typename Point>
      Stream& operator<<(Stream&& out, AreaBase<Point> const& area)
      {
        out << area.Min() << " - " << area.Max();
        return out;
      }

      /// Volume delimited by two points
      template <typename Point = Point3D<double>>
      class Volume : public Area<Point> {
        using Base_t = Area<Point>;

      public:
        using Point_t = typename Base_t::Point_t;
        using Volume_t = Volume<Point_t>;

        Volume() = default;

        Volume(Point_t const& a, Point_t const& b) : Base_t(a, b)
        {
          Base_t::set_sorted(Base_t::min.z, Base_t::max.z, a.z, b.z);
        }

        auto DeltaZ() const { return Base_t::Max().z - Base_t::Min().z; }
        bool isEmptyZ() const { return DeltaZ() == 0; }
        unsigned int nonEmptyDims() const { return Base_t::nonEmptyDims() + (isEmptyZ() ? 0 : 1); }
        bool isNullZ() const { return Base_t::Max().z < Base_t::Min().z; }
        bool isNull() const { return Base_t::isNull() || isNullZ(); }
        bool isEmpty() const { return nonEmptyDims() == 0; }
        bool isLine() const { return nonEmptyDims() == 1; }
        bool isPlane() const { return nonEmptyDims() == 2; }
        bool isVolume() const { return nonEmptyDims() == 3; }

        auto thinnestSize() const { return std::min(DeltaZ(), Base_t::thinnestSize()); }

        /// Returns the index of the thinnest side (0 is x, 1 is y)
        unsigned int thinnestSide() const
        {
          return (DeltaZ() < Base_t::thinnestSize()) ? 2 : Base_t::thinnestSide();
        }

        void IncludePoint(Point_t const& point)
        {
          Base_t::IncludePoint(point);
          Base_t::set_min_max(Base_t::min.z, Base_t::max.z, point.z);
        }

        void Include(Volume_t const& volume)
        {
          IncludePoint(volume.min);
          IncludePoint(volume.max);
        }

        void Intersect(Volume_t const& volume)
        {
          Base_t::Intersect(volume);
          Base_t::set_max(Base_t::min.z, volume.min.z);
          Base_t::set_min(Base_t::max.z, volume.max.z);
          if (isNullZ()) throw typename Base_t::NullIntersection("null z dimension");
        }

      }; // class Volume<>

      /// @}

      /// @{
      /// @name Dimension-bounded volumes

      /// Definition of a range along one dimension.
      template <typename Data = double>
      struct Range {
        using Data_t = Data;         ///< Numeric type for boundaries.
        using Range_t = Range<Data>; ///< This range type.

        Data_t lower = 1.0; ///< Starting coordinate.
        Data_t upper = 0.0; ///< Ending coordinate.

        /// Default constructor: empty range.
        Range() = default;

        /// Constructor from lower and upper bounds.
        Range(Data_t lower, Data_t upper, bool doSort = false) : lower(lower), upper(upper)
        {
          if (doSort) sort();
        }

        /// Returns whether the range is empty.
        bool isNull() const { return lower >= upper; }

        /// Returns the distance between upper and lower bounds.
        Data_t length() const { return std::max(upper - lower, Data_t(0.0)); }

        /// Returns whether the specified value is within the range.
        bool contains(Data_t v) const { return (v >= lower) && (v <= upper); }

        /// Returns whether the specified range overlaps this range.
        bool overlaps(Range_t const& r) const;

        /// Returns a value that, added to v, makes it fall within a margin in
        /// the range.
        Data_t delta(Data_t v, Data_t margin = 0.0) const;

        /// Extends the range to include the specified point.
        void extendToInclude(Data_t);

        /// Extends the range to include the specified point.
        void extendToInclude(Range_t const& r);

        /// Shortens this range to include only points also in `r`.
        void intersect(Range_t const& r);

      private:
        /// Ensures order of boundaries. Corrupts invalid ranges.
        void sort()
        {
          if (lower > upper) std::swap(lower, upper);
        }

        /// Resets this range to be empty (that is, like default-constructed).
        void makeNull() { *this = Range_t{}; }

      }; // struct Range<>

      /// Prints the specified range to a stream: "( lower -- upper )".
      template <typename Stream, typename Data>
      Stream& operator<<(Stream&& out, Range<Data> const& range);

      /**
       * @brief Definition of a rectangle from dimension ranges.
       * @tparam Data numerical type for boundary coordinates
       * @see Range, Area
       *
       * This object defines a 2D area (rectangle) as a list of one range for
       * each dimension. Dimensions are called "width" and "depth".
       *
       * If the use case asks for point-driven area rather than a
       * dimension-driven area, use `Area` instead.
       *
       */
      template <typename Data = double>
      struct Rectangle {
        using Data_t = Data;                 ///< Numerical type for boundaries.
        using Rectangle_t = Rectangle<Data>; ///< This type.
        using Range_t = Range<Data_t>;       ///< Type for dimension boundaries.

        Range_t width; ///< Range along width direction.
        Range_t depth; ///< Range along depth direction.

        /// Default constructor: an empty rectangle.
        Rectangle() = default;

        /// Constructor from width and depth ranges.
        Rectangle(Range_t const& width, Range_t const& depth) : width(width), depth(depth) {}

        /// Returns whether the rectangle has null area.
        bool isNull() const { return width.isNull() || depth.isNull(); }

        /// Returns whether the specified point is in the area.
        bool contains(Data_t w, Data_t d) const { return width.contains(w) && depth.contains(d); }

        /// Returns whether this and the specified rectangle overlap.
        bool overlaps(Rectangle_t const& r) const;

        /// Extends the range to include the specified point.
        void extendToInclude(Rectangle_t const& r);

      }; // Rectangle<>

      /// Prints the specified rectangle to a stream: "w=Wrange d=Drange".
      template <typename Stream, typename Data>
      Stream& operator<<(Stream&& out, Rectangle<Data> const& rect);

      /// @}

    } // namespace simple_geo

    /// @}
    // --- END Geometry group --------------------------------------------------

  } // namespace util
} // namespace lar

//==============================================================================
//--- Template implementation
//---
//------------------------------------------------------------------------------
//---  geo::PlaneGeo::Range<>
//---
template <typename Data>
auto lar::util::simple_geo::Range<Data>::delta(Data_t v, Data_t margin /* = 0 */) const -> Data_t
{

  if (v < (lower + margin)) return lower + margin - v; // always positive
  if (v > (upper - margin)) return upper - margin - v; // always negative
  return 0.0;                                          // always zero

} // lar::util::simple_geo::Range<Data>::delta()

//------------------------------------------------------------------------------
template <typename Data>
void lar::util::simple_geo::Range<Data>::extendToInclude(Data_t v)
{

  if (lower > upper)
    lower = upper = v;
  else if (lower > v)
    lower = v;
  else if (upper < v)
    upper = v;

} // lar::util::simple_geo::Range<Data>::extendToInclude()

//------------------------------------------------------------------------------
template <typename Data>
void lar::util::simple_geo::Range<Data>::extendToInclude(Range_t const& r)
{

  if (r.isNull()) return;
  if (isNull()) {
    *this = r;
    return;
  }
  extendToInclude(r.lower);
  extendToInclude(r.upper);

} // lar::util::simple_geo::Range<Data>::extendToInclude()

//------------------------------------------------------------------------------
template <typename Data>
void lar::util::simple_geo::Range<Data>::intersect(Range_t const& r)
{
  // this implementation explicitly makes the range default-constructed-null
  // if the intersection results empty
  if (isNull()) return;
  if (r.isNull()) {
    makeNull();
    return;
  }
  if (lower < r.lower) lower = r.lower;
  if (upper > r.upper) upper = r.upper;
  if (lower > upper) makeNull();
} // lar::util::simple_geo::Range<Data>::intersect()

//------------------------------------------------------------------------------
template <typename Data>
bool lar::util::simple_geo::Range<Data>::overlaps(Range_t const& r) const
{
  if (isNull() || r.isNull()) return false;
  return (r.lower < upper) && (lower < r.upper);
} // lar::util::simple_geo::Range<Data>::overlaps()

//------------------------------------------------------------------------------
template <typename Stream, typename Data>
Stream& lar::util::simple_geo::operator<<(Stream&& out, Range<Data> const& range)
{
  out << "( " << range.lower << " -- " << range.upper << " )";
  return out;
} // operator<< (Range)

//------------------------------------------------------------------------------
template <typename Data>
void lar::util::simple_geo::Rectangle<Data>::extendToInclude(Rectangle_t const& r)
{
  width.extendToInclude(r.width);
  depth.extendToInclude(r.depth);
} // lar::util::simple_geo::Rectangle<Data>::extendToInclude()

//------------------------------------------------------------------------------
template <typename Data>
bool lar::util::simple_geo::Rectangle<Data>::overlaps(Rectangle_t const& r) const
{
  if (isNull() || r.isNull()) return false;
  return width.overlap(r.width) && depth.overlap(r.depth);
} // lar::util::simple_geo::Rectangle<Data>::overlaps()

//------------------------------------------------------------------------------
template <typename Stream, typename Data>
Stream& lar::util::simple_geo::operator<<(Stream&& out, Rectangle<Data> const& rect)
{
  out << "w=" << rect.width << " d=" << rect.depth;
  return out;
} // operator<< (Rectangle)

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_SIMPLEGEO_H
