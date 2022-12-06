////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/TPCGeo.h
/// \brief Encapsulate the construction of a single detector plane
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARCOREALG_GEOMETRY_TPCGEO_H
#define LARCOREALG_GEOMETRY_TPCGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT libraries
#include "TGeoMatrix.h"
#include "TGeoVolume.h"

// C/C++ standard library
#include <set>
#include <vector>

class TGeoNode;

namespace geo {

  //......................................................................
  /// @brief Geometry information for a single TPC.
  /// @ingroup Geometry
  class TPCGeo : public BoxBoundedGeo {
  public:
    using ID_t = TPCID;
    using PlaneCollection_t = std::vector<geo::PlaneGeo>;
    using GeoNodePath_t = geo::WireGeo::GeoNodePath_t;

    /// Type returned by `IterateElements()`.
    using ElementIteratorBox = PlaneCollection_t const&;

    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     *
     * These types represents points and displacement vectors in the reference
     * frame defined in the TPC geometry box from the GDML geometry
     * description.
     *
     * No alias is explicitly defined for the LArSoft global vector types,
     * `geo::Point_t` and `geo::Vector_t`.
     *
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different
     * instances of `geo::TPCGeo` have the same type but are not compatible.
     */

    /// Tag for vectors in the "local" GDML coordinate frame of the TPC.
    struct TPCGeoCoordinatesTag {};

    /// Type of points in the local GDML TPC frame.
    using LocalPoint_t = geo::Point3DBase_t<TPCGeoCoordinatesTag>;

    /// Type of displacement vectors in the local GDML TPC frame.
    using LocalVector_t = geo::Vector3DBase_t<TPCGeoCoordinatesTag>;

    ///@}

    /// Type of 2D vector projection on a plane.
    using Projection_t = geo::PlaneGeo::WidthDepthProjection_t;

    /// Data structure with plane and drift projections of a 3D vector.
    using DecomposedVector_t = geo::PlaneGeo::WDDecomposedVector_t;

    // Construct a representation of a single plane of the detector
    TPCGeo(TGeoNode const* node, geo::TransformationMatrix&& trans, PlaneCollection_t&& planes);

    /// @{
    /// @name TPC properties

    /// Half width (associated with x coordinate) of active TPC volume [cm].
    double ActiveHalfWidth() const { return fActiveHalfWidth; }
    /// Width (associated with x coordinate) of active TPC volume [cm].
    double ActiveWidth() const { return 2.0 * ActiveHalfWidth(); }
    /// Half height (associated with y coordinate) of active TPC volume [cm].
    double ActiveHalfHeight() const { return fActiveHalfHeight; }
    /// Height (associated with y coordinate) of active TPC volume [cm].
    double ActiveHeight() const { return 2.0 * ActiveHalfHeight(); }
    /// Length (associated with z coordinate) of active TPC volume [cm].
    double ActiveLength() const { return fActiveLength; }
    /// Length (associated with z coordinate) of active TPC volume [cm].
    double ActiveHalfLength() const { return fActiveLength / 2.0; }
    /// Width is associated with x coordinate [cm].
    double HalfWidth() const { return fHalfWidth; }
    /// Width is associated with x coordinate [cm].
    double Width() const { return 2.0 * HalfWidth(); }
    /// Height is associated with y coordinate [cm].
    double HalfHeight() const { return fHalfHeight; }
    /// Height is associated with y coordinate [cm].
    double Height() const { return 2.0 * HalfHeight(); }
    /// Length is associated with z coordinate [cm].
    double Length() const { return fLength; }
    /// Length is associated with z coordinate [cm].
    double HalfLength() const { return fLength / 2.0; }
    double ActiveMass() const { return fActiveVolume->Weight(); }
    const TGeoVolume* ActiveVolume() const { return fActiveVolume; }
    const TGeoVolume* TotalVolume() const { return fTotalVolume; }

    /// Returns the direction `Width()` is measured on.
    decltype(auto) WidthDir() const { return fWidthDir; }

    /// Returns the direction `Height()` is measured on.
    decltype(auto) HeightDir() const { return fHeightDir; }

    /// Returns the direction `Length()` is measured on.
    decltype(auto) LengthDir() const { return fLengthDir; }

    /// Returns an enumerator value describing the drift direction.
    DriftDirection_t DriftDirection() const { return fDriftDirection; }

    /// Returns the direction of the drift (vector pointing toward the planes).
    Vector_t DriftDir() const { return fDriftDir; }

    /// Drift distance is defined as the distance between the last anode plane
    /// and the opposite face of the TPC, in centimeters.
    double DriftDistance() const { return ComputeDriftDistance(); }

    /// @}

    /// @{
    /// @name Plane access

    //@{
    /// Number of planes in this tpc
    unsigned int Nplanes() const { return fPlanes.size(); }
    unsigned int NElements() const { return fPlanes.size(); }
    //@}

    //@{
    /**
     * @brief Returns whether a plane with index iplane is present in this TPC
     * @param iplane index of plane in this TPC
     * @return whether the plane with index iplane is present in this TPC
     */
    bool HasPlane(unsigned int iplane) const { return iplane < Nplanes(); }
    bool HasElement(unsigned int iplane) const { return HasPlane(iplane); }
    //@}

    //@{
    /**
     * @brief Returns whether the plane in planeid is present in this TPC
     * @param planeid full plane ID
     * @return whether the plane in planeid is present in this TPC
     *
     * The cryostat and TPC numbers in planeid are ignored, as it is ignored
     * whether planeid is invalid.
     */
    bool HasPlane(geo::PlaneID const& planeid) const { return HasPlane(planeid.Plane); }
    bool HasElement(geo::PlaneID const& planeid) const { return HasPlane(planeid); }
    //@}

    /// Return the plane in the tpc with View_t view.
    PlaneGeo const& Plane(geo::View_t view) const;

    /// Return the iplane'th plane in the TPC.
    /// @throws cet::exception (category "PlaneOutOfRange")  if no such plane
    PlaneGeo const& Plane(unsigned int iplane) const;

    //@{
    /**
     * @brief Returns the plane in planeid from this TPC
     * @param planeid full plane ID
     * @return a constant reference to the plane in planeid
     * @throws cet::exception (category "PlaneOutOfRange") if no such plane
     *
     * The cryostat and TPC numbers in planeid are ignored, as it is ignored
     * whether planeid is invalid.
     */
    const PlaneGeo& Plane(PlaneID const& planeid) const { return Plane(planeid.Plane); }
    const PlaneGeo& GetElement(PlaneID const& planeid) const { return Plane(planeid); }
    //@}

    /**
     * @brief Returns the plane number iplane from this TPC
     * @param iplane the number of local plane
     * @return a constant pointer to the plane, or nullptr if it does not exist
     */
    PlaneGeo const* PlanePtr(unsigned int iplane) const
    {
      return HasPlane(iplane) ? &(fPlanes[iplane]) : nullptr;
    }

    //@{
    /**
     * @brief Returns the plane in planeid from this TPC
     * @param planeid full plane ID
     * @return a constant pointer to the plane, or nullptr if it does not exist
     *
     * The cryostat and TPC numbers in planeid are ignored, as it is ignored
     * whether planeid is invalid.
     */
    PlaneGeo const* PlanePtr(PlaneID const& planeid) const { return PlanePtr(planeid.Plane); }
    PlaneGeo const* GetElementPtr(PlaneID const& planeid) const { return PlanePtr(planeid); }
    //@}

    /// Returns the wire plane with the smallest surface
    geo::PlaneGeo const& SmallestPlane() const;

    /// Returns the first wire plane (the closest to TPC center).
    geo::PlaneGeo const& FirstPlane() const { return fPlanes[0]; }

    /// Returns the last wire plane (the farther from TPC center).
    geo::PlaneGeo const& LastPlane() const { return fPlanes[Nplanes() - 1]; }

    /// @brief Returns the largest number of wires among the planes in this TPC
    unsigned int MaxWires() const;

    // @{
    /**
     * @brief Returns an object for iterating through all `geo::PlaneGeo`.
     * @return an (undisclosed) object for iterating through all `geo::PlaneGeo`
     *
     * For example, this snippet computes `MaxWires()` of `TPC` (a `TPCGeo`):
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * unsigned int maxWires = 0U;
     * for (geo::PlaneGeo const& plane: TPC.IteratePlanes())
     *   maxWires = std::max(maxWires, plane.Nwires());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The resulting sequence exposes the planes within the TPC in their
     * ID order, from plane `0` to `Nplanes() - 1`.
     */
    ElementIteratorBox IterateElements() const;
    ElementIteratorBox IteratePlanes() const { return IterateElements(); }
    // @}

    /// Returns a set of all views covered in this TPC.
    std::set<geo::View_t> Views() const;

    /// @}

    /// @{
    /// @name TPC geometry properties

    /// Returns the center of the TPC volume in world coordinates [cm]
    Point_t GetCenter() const { return toWorldCoords(LocalPoint_t{}); }

    /// Returns the center of the TPC active volume in world coordinates [cm]
    Point_t GetActiveVolumeCenter() const { return fActiveCenter; }

    /// Returns the center of the active volume face opposite to the wire planes
    /// [cm]
    Point_t GetCathodeCenter() const { return GetCathodeCenterImpl(); }

    /// Returns the center of the active TPC volume side facing negative _z_.
    geo::Point_t GetFrontFaceCenter() const { return GetFrontFaceCenterImpl(); }

    /// Returns the bounding box of this TPC.
    geo::BoxBoundedGeo const& BoundingBox() const { return *this; }

    /// Returns the box of the active volume of this TPC.
    geo::BoxBoundedGeo const& ActiveBoundingBox() const { return fActiveBox; }

    double Plane0Pitch(unsigned int p) const;
    double PlanePitch(unsigned int p1 = 0, unsigned int p2 = 1) const;
    double WirePitch(unsigned plane = 0) const;

    /// Returns the identifier of this TPC
    geo::TPCID const& ID() const { return fID; }

    /// @}

    /// @{
    /// @name Projection on a wire plane
    ///
    /// These methods deal with projection of points and vectors on a plane,
    /// using a geometric reference base which is not dependent on the wire
    /// direction. Technically, the objects are projected on the reference
    /// plane, that happens to be the first wire plane.
    /// In practice, the important bit is that all entities are consistently
    /// decomposed into a drift component and a remaining component. If which
    /// plane the projection happens on should matter, then geo::PlaneGeo can be
    /// used directly (but note that geo::PlaneGeo defines _two_ different
    /// frames, and the names of the frame equivalent to the one used here are
    /// different).
    ///

    /// Returns the plane used for reference by projection methods.
    geo::PlaneGeo const& ReferencePlane() const { return FirstPlane(); }

    /// Returns the ID of the plane used for reference by projection methods.
    geo::PlaneID const& ReferencePlaneID() const { return ReferencePlane().ID(); }

    //@{
    /**
     * @brief Return the direction of reference plane width.
     * @tparam Vector type of vector to return (current default: `TVector3`)
     * @see RefDepthDir(), DriftDir()
     *
     * The precise definition of the vector is arbitrary, but it is defined
     * to lie on the wire plane and so that RefWidthDir(), RefDepthDir() and
     * a vector opposite to DriftDir() make a orthonormal base.
     * That base (width, depth, normal) is guaranteed to be positive defined.
     */
    Vector_t RefWidthDir() const { return ReferencePlane().WidthDir(); }
    //@}

    //@{
    /**
     * @brief Return the direction of reference plane depth.
     * @tparam Vector type of vector to return (current default: `TVector3`)
     * @see RefWidthDir(), DriftDir()
     *
     * The precise definition of the vector is arbitrary, but it is defined
     * to lie on the wire plane and so that RefWidthDir(), RefDepthDir() and
     * a vector opposite to DriftDir() make a orthonormal base.
     * That base (width, depth, normal) is guaranteed to be positive defined.
     */
    Vector_t RefDepthDir() const { return ReferencePlane().DepthDir(); }
    //@}

    //@{
    /**
     * @brief Returns the distance of the specified point from reference plane.
     * @param point a point in world coordinates [cm]
     * @return the signed distance from the plane
     *
     * The distance is defined positive if the point lies in the inner side
     * with respect to the reference plane.
     *
     * It should match the drift distance from that plane.
     */
    double DistanceFromReferencePlane(geo::Point_t const& point) const
    {
      return ReferencePlane().DistanceFromPlane(point);
    }
    //@}

    //@{
    /**
     * @brief Shifts the position of an electron drifted by a distance.
     * @param position _(modified)_ the position of the electron
     * @param distance drift distance to shift the electron by [cm]
     * @see geo::Plane::DriftPoint(), ReferencePlane()
     *
     * This operation is delegated to the reference plane
     * (see `geo::Plane::DriftPoint()`).
     */
    void DriftPoint(geo::Point_t& position, double distance) const
    {
      ReferencePlane().DriftPoint(position, distance);
    }
    //@}

    //@{
    /**
     * @brief Decomposes a 3D point in two components.
     * @param point the point to be decomposed
     * @return the two components of point, on the plane and orthogonal to it
     *
     * The point is decomposed in:
     *
     * 1. a component orthogonal to the reference plane, expressed as a signed
     *    real number
     * 2. a component lying on the reference plane, expressed as a 2D vector
     *
     * The distance is from the reference plane
     * (`DistanceFromReferencePlane()`).
     * The projection on the plane is obtained following the same convention
     * as `PointProjection()`.
     */
    DecomposedVector_t DecomposePoint(geo::Point_t const& point) const
    {
      return ReferencePlane().DecomposePointWidthDepth(point);
    }
    //@}

    //@{
    /**
     * @brief Returns the reference point used by `PointProjection()`.
     * @tparam Point type of point to be returned
     *
     * The returned point is such that its decomposition results in a null
     * projection and a 0 distance from the plane.
     */
    Point_t ProjectionReferencePoint() const { return ReferencePlane().GetCenter(); }
    //@}

    //@{
    /**
     * @brief Returns the projection of the specified point on the plane.
     * @param point the 3D point to be projected, in world coordinates
     * @return a 2D vector representing the projection of point on the plane
     *
     * The returned vector is a 2D vector expressing the projection of the point
     * (from world coordinates) on the reference plane.
     * The vector is expressed as @f$ ( w, d ) @f$, components following the
     * width direction (`RefWidthDir()`) and the depth direction
     * (`RefDepthDir()`) respectively. The origin point is returned by
     * `ProjectionReferencePoint()`.
     * All coordinates are in centimeters.
     */
    Projection_t Projection(geo::Point_t const& point) const
    {
      return ReferencePlane().PointWidthDepthProjection(point);
    }
    Projection_t PointProjection(geo::Point_t const& point) const { return Projection(point); }
    //@}

    //@{
    /**
     * @brief Returns the projection of the specified vector on the plane.
     * @param v the 3D vector to be projected, in world units [cm]
     * @return a 2D vector representing the projection of v on the plane
     *
     * The returned vector is a 2D vector expressing the projection of the
     * vector (from world units) on the reference plane.
     * The vector is expressed as @f$ ( w, d ) @f$, components following the
     * width direction (`RefWidthDir()`) and the depth direction
     * (`RefDepthDir()`) respectively.
     * All coordinates are in centimeters.
     */
    Projection_t Projection(geo::Vector_t const& v) const
    {
      return ReferencePlane().VectorWidthDepthProjection(v);
    }
    Projection_t VectorProjection(geo::Vector_t const& v) const { return Projection(v); }
    //@}

    //@{
    /**
     * @brief Returns the 3D vector from composition of projection and distance.
     * @tparam Point type of point to be returned
     * @param decomp decomposed point
     * @return the 3D vector from composition of projection and distance
     * @see DecomposePoint(), ComposePoint(double, Projection_t const&)
     *
     * This is the "inverse" operation respect to `DecomposePoint()`.
     * The argument stores the two components, orthogonal and parallel to the
     * plane, and compose them into a 3D point which departs from the reference
     * point (`ProjectionReferencePoint()`) by those components.
     * See `ComposePoint(double, Projection_t const&)` for more details.
     */
    Point_t ComposePoint(DecomposedVector_t const& decomp) const
    {
      return ReferencePlane().ComposePoint(decomp);
    }
    //@}

    //@{
    /**
     * @brief Returns the 3D point from composition of projection and distance.
     * @tparam Point type of point to be returned
     * @param distance distance of the target point from the reference plane
     * @param proj projection of the target point on the reference plane
     * @return the 3D point from composition of projection and distance
     * @see DecomposePoint()
     *
     * The returned point is the reference point, translated by two 3D vectors:
     *
     * 1. a vector parallel to the plane normal, with norm the input distance
     * 2. a vector lying on the plane, whose projection via `PointProjection()`
     *    gives the input projection
     *
     * The choice of the projection reference point embodies the same convention
     * used in `PointProjection()` and `DecomposePoint()`.
     *
     */
    Point_t ComposePoint(double distance, Projection_t const& proj) const
    {
      return ReferencePlane().ComposePoint(distance, proj);
    }
    //@}

    /// @}

    /// @{
    /// @name Coordinate transformation

    /// Transform point from local TPC frame to world frame.
    geo::Point_t toWorldCoords(LocalPoint_t const& local) const
    {
      return fTrans.toWorldCoords(local);
    }

    /// Transform direction vector from local to world.
    geo::Vector_t toWorldCoords(LocalVector_t const& local) const
    {
      return fTrans.toWorldCoords(local);
    }

    /// Transform point from world frame to local TPC frame.
    LocalPoint_t toLocalCoords(geo::Point_t const& world) const
    {
      return fTrans.toLocalCoords(world);
    }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(geo::Vector_t const& world) const
    {
      return fTrans.toLocalCoords(world);
    }

    /// @}

    /**
     * @brief Returns the expected drift direction based on geometry
     *
     * The return value is coded as follow:
     *
     * * +1: positive x
     * * +2: positive y
     * * +3: positive z
     * * -1: negative x
     * * -2: negative y
     * * -3: negative z
     * *  0: other (or algorithm failed)
     *
     * The current implementation is based on the assumption that electrons in
     * the middle of TPC will drift toward the wire planes, and it "never
     * fails".
     */
    short int DetectDriftDirection() const;

    /// Apply sorting to the PlaneGeo objects
    void SortSubVolumes(geo::GeoObjectSorter const& sorter);

    /// Performs all updates after cryostat has sorted TPCs
    void UpdateAfterSorting(geo::TPCID tpcid);

    /**
     * @brief Prints information about this TPC.
     * @tparam Stream type of output stream to use
     * @param out stream to send the information to
     * @param indent prepend each line with this string
     * @param verbosity amount of information printed
     *
     * Note that the first line out the output is _not_ indented.
     *
     * Verbosity levels
     * -----------------
     *
     * * 0: only TPC ID
     * * 1 _(default)_: also center and size
     * * 2: also drift direction, cathode position and number of planes
     * * 3: also maximum number of wires per plane
     * * 4: also information on main direction
     * * 5: also information on bounding box
     *
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintTPCInfo(Stream&& out, std::string indent = "", unsigned int verbosity = 1) const;

    /**
     * @brief Returns a string with information about this TPC.
     * @see `PrintTPCInfo()`
     *
     * Arguments and provided information are the same as in `PrintTPCInfo()`.
     */
    std::string TPCInfo(std::string indent = "", unsigned int verbosity = 1) const;

    /// Maximum verbosity supported by `PrintTPCInfo()`.
    static constexpr unsigned int MaxVerbosity = 6;

    /**
     * @brief Returns whether the specified coordinate is in a range
     * @param c the coordinate
     * @param min lower boundary of the range
     * @param max upper boundary of the range
     * @param wiggle expansion factor for the range
     * @return whether the specified coordinate is in a range
     *
     * If the wiggle is larger than 1, the range is expanded by the wiggle factor.
     * If the wiggle is less than 1, the range is shrinked.
     */
    static bool CoordinateContained(double c, double min, double max, double wiggle = 1.)
    {
      return (c >= (min > 0 ? min / wiggle : min * wiggle)) &&
             (c <= (max < 0 ? max / wiggle : max * wiggle));
    } // CoordinateContained()

    static bool CoordinateContained(double c, double const* range, double wiggle = 1.)
    {
      return CoordinateContained(c, range[0], range[1], wiggle);
    }

  private:
    void FindPlane(GeoNodePath_t& path, size_t depth);
    void MakePlane(GeoNodePath_t& path, size_t depth);

  private:
    using LocalTransformation_t =
      geo::LocalTransformationGeo<ROOT::Math::Transform3D, LocalPoint_t, LocalVector_t>;

    LocalTransformation_t fTrans; ///< TPC-to-world transformation.

    PlaneCollection_t fPlanes;        ///< List of planes in this plane.
    TGeoVolume* fActiveVolume;        ///< Active volume of LAr, called volTPCActive in GDML file.
    TGeoVolume* fTotalVolume;         ///< Total volume of TPC, called volTPC in GDML file.
    DriftDirection_t fDriftDirection; ///< Direction of the electron drift in the TPC.
    std::vector<double> fPlane0Pitch; ///< Pitch between planes.
    std::vector<std::vector<double>> fPlaneLocation; ///< xyz locations of planes in the TPC.
    geo::Point_t fActiveCenter; ///< Center of the active volume, in world coordinates [cm].

    double fActiveHalfWidth;  ///< Half width of active volume.
    double fActiveHalfHeight; ///< Half height of active volume.
    double fActiveLength;     ///< Length of active volume.
    double fHalfWidth;        ///< Half width of total volume.
    double fHalfHeight;       ///< Half height of total volume.
    double fLength;           ///< Length of total volume.

    geo::Vector_t fWidthDir;  ///< Direction width refers to.
    geo::Vector_t fHeightDir; ///< Direction height refers to.
    geo::Vector_t fLengthDir; ///< Direction length refers to.
    geo::Vector_t fDriftDir;  ///< Direction electrons drift along.

    geo::BoxBoundedGeo fActiveBox; ///< Box of the active volume.

    geo::TPCID fID; ///< ID of this TPC.

    /// Index of the plane for each view (InvalidID if none).
    std::vector<geo::PlaneID::PlaneID_t> fViewToPlaneNumber;

    /// Recomputes the drift direction; needs planes to have been initialised.
    void ResetDriftDirection();

    /// Computes the distance between the cathode and the last wire plane
    /// (last respect to the sorting order).
    double ComputeDriftDistance() const;

    /// Refills the plane vs. view cache of the TPC.
    void UpdatePlaneViewCache();

    /// Updates plane cached information.
    void UpdatePlaneCache();

    /// Recomputes the TPC boundary.
    void InitTPCBoundaries();

    /// Sorts (in place) the specified `PlaneGeo` objects by drift distance.
    void SortPlanes(std::vector<geo::PlaneGeo>&) const;

    geo::Point_t GetFrontFaceCenterImpl() const;
    geo::Point_t GetCathodeCenterImpl() const;
  };
}

//------------------------------------------------------------------------------
//--- template implementation
//---
//------------------------------------------------------------------------------
template <typename Stream>
void geo::TPCGeo::PrintTPCInfo(Stream&& out,
                               std::string indent /* = "" */,
                               unsigned int verbosity /* = 1 */
                               ) const
{

  //----------------------------------------------------------------------------
  out << "TPC " << std::string(ID());

  if (verbosity-- <= 0) return; // 0

  //----------------------------------------------------------------------------
  out << " (" << Width() << " x " << Height() << " x " << Length() << ") cm^3 at " << GetCenter();

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------

  out << "\n"
      << indent << "drift direction " << DriftDir() << " from cathode around " << GetCathodeCenter()
      << " through " << DriftDistance() << " cm toward " << Nplanes() << " wire planes";

  if (verbosity-- <= 0) return; // 2

  //----------------------------------------------------------------------------
  out << "\n" << indent << "maximum wires on any plane: " << MaxWires();

  if (verbosity-- <= 0) return; // 3

  //----------------------------------------------------------------------------
  out << "\n"
      << indent << "active volume (" << ActiveWidth() << " x " << ActiveHeight() << " x "
      << ActiveLength() << ") cm^3, front face at " << GetFrontFaceCenter() << " cm;"
      << "\n"
      << indent << "main directions:"
      << " width " << WidthDir() << " height " << HeightDir() << " length " << LengthDir();

  if (verbosity-- <= 0) return; // 4

  //----------------------------------------------------------------------------
  // print also the containing box
  geo::BoxBoundedGeo const& box = BoundingBox();
  out << "\n" << indent << "bounding box: " << box.Min() << " -- " << box.Max();

  //----------------------------------------------------------------------------
  // print also the active box
  geo::BoxBoundedGeo const& activeBox = ActiveBoundingBox();
  out << "\n" << indent << "active volume box: " << activeBox.Min() << " -- " << activeBox.Max();

  //----------------------------------------------------------------------------
} // geo::TPCGeo::PrintTPCInfo()

#endif // LARCOREALG_GEOMETRY_TPCGEO_H
////////////////////////////////////////////////////////////////////////
