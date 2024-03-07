////////////////////////////////////////////////////////////////////////
/// @file  larcorealg/Geometry/TPCGeo.h
/// @brief Encapsulate the construction of a single detector plane
/// @ingroup Geometry
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_TPCGEO_H
#define LARCOREALG_GEOMETRY_TPCGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
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

    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     *
     * These types represents points and displacement vectors in the reference frame
     * defined in the TPC geometry box from the GDML geometry description.
     *
     * No alias is explicitly defined for the LArSoft global vector types, `geo::Point_t`
     * and `geo::Vector_t`.
     *
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different instances of
     * `geo::TPCGeo` have the same type but are not compatible.
     */

    /// Tag for vectors in the "local" GDML coordinate frame of the TPC.
    struct TPCGeoCoordinatesTag {};

    /// Type of points in the local GDML TPC frame.
    using LocalPoint_t = Point3DBase_t<TPCGeoCoordinatesTag>;

    /// Type of displacement vectors in the local GDML TPC frame.
    using LocalVector_t = Vector3DBase_t<TPCGeoCoordinatesTag>;

    ///@}

    // Construct a representation of a single plane of the detector
    TPCGeo(TGeoNode const* tpc_node,
           std::size_t hash_value,
           TransformationMatrix&& trans,
           DriftAxis driftAxis,
           double driftDistance);

    /// @{
    /// @name TPC properties

    /**
     * @brief Returns the expected drift direction based on geometry
     *
     * The current implementation is based on the assumption that electrons in the middle
     * of TPC will drift toward the wire planes, and it "never fails".
     */
    DriftAxis DriftAxisWithSign() const { return fDriftAxis; }
    geo::DriftSign DriftSign() const { return DriftAxisWithSign().sign; }

    Point_t GetCathodeCenter() const;

    /// Returns the direction of the drift (vector pointing toward the planes).
    Vector_t DriftDir() const { return fDriftDir; }

    /// Drift distance is defined as the distance between the anode and the cathode, in centimeters.
    double DriftDistance() const;

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

    /// @}

    /// @{
    /// @name TPC geometry properties

    /// Returns the center of the TPC volume in world coordinates [cm]
    Point_t GetCenter() const { return toWorldCoords(LocalPoint_t{}); }

    /// Returns the center of the TPC active volume in world coordinates [cm]
    Point_t GetActiveVolumeCenter() const { return fActiveCenter; }

    /// Returns the center of the active TPC volume side facing negative _z_.
    Point_t GetFrontFaceCenter() const;

    /// Returns the bounding box of this TPC.
    BoxBoundedGeo const& BoundingBox() const { return *this; }

    /// Returns the box of the active volume of this TPC.
    BoxBoundedGeo const& ActiveBoundingBox() const { return fActiveBox; }

    /// Returns the identifier of this TPC
    TPCID const& ID() const { return fID; }

    /// @}

    /// @{
    /// @name Coordinate transformation

    /// Transform point from local TPC frame to world frame.
    Point_t toWorldCoords(LocalPoint_t const& local) const { return fTrans.toWorldCoords(local); }

    /// Transform direction vector from local to world.
    Vector_t toWorldCoords(LocalVector_t const& local) const { return fTrans.toWorldCoords(local); }

    /// Transform point from world frame to local TPC frame.
    LocalPoint_t toLocalCoords(Point_t const& world) const { return fTrans.toLocalCoords(world); }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(Vector_t const& world) const { return fTrans.toLocalCoords(world); }

    /// @}

    /// Performs all updates after cryostat has sorted TPCs
    void UpdateAfterSorting(TPCID tpcid);

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
     * The constant `MaxVerbosity` is set to the highest supported verbosity level.
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
     * If the wiggle is larger than 1, the range is expanded by the wiggle factor.  If the
     * wiggle is less than 1, the range is shrinked.
     */
    static bool CoordinateContained(double c, double min, double max, double wiggle = 1.)
    {
      return (c >= (min > 0 ? min / wiggle : min * wiggle)) &&
             (c <= (max < 0 ? max / wiggle : max * wiggle));
    }

    static bool CoordinateContained(double c, double const* range, double wiggle = 1.)
    {
      return CoordinateContained(c, range[0], range[1], wiggle);
    }

    std::size_t Hash() const { return fHash; }
    static TGeoNode const* NodeForActiveVolume(TGeoNode const* tpc);

  private:
    using LocalTransformation_t =
      LocalTransformationGeo<ROOT::Math::Transform3D, LocalPoint_t, LocalVector_t>;

    std::size_t fHash;            ///< Uniquely identifies TPC before sorting has been performed.
    LocalTransformation_t fTrans; ///< TPC-to-world transformation.

    DriftAxis fDriftAxis;
    Vector_t fDriftDir;
    double fDriftDistance;

    TGeoVolume* fActiveVolume{nullptr}; ///< Active volume of LAr, called volTPCActive in GDML file.
    TGeoVolume* fTotalVolume{nullptr};  ///< Total volume of TPC, called volTPC in GDML file.
    Point_t fActiveCenter;              ///< Center of the active volume, in world coordinates [cm].

    double fActiveHalfWidth;  ///< Half width of active volume.
    double fActiveHalfHeight; ///< Half height of active volume.
    double fActiveLength;     ///< Length of active volume.
    double fHalfWidth;        ///< Half width of total volume.
    double fHalfHeight;       ///< Half height of total volume.
    double fLength;           ///< Length of total volume.

    Vector_t fWidthDir{Xaxis()};  ///< Direction width refers to.
    Vector_t fHeightDir{Yaxis()}; ///< Direction height refers to.
    Vector_t fLengthDir{Zaxis()}; ///< Direction length refers to.

    BoxBoundedGeo fActiveBox; ///< Box of the active volume.

    TPCID fID; ///< ID of this TPC.

    /// Recomputes the TPC boundary.
    void InitTPCBoundaries();
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

  if (verbosity-- <= 0) return; // 2

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
  BoxBoundedGeo const& box = BoundingBox();
  out << "\n" << indent << "bounding box: " << box.Min() << " -- " << box.Max();

  //----------------------------------------------------------------------------
  // print also the active box
  BoxBoundedGeo const& activeBox = ActiveBoundingBox();
  out << "\n" << indent << "active volume box: " << activeBox.Min() << " -- " << activeBox.Max();

  //----------------------------------------------------------------------------
} // geo::TPCGeo::PrintTPCInfo()

#endif // LARCOREALG_GEOMETRY_TPCGEO_H
////////////////////////////////////////////////////////////////////////
