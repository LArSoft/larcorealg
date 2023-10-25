////////////////////////////////////////////////////////////////////////
/// @file  larcorealg/Geometry/AuxDetSensitiveGeo.h
/// @brief Encapsulate the geometry of the sensitive portion of an auxiliary detector
/// @ingroup Geometry
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_AUXDETSENSITIVEGEO_H
#define LARCOREALG_GEOMETRY_AUXDETSENSITIVEGEO_H

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/Geometry/GeoVectorLocalTransformation.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// ROOT libraries
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/Transform3D.h"

// C/C++ standard libraries
#include <string>
#include <type_traits>

class TGeoNode;
class TGeoVolume;

namespace geo {

  /// @ingroup Geometry
  class AuxDetSensitiveGeo {
  public:
    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     *
     * These types represents points and displacement vectors in the reference
     * frame defined in the auxiliary detector geometry box from the GDML
     * geometry description.
     *
     * No alias is explicitly defined for the LArSoft global vector types,
     * `geo::Point_t` and `geo::Vector_t`.
     *
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different
     * instances of `geo::AuxDetSensitiveGeo` have the same type but are not
     * compatible.
     */

    /// Tag for vectors in the "local" GDML coordinate frame of the
    /// auxiliary detector.
    struct AuxDetSensitiveGeoCoordinatesTag {};

    /// Type of points in the local GDML auxiliary detector frame.
    using LocalPoint_t = Point3DBase_t<AuxDetSensitiveGeoCoordinatesTag>;

    /// Type of displacement vectors in the local GDML auxiliary detector frame.
    using LocalVector_t = Vector3DBase_t<AuxDetSensitiveGeoCoordinatesTag>;

    ///@}

    AuxDetSensitiveGeo(TGeoNode const* node, TransformationMatrix&& trans);

    /**
     * @brief Returns the geometric center of the sensitive volume.
     * @param localz (default: `0`) distance from the center along the length
     *               of the volume (z) [cm]
     * @return the geometric center of the sensitive volume [cm]
     */
    Point_t GetCenter(double localz = 0.0) const;

    /// Returns the unit normal vector to the detector.
    Vector_t GetNormalVector() const;

    //box geometry
    double Length() const { return fLength; }
    double HalfLength() const { return Length() / 2.0; }
    double HalfWidth1() const { return fHalfWidth1; }
    double HalfWidth2() const { return fHalfWidth2; }
    double HalfCenterWidth() const { return (HalfWidth1() + HalfWidth2()) / 2.0; }
    double HalfHeight() const { return fHalfHeight; }
    const TGeoVolume* TotalVolume() const { return fTotalVolume; }

    //@{
    /// Returns the distance of `point` from the center of the detector.
    Length_t DistanceToPoint(Point_t const& point) const { return (point - GetCenter()).R(); }
    Length_t DistanceToPoint(double const* point) const;
    //@}

    /// @{
    /// @name Coordinate transformation

    /// Transform point from local auxiliary detector frame to world frame.
    Point_t toWorldCoords(LocalPoint_t const& local) const { return fTrans.toWorldCoords(local); }

    /// Transform direction vector from local to world.
    Vector_t toWorldCoords(LocalVector_t const& local) const { return fTrans.toWorldCoords(local); }

    /// Transform point from world frame to local auxiliary detector frame.
    LocalPoint_t toLocalCoords(Point_t const& world) const { return fTrans.toLocalCoords(world); }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(Vector_t const& world) const { return fTrans.toLocalCoords(world); }

    /// @}

    /**
     * @brief Prints information about this auxiliary sensitive detector.
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
     * * 0 _(default)_: only center
     * * 1: also size
     * * 2: also normal direction
     *
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintAuxDetInfo(Stream&& out, std::string indent = "", unsigned int verbosity = 0) const;

    /**
     * @brief Returns a string with auxiliary sensitive detector information.
     * @see `PrintAuxDetInfo()`
     *
     * The arguments and provided information are the same as in
     * `PrintAuxDetInfo()`.
     */
    std::string AuxDetInfo(std::string indent = "", unsigned int verbosity = 1) const;

    /// Maximum verbosity supported by `PrintAuxDetInfo()`.
    static constexpr unsigned int MaxVerbosity = 2;

  private:
    using LocalTransformation_t =
      LocalTransformationGeo<ROOT::Math::Transform3D, LocalPoint_t, LocalVector_t>;

    LocalTransformation_t fTrans;   ///< Auxiliary detector-to-world transformation.
    const TGeoVolume* fTotalVolume; ///< Total volume of AuxDet, called vol*
    double fLength;                 ///< length of volume, along z direction in local
    double fHalfWidth1;             ///< 1st half width of volume, at -z/2 in local coordinates
    double fHalfWidth2;             ///< 2nd half width (width1==width2 for boxes), at +z/2
    double fHalfHeight;             ///< half height of volume

    /// Extracts the size of the detector from the geometry information.
    void InitShapeSize();

  }; // class AuxDetSensitiveGeo

  static_assert(std::is_move_assignable_v<AuxDetSensitiveGeo>);
  static_assert(std::is_move_constructible_v<AuxDetSensitiveGeo>);

} // namespace geo

//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::AuxDetSensitiveGeo::PrintAuxDetInfo(Stream&& out,
                                              std::string indent /* = "" */,
                                              unsigned int verbosity /* = 0 */
                                              ) const
{

  //----------------------------------------------------------------------------
  out << "centered at " << GetCenter() << " cm";

  if (verbosity-- <= 0) return; // 0

  //----------------------------------------------------------------------------
  lar::util::RealComparisons<double> coordIs(1e-4);
  out << ", size ( " << (2.0 * HalfWidth1());
  if (coordIs.nonEqual(HalfWidth1(), HalfWidth2())) out << "/" << (2.0 * HalfWidth2());
  out << " x " << (2.0 * HalfHeight()) << " x " << Length() << " ) cm";

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------
  out << ", normal facing " << GetNormalVector();

  //  if (verbosity-- <= 0) return; // 2

  //----------------------------------------------------------------------------

} // geo::AuxDetSensitiveGeo::PrintAuxDetInfo()

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_AUXDETSENSITIVEGEO_H
