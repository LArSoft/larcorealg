////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/AuxDetGeo.h
/// \brief Encapsulate the geometry of an auxiliary detector
/// \ingroup Geometry
///
/// \author  miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_AUXDETGEO_H
#define LARCOREALG_GEOMETRY_AUXDETGEO_H

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/GeoVectorLocalTransformation.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// ROOT libraries
#include "TGeoVolume.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"

// C/C++ libraries
#include <string>
#include <vector>
#include <type_traits>

class TGeoNode;

namespace geo {

  class GeoObjectSorter;

  /// \ingroup Geometry
  class AuxDetGeo {
  public:

    /// Type of list of sensitive volumes.
    using AuxDetSensitiveList_t = std::vector<geo::AuxDetSensitiveGeo>;

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
     * instances of `geo::AuxDetGeo` have the same type but are not compatible.
     */

    /// Tag for vectors in the "local" GDML coordinate frame of the
    /// auxiliary detector.
    struct AuxDetGeoCoordinatesTag {};

    /// Type of points in the local GDML auxiliary detector frame.
    using LocalPoint_t = geo::Point3DBase_t<AuxDetGeoCoordinatesTag>;

    /// Type of displacement vectors in the local GDML auxiliary detector frame.
    using LocalVector_t = geo::Vector3DBase_t<AuxDetGeoCoordinatesTag>;

    ///@}

    AuxDetGeo(
      TGeoNode const& node, geo::TransformationMatrix&& trans,
      AuxDetSensitiveList_t&& sensitive
      );


    /**
     * @brief Return the center position of an AuxDet.
     * @param xyz _(output)_ the returned location: `{ x, y, z }` [cm]
     * @param localz (default: `0`) distance along the length of the volume (z)
     *               [cm]
     */
    void GetCenter(double* xyz, double localz=0.0)    const;

    /**
     * @brief Returns the geometric center of the sensitive volume.
     * @param localz (default: `0`) distance from the center along the length
     *               of the volume (z) [cm]
     * @return the geometric center of the sensitive volume [cm]
     */
    geo::Point_t GetCenter(double localz = 0.0) const;

    /// Returns the unit normal vector to the detector.
    geo::Vector_t GetNormalVector() const;

    /// Fills the unit normal vector to the detector.
    void GetNormalVector(double* xyzDir) const;


    /// @{
    /// @name Box geometry
    double Length()                  const { return fLength;      }
    double HalfWidth1()      	       const { return fHalfWidth1;  }
    double HalfWidth2()      	       const { return fHalfWidth2;  }
    double HalfHeight()      	       const { return fHalfHeight;  }
    const TGeoVolume* TotalVolume()  const { return fTotalVolume; }
    /// @}

    //@{
    /// Returns the distance of `point` from the center of the detector.
    geo::Length_t DistanceToPoint(geo::Point_t const& point) const
      { return (point - GetCenter()).R(); }
    geo::Length_t DistanceToPoint(double const* point) const;
    //@}

    /// @{
    /// @name Coordinate transformation

    /// Transform point from local auxiliary detector frame to world frame.
    void LocalToWorld(const double* auxdet, double* world) const
      { fTrans.LocalToWorld(auxdet, world); }

    /// Transform point from local auxiliary detector frame to world frame.
    geo::Point_t toWorldCoords(LocalPoint_t const& local) const
      { return fTrans.toWorldCoords(local); }

    /// Transform direction vector from local to world.
    void LocalToWorldVect(const double* auxdet, double* world) const
      { fTrans.LocalToWorldVect(auxdet, world); }

    /// Transform direction vector from local to world.
    geo::Vector_t toWorldCoords(LocalVector_t const& local) const
      { return fTrans.toWorldCoords(local); }

    /// Transform point from world frame to local auxiliary detector frame.
    void WorldToLocal(const double* world, double* auxdet) const
      { fTrans.WorldToLocal(world, auxdet); }

    /// Transform point from world frame to local auxiliary detector frame.
    LocalPoint_t toLocalCoords(geo::Point_t const& world) const
      { return fTrans.toLocalCoords(world); }

    /// Transform direction vector from world to local.
    void WorldToLocalVect(const double* world, double* auxdet) const
      { fTrans.WorldToLocalVect(world, auxdet); }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(geo::Vector_t const& world) const
      { return fTrans.toLocalCoords(world); }

    /// @}


    std::string Name() const { return fTotalVolume->GetName(); }

    /// @{
    /// @name Access to the sensitive volumes in the detector

    //@{
    std::size_t FindSensitiveVolume(geo::Point_t const& point)  const;
    /// @deprecated Use the version with `geo::Point_t` argument instead
    std::size_t FindSensitiveVolume(double const worldLoc[3])  const;
    //@}
    //@{
    AuxDetSensitiveGeo const& PositionToSensitiveVolume
      (geo::Point_t const& point, size_t& sv) const;
    /// @deprecated Use the version with `geo::Point_t` argument instead
    AuxDetSensitiveGeo const& PositionToSensitiveVolume
      (double const worldLoc[3], size_t& sv) const;
    //@}
    AuxDetSensitiveGeo const& SensitiveVolume(size_t sv) const { return fSensitive[sv];   }
    size_t                    NSensitiveVolume()         const { return fSensitive.size(); }

    /// @}

    void                      SortSubVolumes(geo::GeoObjectSorter const& sorter);


    /**
     * @brief Prints information about this auxiliary detector.
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
     * * 0: only detector name
     * * 1 _(default)_: also center
     * * 2: also size
     * * 3: also number of sensitive detectors
     * * 4: also normal direction
     *
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintAuxDetInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1) const;

    /**
     * @brief Returns a string with auxiliary detector information.
     * @see `PrintAuxDetInfo()`
     *
     * The arguments and provided information are the same as in
     * `PrintAuxDetInfo()`.
     */
    std::string AuxDetInfo
      (std::string indent = "", unsigned int verbosity = 1) const;

    /// Maximum verbosity supported by `PrintAuxDetInfo()`.
    static constexpr unsigned int MaxVerbosity = 4;


  private:

    using LocalTransformation_t = geo::LocalTransformationGeo
      <geo::TransformationMatrix, LocalPoint_t, LocalVector_t>;

    const TGeoVolume*     	         fTotalVolume; ///< Total volume of AuxDet, called vol*
    LocalTransformation_t 	         fTrans;       ///< Auxiliary detector-to-world transformation.
    double                	         fLength;      ///< length of volume, along z direction in local
    double                	         fHalfWidth1;  ///< 1st half width of volume, at -z/2 in local coordinates
    double                	         fHalfWidth2;  ///< 2nd half width (width1==width2 for boxes), at +z/2
    double                	         fHalfHeight;  ///< half height of volume
    std::vector<AuxDetSensitiveGeo> fSensitive;   ///< sensitive volumes in the detector

    /// Extracts the size of the detector from the geometry information.
    void InitShapeSize();
  }; // class AuxDetGeo

  static_assert(std::is_move_assignable_v<geo::AuxDetGeo>);
  static_assert(std::is_move_constructible_v<geo::AuxDetGeo>);

}


//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::AuxDetGeo::PrintAuxDetInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1 */
) const {

  //----------------------------------------------------------------------------
  out << "\"" << Name() << "\"";

  if (verbosity-- <= 0) return; // 0

  //----------------------------------------------------------------------------
  out << " centered at " << GetCenter() << " cm";

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------
  lar::util::RealComparisons<double> coordIs(1e-4);
  out << ", size ( " << (2.0 * HalfWidth1());
  if (coordIs.nonEqual(HalfWidth1(), HalfWidth2()))
    out << "/" << (2.0 * HalfWidth2());
  out << " x " << (2.0 * HalfHeight()) << " x " << Length() << " ) cm";

  if (verbosity-- <= 0) return; // 2

  //----------------------------------------------------------------------------
  out << "\n" << indent
    << "with ";
  switch (NSensitiveVolume()) {
    case 0:  out << "no sensitive volume"; break;
    case 1:  out << "1 sensitive volume"; break;
    default: out << NSensitiveVolume() << " sensitive volumes"; break;
  } // switch

  if (verbosity-- <= 0) return; // 3

  //----------------------------------------------------------------------------
  out << ", normal facing " << GetNormalVector();

//  if (verbosity-- <= 0) return; // 4

  //----------------------------------------------------------------------------

} // geo::AuxDetGeo::PrintAuxDetInfo()


//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_AUXDETGEO_H
