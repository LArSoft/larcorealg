////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/CryostatGeo.h
/// \brief Encapsulate the construction of a single cyostat
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARCOREALG_GEOMETRY_CRYOSTATGEO_H
#define LARCOREALG_GEOMETRY_CRYOSTATGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeoVectorLocalTransformation.h" // for LocalT...
#include "larcorealg/Geometry/LocalTransformationGeo.h"       // for LocalT...
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/WireGeo.h"           // for WireGeo
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t

// ROOT libraries
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/Transform3D.h"
#include "TGeoVolume.h"

// C/C++ standard libraries
#include <string>
#include <vector>

// forward declarations
class TGeoNode;

namespace geo {

  class GeoObjectSorter;

  //......................................................................
  /// @brief Geometry information for a single cryostat.
  /// @ingroup Geometry
  class CryostatGeo : public BoxBoundedGeo {

  public:
    using ID_t = CryostatID;

    /// Type used internally to store the TPCs.
    using TPCList_t = std::vector<TPCGeo>;

    /// Type used internally to store the optical detectors.
    using OpDetList_t = std::vector<OpDetGeo>;

    using GeoNodePath_t = WireGeo::GeoNodePath_t;

    /// Type returned by `IterateElements()`.
    using ElementIteratorBox = TPCList_t const&;

    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     *
     * These types represents points and displacement vectors in the reference
     * frame defined in the cryostat geometry box from the GDML geometry
     * description.
     *
     * No alias is explicitly defined for the LArSoft global vector types,
     * `geo::Point_t` and `geo::Vector_t`.
     *
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different
     * instances of `geo::CryostatGeo` have the same type but are not
     * compatible.
     */

    /// Tag for vectors in the "local" GDML coordinate frame of the cryostat.
    struct CryostatGeoCoordinatesTag {};

    /// Type of points in the local GDML cryostat frame.
    using LocalPoint_t = Point3DBase_t<CryostatGeoCoordinatesTag>;

    /// Type of displacement vectors in the local GDML cryostat frame.
    using LocalVector_t = Vector3DBase_t<CryostatGeoCoordinatesTag>;

    ///@}

    /// Construct a representation of a single cryostat of the detector.
    CryostatGeo(TGeoNode const* node,
                TransformationMatrix&& trans,
                TPCList_t&& TPCs,
                OpDetList_t&& OpDets);

    /// @{
    /// @name Cryostat geometry information

    /// Half width of the cryostat [cm]
    double HalfWidth() const;
    /// Half height of the cryostat [cm]
    double HalfHeight() const;
    /// Half height of the cryostat [cm]
    double HalfLength() const;
    /// Full width of the cryostat [cm]
    double Width() const { return 2. * HalfWidth(); }
    /// Full height of the cryostat [cm]
    double Height() const { return 2. * HalfHeight(); }
    /// Length of the cryostat [cm]
    double Length() const { return 2. * HalfLength(); }
    /// Mass of the cryostat
    double Mass() const { return fVolume->Weight(); }
    /// Pointer to ROOT's volume descriptor
    const TGeoVolume* Volume() const { return fVolume; }

    /// @brief Returns boundaries of the cryostat (in centimetres).
    /// @return boundaries in a geo::BoxBoundedGeo
    BoxBoundedGeo const& Boundaries() const { return BoundingBox(); }

    /// @brief Fills boundaries of the cryostat (in centimetres).
    /// @param boundaries filled as: [0] -x [1] +x [2] -y [3] +y [4] -z [5] +z
    void Boundaries(double* boundaries) const;

    /// Returns the geometrical center of the cryostat.
    Point_t GetCenter() const { return Boundaries().Center(); }

    /// Returns the bounding box of this cryostat.
    BoxBoundedGeo const& BoundingBox() const { return static_cast<BoxBoundedGeo const&>(*this); }

    /// Returns the identifier of this cryostat.
    CryostatID const& ID() const { return fID; }

    /**
     * @brief Prints information about this cryostat.
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
     * * 0: only cryostat ID
     * * 1 _(default)_: also center and size
     * * 2: also number of TPCs, optical detectors, and maximum wires per plane
     * *    and of planes for TPC
     * * 3: also information on bounding box
     *
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintCryostatInfo(Stream&& out, std::string indent = "", unsigned int verbosity = 1) const;

    /**
     * @brief Returns a string with cryostat information.
     * @see `PrintCryostatInfo()`
     *
     * The arguments and provided information are the same as in
     * `PrintCryostatInfo()`.
     */
    std::string CryostatInfo(std::string indent = "", unsigned int verbosity = 1) const;

    /// Maximum verbosity supported by `PrintCryostatInfo()`.
    static constexpr unsigned int MaxVerbosity = 3;

    /// @}

    // BEGIN TPC access --------------------------------------------------------
    /// @{
    /// @name TPC access

    /// Number of TPCs in this cryostat.
    unsigned int NTPC() const { return fTPCs.size(); }
    /// Alias for `NTPC()`.
    unsigned int NElements() const { return fTPCs.size(); }

    /**
     * @brief Returns whether a TPC with index itpc is present in this cryostat.
     * @param itpc index of TPC in this cryostat
     * @return whether the TPC with index itpc is present in this cryostat
     */
    bool HasTPC(unsigned int itpc) const { return itpc < NTPC(); }

    /// Alias for `HasTPC()`.
    bool HasElement(unsigned int itpc) const { return HasTPC(itpc); }

    /**
     * @brief Returns whether the TPC in tpcid is present in this cryostat
     * @param tpcid full TPC ID
     * @return whether the TPC in tpcid is present in this cryostat
     *
     * The cryostat number in tpcid is ignored, as it is ignored whether tpcid
     * is invalid.
     */
    bool HasTPC(TPCID const& tpcid) const { return HasTPC(tpcid.TPC); }
    /// Alias for `HasTPC(geo::TPCID const&)`
    bool HasElement(TPCID const& tpcid) const { return HasTPC(tpcid); }

    /// @brief Return the itpc'th TPC in the cryostat.
    /// @throws cet::exception (category "TPCOutOfRange") if no such TPC
    TPCGeo const& TPC(unsigned int itpc) const;

    /**
     * @brief Returns the TPC in tpcid from this cryostat
     * @param tpcid full TPC ID
     * @return a constant reference to the TPC in tpcid
     * @throws cet::exception (category "TPCOutOfRange") if no such TPC
     *
     * The cryostat number in tpcid is ignored, as it is ignored whether tpcid
     * is invalid.
     */
    TPCGeo const& TPC(TPCID const& tpcid) const { return TPC(tpcid.TPC); }
    /// Alias for `TPC()`.
    TPCGeo const& GetElement(TPCID const& tpcid) const { return TPC(tpcid); }

    /**
     * @brief Returns an object suitable for iterating through all TPCs.
     * @see `IterateTPCs()`
     *
     * The returned value can be used in a range-for loop like:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (geo::TPCGeo const& tpc: cryo.IterateElements()) { ... }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The resulting sequence exposes the TPCs within the cryostat in their
     * ID order, from TPC `0` to `NTPC() - 1`.
     *
     * This method is designed for templated code, where the object
     * `obj.IterateElements()` may be a `geo::CryostatGeo` or some other one.
     * For non-template code, prefer `IterateTPCs()` for clarity.
     */
    ElementIteratorBox IterateElements() const;

    /**
     * @brief Returns an object suitable for iterating through all TPCs.
     * @see `IterateElements()`
     *
     * The returned value can be used in a range-for loop like:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (geo::TPCGeo const& tpc: cryo.IterateTPCs()) { ... }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The resulting sequence exposes the TPCs within the cryostat in their
     * ID order, from TPC `0` to `NTPC() - 1`.
     *
     * A version of this functionality designed for template code is provided
     * under the generic name `IterateElements()`.
     */
    ElementIteratorBox IterateTPCs() const { return IterateElements(); }

    /**
     * @brief Returns the TPC number itpc from this cryostat.
     * @param itpc the number of local TPC
     * @return a constant pointer to the TPC, or nullptr if it does not exist
     */
    TPCGeo const* TPCPtr(unsigned int itpc) const
    {
      return HasTPC(itpc) ? &(fTPCs[itpc]) : nullptr;
    }

    /**
     * @brief Returns the TPC in tpcid from this cryostat.
     * @param tpcid full TPC ID
     * @return a constant pointer to the TPC, or nullptr if it does not exist
     *
     * The cryostat number in tpcid is ignored, as it is ignored whether tpcid
     * is invalid.
     */
    TPCGeo const* TPCPtr(TPCID const& tpcid) const { return TPCPtr(tpcid.TPC); }
    /// Alias for `TPCPtr()`.
    TPCGeo const* GetElementPtr(TPCID const& tpcid) const { return TPCPtr(tpcid); }

    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param point 3D coordinates of the point (world reference frame)
     * @param wiggle a small factor (like 1+epsilon) to avoid rounding errors
     * @return the ID of the TPC at the specified point (invalid ID if none)
     */
    TPCID PositionToTPCID(Point_t const& point, double wiggle) const;

    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param point the location (world reference frame)
     * @param wiggle a small factor (like 1+epsilon) to avoid rounding errors
     * @return the ID of the TPC at the specified point (invalid ID if none)
     */
    TPCGeo const& PositionToTPC(Point_t const& point, double wiggle) const;
    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param wiggle a small factor (like 1+epsilon) to avoid rounding errors
     * @return the ID of the TPC at the specified point (invalid ID if none)
     */
    TPCGeo const& PositionToTPC(double const worldLoc[3], double wiggle) const
    {
      return PositionToTPC(vect::makePointFromCoords(worldLoc), wiggle);
    }

    /**
     * @brief Returns a pointer to the TPC at specified location.
     * @param point position in space [cm]
     * @param wiggle a small factor (like 1+&epsilon;) to avoid rounding errors
     * @return a pointer to the `geo::TPCGeo` at `point` (`nullptr` if none)
     */
    TPCGeo const* PositionToTPCptr(Point_t const& point, double wiggle) const;

    /// Returns the largest number of planes among the TPCs in this cryostat
    unsigned int MaxPlanes() const;

    /// Returns the largest number of wires among the TPCs in this cryostat
    unsigned int MaxWires() const;

    /// @}
    // END TPC access ----------------------------------------------------------

    // BEGIN Optical detector access -------------------------------------------
    /// @{
    /// @name Optical detector access

    /// Number of optical detectors in this TPC
    unsigned int NOpDet() const { return fOpDets.size(); }

    /// Return the iopdet'th optical detector in the cryostat
    OpDetGeo const& OpDet(unsigned int iopdet) const;

    /// Returns the index of the optical detector in this cryostat closest to
    /// `point`.
    unsigned int GetClosestOpDet(Point_t const& point) const;
    /// @see `GetClosestOpDet(geo::Point_t const&) const`
    unsigned int GetClosestOpDet(double const* point) const;

    /// Returns the optical detector det in this cryostat nearest to `point`.
    /// If there are no optical detectors, `nullptr` is returned.
    OpDetGeo const* GetClosestOpDetPtr(Point_t const& point) const;

    /// Get name of opdet geometry element
    std::string OpDetGeoName() const { return fOpDetGeoName; }

    /// @}
    // END Optical detector access ---------------------------------------------

    // BEGIN Coordinate transformation -----------------------------------------
    /// @{
    /// @name Coordinate transformation

    /// Transform point from local cryostat frame to world frame.
    Point_t toWorldCoords(LocalPoint_t const& local) const { return fTrans.toWorldCoords(local); }

    /// Transform direction vector from local to world.
    Vector_t toWorldCoords(LocalVector_t const& local) const { return fTrans.toWorldCoords(local); }

    /// Transform point from world frame to local cryostat frame.
    LocalPoint_t toLocalCoords(Point_t const& world) const { return fTrans.toLocalCoords(world); }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(Vector_t const& world) const { return fTrans.toLocalCoords(world); }

    /// @}
    // END Coordinate transformation -------------------------------------------

    /// Method to sort TPCGeo objects
    void SortSubVolumes(Compare<TPCGeo> tpc_cmp, Compare<OpDetGeo> opdet_cmp);

    /// Performs all needed updates after geometry has sorted the cryostats
    void UpdateAfterSorting(CryostatID cryoid);

  private:
    void FindTPC(std::vector<const TGeoNode*>& path, unsigned int depth);
    void MakeTPC(std::vector<const TGeoNode*>& path, int depth);

    void FindOpDet(std::vector<const TGeoNode*>& path, unsigned int depth);
    void MakeOpDet(std::vector<const TGeoNode*>& path, int depth);

    /// Fill the boundary information of the cryostat
    void InitCryoBoundaries();

  private:
    using LocalTransformation_t =
      LocalTransformationGeo<ROOT::Math::Transform3D, LocalPoint_t, LocalVector_t>;

    LocalTransformation_t fTrans; ///< Cryostat-to-world transformation.
    TPCList_t fTPCs;              ///< List of tpcs in this cryostat
    OpDetList_t fOpDets;          ///< List of opdets in this cryostat
    TGeoVolume* fVolume;          ///< Total volume of cryostat, called volCryostat in GDML file
    std::string fOpDetGeoName{"volOpDetSensitive"}; ///< Name of opdet geometry elements in gdml
    CryostatID fID{};                               ///< ID of this cryostat
  };
}

//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::CryostatGeo::PrintCryostatInfo(Stream&& out,
                                         std::string indent /* = "" */,
                                         unsigned int verbosity /* = 1 */
                                         ) const
{
  //----------------------------------------------------------------------------
  out << "Cryostat " << std::string(ID());

  if (verbosity-- <= 0) return; // 0

  //----------------------------------------------------------------------------
  out << " (" << Width() << " x " << Height() << " x " << Length() << ") cm^3 at " << GetCenter();

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------

  out << "\n" << indent << "hosts " << NTPC() << " TPCs and " << NOpDet() << " optical detectors";

  if (verbosity-- <= 0) return; // 2

  //----------------------------------------------------------------------------
  // print also the containing box
  BoxBoundedGeo const& box = BoundingBox();
  out << "\n" << indent << "bounding box: " << box.Min() << " -- " << box.Max();
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_CRYOSTATGEO_H
////////////////////////////////////////////////////////////////////////
