////////////////////////////////////////////////////////////////////////
/// @file  larcorealg/Geometry/OpDetGeo.h
/// @brief Encapsulate the geometry of an optical detector
/// @ingroup Geometry
///
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_OPDETGEO_H
#define LARCOREALG_GEOMETRY_OPDETGEO_H

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_optical_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::OpDetID
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// ROOT libraries
#include "TClass.h"
#include "TGeoBBox.h"
#include "TGeoMatrix.h" // TGeoHMatrix
#include "TGeoSphere.h"
#include "TGeoTube.h"

// C/C++ standard libraries
#include <algorithm> // std::minmax()
#include <array>
#include <cassert>
#include <string>
#include <type_traits> // std::decay_t(), std::is_base_of_v
#include <typeinfo>    // typeid()
#include <vector>

// forward declarations
class TGeoNode;

namespace geo {

  /// @ingroup Geometry
  class OpDetGeo {
  public:
    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     *
     * These types represents points and displacement vectors in the reference
     * frame defined in the optical detector geometry box from the GDML geometry
     * description.
     *
     * No alias is explicitly defined for the LArSoft global vector types,
     * `geo::Point_t` and `geo::Vector_t`.
     *
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different
     * instances of `geo::OpDetGeo` have the same type but are not compatible.
     */

    /// Type of points in the local GDML TPC frame.
    using LocalPoint_t = OpticalPoint_t;

    /// Type of displacement vectors in the local GDML TPC frame.
    using LocalVector_t = OpticalVector_t;

    ///@}

    OpDetGeo(TGeoNode const* node, TransformationMatrix&& trans);

    /// Returns the geometry ID of this optical detector.
    OpDetID const& ID() const { return fID; }

    Point_t const& GetCenter() const { return fCenter; }
    double RMin() const;
    double RMax() const;
    double HalfL() const;
    double HalfW() const;
    double HalfH() const;
    double Length() const { return 2.0 * HalfL(); }
    double Width() const { return 2.0 * HalfW(); }
    double Height() const { return 2.0 * HalfH(); }
    double ThetaZ() const;             ///< returns angle of detector
                                       ///< with respect to z axis
                                       ///< in the Y-Z plane, in radians
    double ThetaZ(bool degrees) const; ///< returns angle of detector
                                       ///< with respect to z axis
                                       ///< in the Y-Z plane
    //@{
    /// Get cos(angle) to normal of this detector - used for solid angle calcs
    double CosThetaFromNormal(Point_t const& point) const;
    //@}
    //@{
    /// Returns the distance of the specified point from detector center [cm]
    double DistanceToPoint(Point_t const& point) const;
    //@}

    /// @{
    /**
     * @name Coordinate transformation
     *
     * Local points and displacement vectors are described by the types
     * `geo::OpDetGeo::LocalPoint_t` and `geo::OpDetGeo::LocalVector_t`,
     * respectively.
     */

    /// Transform point from local optical detector frame to world frame.
    Point_t toWorldCoords(LocalPoint_t const& local) const { return fTrans.toWorldCoords(local); }

    /// Transform direction vector from local to world.
    Vector_t toWorldCoords(LocalVector_t const& local) const { return fTrans.toWorldCoords(local); }

    /// Transform point from world frame to local optical detector frame.
    LocalPoint_t toLocalCoords(Point_t const& world) const { return fTrans.toLocalCoords(world); }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(Vector_t const& world) const { return fTrans.toLocalCoords(world); }

    /// @}

    /// Returns the ROOT object describing the detector geometry.
    const TGeoNode* Node() const { return fOpDetNode; }

    // --- BEGIN -- detector shape ---------------------------------------------
    /// @name Detector shape
    /// @{

    /// Returns the geometry object as `TGeoShape`.
    TGeoShape const* Shape() const { return Node()->GetVolume()->GetShape(); }

    /**
     * @brief Returns whether the detector has the specified shape.
     * @tparam ShapeObj type of ROOT geometry object representing the shape
     * @return whether this detector has the specified shape
     * @see `isShapeLike()`, `isBox()`, `isSphere()`, `isTube()`
     *
     * Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * bool const isSphere = opDet.isShape<TGeoSphere>();
     * bool const isBox = opDet.isShape<TGeoBBox>();
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will have `isSphere` `true` only if the shape of this object is a sphere
     * (`TGeoSphere`), and `isBox` `true` only if the shape of this object is a
     * box (`TGeoBBox`).
     */
    template <typename ShapeObj>
    bool isShape() const;

    /**
     * @brief Returns whether the detector inherits from the specified shape.
     * @tparam ShapeObj type of ROOT geometry object representing the shape
     * @return whether this detector has a shape derived from the specified one
     * @see `isShape()`, `isBox()`, `isSphere()`, `isTube()`
     *
     * Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * bool const isTubeLike = opDet.isShapeLike<TGeoTube>();
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * `isTubeLike` will be `true` if its shape is either a box (`TGeoTube`)
     * or any other whose shape object is derived from `TGeoTube` (including
     * for example a C-shape, half-cylinder).
     */
    template <typename ShapeObj>
    bool isShapeLike() const;

    /// Returns whether the detector shape is a cylinder (`TGeoTube`).
    bool isTube() const { return isShapeLike<TGeoTube>(); }

    /// Returns whether the detector shape is a bar (`TGeoBBox`).
    bool isBar() const { return isShape<TGeoBBox>(); }

    /// Returns whether the detector shape is a hemisphere (`TGeoSphere`).
    bool isSphere() const { return isShape<TGeoSphere>(); }

    /// @}
    // --- END -- detector shape -----------------------------------------------

    /// Performs all updates after cryostat has sorted the optical detectors.
    void UpdateAfterSorting(OpDetID opdetid);

    /**
     * @brief Prints information about this optical detector.
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
     * * 2: also angle from z axis
     *
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintOpDetInfo(Stream&& out, std::string indent = "", unsigned int verbosity = 0) const;

    /**
     * @brief Returns a string with optical detector information
     * @see `PrintOpDetInfo()`
     *
     * Arguments and provided information are the same as in `PrintOpDetInfo()`.
     */
    std::string OpDetInfo(std::string indent = "", unsigned int verbosity = 0) const;

    /// Maximum verbosity supported by `PrintOpDetInfo()`.
    static constexpr unsigned int MaxVerbosity = 2;

  private:
    using LocalTransformation_t =
      LocalTransformationGeo<ROOT::Math::Transform3D, LocalPoint_t, LocalVector_t>;

    LocalTransformation_t fTrans; ///< Optical-detector-to-world transformation.
    const TGeoNode* fOpDetNode;   ///< Pointer to theopdet node
    Point_t fCenter;              ///< Stored geometric center of the optical detector.

    OpDetID fID; ///< Identifier of this optical detector.

    /// Returns the geometry object as `TGeoTube`, `nullptr` if not a tube.
    TGeoTube const* asTube() const { return dynamic_cast<TGeoTube const*>(Shape()); }

    /// Returns the geometry object as `TGeoSphere`, `nullptr` if not a sphere.
    TGeoSphere const* asSphere() const { return dynamic_cast<TGeoSphere const*>(Shape()); }

    /// Returns the geometry object as `TGeoBBox`, `nullptr` if not box-derived.
    TGeoBBox const* asBox() const { return dynamic_cast<TGeoBBox const*>(Shape()); }

  }; // class OpDetGeo

} // namespace geo

//------------------------------------------------------------------------------
//--- template implementation
//---

template <typename ShapeObj>
bool geo::OpDetGeo::isShape() const
{
  static_assert(std::is_base_of_v<TGeoShape, std::decay_t<ShapeObj>>);

  // C++ understanding of the business instead of ROOT's (no strong reason)
  TGeoShape const* shape = Shape(); // needed to convince Clang 7 I really mean it
  return typeid(*shape) == typeid(std::decay_t<ShapeObj>);
}

//------------------------------------------------------------------------------
template <typename ShapeObj>
bool geo::OpDetGeo::isShapeLike() const
{
  static_assert(std::is_base_of_v<TGeoShape, std::decay_t<ShapeObj>>);

  // C++ understanding of the business instead of ROOT's (no strong reason)
  return dynamic_cast<std::decay_t<ShapeObj> const*>(Shape()) != nullptr;
}

//------------------------------------------------------------------------------
template <typename Stream>
void geo::OpDetGeo::PrintOpDetInfo(Stream&& out,
                                   std::string indent /* = "" */,
                                   unsigned int verbosity /* = 0 */
                                   ) const
{

  lar::util::RealComparisons<double> cmp(1e-5);

  //----------------------------------------------------------------------------
  out << "optical detector " << ID() << " centered at " << GetCenter() << " cm";

  if (verbosity-- <= 0) return; // 0

  //----------------------------------------------------------------------------
  if (isTube()) {
    out << ", radius: " << RMax() << " cm";
    if (cmp.nonZero(RMin())) out << " (inner: " << RMin() << " cm)";
    out << ", length: " << Length() << " cm";
  }
  else if (isBar()) {
    out << ", bar size " << Width() << " x " << Height() << " x " << Length() << " cm";
  }
  else if (TGeoSphere const* sphere = asSphere(); sphere) {
    assert(isSphere());
    auto const [th1, th2] = std::minmax({sphere->GetTheta1(), sphere->GetTheta2()});
    out << ", ";
    // some information out of the interface
    if (cmp.zero(th1) && cmp.equal(th2, 180.0))
      out << "spherical";
    else if ((cmp.zero(th1) && cmp.equal(th2, 90.0)) ||
             (cmp.equal(th1, 90.0) && cmp.equal(th2, 180.0))) {
      out << "hemispherical";
    }
    else
      out << "spherical portion (" << th1 << " -> " << th2 << " degree)";
    out << " with external radius " << RMax() << " cm";
  }
  else
    out << ", shape: '" << Shape()->IsA()->GetName() << "'";

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------
  out << ", theta(z): " << ThetaZ() << " rad";

  //  if (verbosity-- <= 0) return; // 2

  //----------------------------------------------------------------------------

} // geo::OpDetGeo::PrintOpDetInfo()

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_OPDETGEO_H
