////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/OpDetGeo.h
/// \brief Encapsulate the geometry of an optical detector
/// \ingroup Geometry
///
/// \author  bjpjones@mit.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_OPDETGEO_H
#define LARCOREALG_GEOMETRY_OPDETGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/details/OpDetGeoSpecs.h"
#include "larcorealg/Geometry/Decomposer.h" // geo::AffinePlaneBase<>
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_optical_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::OpDetID

// ROOT libraries
#include "TGeoMatrix.h" // TGeoHMatrix
#include "TGeoTube.h"
#include "TGeoSphere.h"
#include "TGeoBBox.h"
#include "TClass.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <array>
#include <variant>
#include <algorithm> // std::minmax()
#include <typeinfo> // typeid()
#include <type_traits> // std::decay_t(), std::is_base_of_v
#include <cassert>


// forward declarations
class TGeoNode;

namespace geo {
  
  /**
   * @brief Geometry of the sensitive part of an optical detector.
   * @ingroup Geometry
   * 
   * Note that the definition of "sensitive part" may be loosened for
   * convenience sake. For example, while in a detector made of a silicon
   * photomultiplier (SiPM) with a wavelength guide the SiPM is, strictly,
   * the sensitive part, It may be more convenient to define the guide volume
   * as the sensitive one. Likewise in a Arapuca configuration the window with
   * the dichroic filter is likely the most convenient volume to be defined
   * "sensitive".
   * 
   * The geometry defines a few quantities in a non-trivial way.
   *  * _Normal direction_ is the direction the detector is most sensitive
   *    to; this geometry does not natively support the concept of double-faced
   *    detectors.
   *  * _Center_ is a reference point of the detector; see `GetCenter()` for
   *    the details of its definition.
   *  * _Reference directions_ are directions aligned with the wire planes the
   *    optical detector faces or more in general with the cryostat ones.
   * 
   * A few shapes are explicitly supported.
   *  * _Cylinder_ ("tube"): the detector is treated as a thin disc, reporting as
   *    center the one of the circular face pointed by the normal direction and
   *    verse; a full cylinder is always assumed and its axis is required to be
   *    aligned with the normal direction.
   *  * _Sphere_: complete and incomplete spheres are supported. An incomplete
   *    sphere is a `TGeoSphere` with `GetTheta1()` and `GetTheta2()` not set to
   *    0 and 180 degrees respectively. Only spheres with either one of these
   *    two angles at its limit value are supported. The normal is required to
   *    be aligned with the local _z_ of the sphere (that is the axis the two
   *    theta angles are defined on).
   *  * _Box_: the sides of the box are required to be aligned with the
   *    reference directions.
   *  * _Others_: if the ROOT shape object (`TGeoShape`) is derived from
   *    `TGeoBBox`, the shape is treated as it were that box. Otherwise,
   *    some features of `geo::OpDetGeo` will not be available
   *    Center and normal direction are guaranteed to be available though.
   * 
   * 
   * Initialization
   * ---------------
   * 
   * An `geo::OpDetGeo` is initialized in two steps:
   *  1. construction assigns the aspects that are intrinsic to the detector
   *     itself;
   *  2. an update (`UpdateAfterSorting()`) assigns the aspects which involve
   *     the context surrounding the detector (e.g. ID, directions relative to
   *     the TPC, etc.).
   * 
   * Before the two-step initialization is completed the object should be
   * considered in a non-initialized state and only initialization code (e.g.
   * a `geo::GeometryBuilder` object) should use it while in that state.
   * 
   */
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

    /// Type of points in the local GDML optical detector frame.
    using LocalPoint_t = geo::OpticalPoint_t;

    /// Type of displacement vectors in the local GDML optical detector frame.
    using LocalVector_t = geo::OpticalVector_t;

    ///@}

    /**
     * @brief Constructor: places detector in space and determines its shape.
     * @param node ROOT node representing the sensitive volume of this detector
     * @param trans transformation from local to world coordinates
     * 
     * The specified `node` is used to extract the location of the optical
     * detector (nominal center and size).
     */
    OpDetGeo(TGeoNode const& node, geo::TransformationMatrix&& trans);

    /// Returns the geometry ID of this optical detector.
    geo::OpDetID const& ID() const { return fID; }

    /**
     * @brief Return the center position of the detector.
     * @param[out] xyz array where the position is filled, as `{ x, y, z }` [cm]
     * @param offset (default: 0) shift from the center along `LengthDir()` [cm]
     * 
     * The position is expressed in world coordinates, and in centimeters.
     * 
     * @deprecated unconditionally if `localz` is `0`; and either way the use
     *             of `geo::Point_t` interface (`GetCenter()`) is recommended
     *             instead: `opDet.GetCenter() + offset * opDet.LengthDir()`.
     */
    void GetCenter (double* xyz, double offset = 0.0) const;
    /**
     * @brief Returns the center of the _volume_ of the detector.
     * @return the center of the volume (world coordinates, centimeters)
     * 
     * This point is always within the volume of the detector.
     * It is defined as the geometric center of the solid for boxes, cylinders
     * and for spheres larger than an hemisphere. For sphere sections smaller
     * than an hemisphere, this is instead the center of the flat circular face.
     */
    geo::Point_t const& GetCenter() const { return fSpecs.center; }
    /// Returns the central point of the _surface_ of the detector.
    /// @return the center of the solid (world coordinates, centimeters)
    geo::Point_t GetDetectorTip() const
      { return GetCenter() + fSpecs.centerToTip * LengthDir(); }
    /// Returns the center of the geometric solid the detector is modeled after
    /// @return the center of the solid (world coordinates, centimeters)
    geo::Point_t const& GetShapeCenter() const { return fBoxCenter; }
    
    /**
     * @brief Returns the inner radius of the solid describing the detector.
     * @return the inner radius of the solid [cm]
     * @throw std::bad_cast if the shape does not support this concept
     * 
     * For sphere- and cylinder-based shapes, this method returns the inner
     * radius of the detector shape. For a cylinder (disc) this is likely `0`,
     * while for a sphere it will depend of the thickness of the spherical
     * surface (for example in case of photomultiplier tubes, the difference
     * `Rmax() - Rmin()` may match the thickness of its glass).
     */
    double RMin() const;
    
    /**
     * @brief Returns the outer radius of the solid describing the detector.
     * @return the outer radius of the solid [cm]
     * @throw std::bad_cast if the shape does not support this concept
     * 
     * For sphere- and cylinder-based shapes, this method returns the outer
     * radius of the detector shape, that is the same as `Width()` and
     * `Height()` (except for partial spheres smaller than a hemisphere).
     * For other shapes (box), this method should not be called and it will
     * throw an exception otherwise.
     */
    double RMax() const;
    
    
    // --- BEGIN -- Directions and sizes ---------------------------------------
    /// @name Directions and sizes
    /// @{
    /**
     * @brief Returns the "length" of the detector from its center.
     * @see `Length()`, `LengthDir()`, `HalfW()`, `HalfH()`
     * 
     * This is a rough estimate of the distance between the center of the
     * detector and the border toward the normal direction and verse
     * (`LengthDir()`, also `NormalDir()`).
     * 
     * The details of the computation differ according to the shape of the
     * detector.
     * * _Sphere_: distance from the center to the tip of the sensitive surface;
     *   the center is defined as in `GetCenter()`. Note that this may be
     *   different than half the `Length()`.
     * * _Tube_: distance from the center of the tube to the sensitive disc.
     * * _Box_: distance from the center of the box to the sensitive side.
     */
    double HalfL() const { return Length() / 2.0; }
    
    /**
     * @brief Returns the "width" of the detector from its center.
     * @see `Width()`, `WidthDir()`, `HalfH()`, `HalfL()`
     * 
     * This is a rough estimate of the distance between the center of the
     * detector and the borders along the width direction (`WidthDir()`).
     * Supported shapes (sphere, tube, box) are assumed symmetric on this
     * direction so that it is irrelevant which of the two borders toward
     * the two possible verses of the directions is being considered.
     */
    double HalfW() const { return Width() / 2.0; }
    
    /**
     * @brief Returns the "height" of the detector from its center.
     * @see `Height()`, `HeightDir()`, `HalfW()`, `HalfL()`
     * 
     * This is a rough estimate of the distance between the center of the
     * detector and the borders along the height direction (`HeightDir()`).
     * Supported shapes (sphere, tube, box) are assumed symmetric on this
     * direction so that it is irrelevant which of the two borders toward
     * the two possible verses of the directions is being considered.
     */
    double HalfH() const { return Height() / 2.0; }
    
    /**
     * @brief Returns the full length of the detector [cm]
     * @see `LengthDir()`, `HalfL()`, `Width()`, `Height()`
     * 
     * The length is defined along the "normal" direction (`LengthDir()`).
     * It is defined as twice the half length from `HalfL()`, with the following
     * exception:
     *  * spherical detector: if the sphere is incomplete, the length is the
     *    distance tip-to-tip; for example, in a hemispherical PMT surface,
     *    the `Length()` is the same as the radius (`Rmax()`), and also has the
     *    same value as `HalfL()`.
     */
    double Length() const { return fSpecs.sides.length; }
    
    /**
     * @brief Returns the full width of the detector [cm]
     * @see `WidthDir()`, `HalfW()`, `Height()`, `Length()`
     * 
     * The width is defined along the width direction (`WidthDir()`).
     * It is defined as twice the half width from `HalfW()`:
     *  * Cylindrical detector: diameter of the base.
     *  * Spherical detector: diameter of the sphere, or smaller if the sphere
     *    is partial and smaller than a hemisphere.
     *  * Box: full side on width direction.
     */
    double Width() const { return fSpecs.sides.width; }
    
    /**
     * @brief Returns the full height of the detector [cm]
     * @see `HeightDir()`, `HalfH()`, `Width()`, `Length()`
     * 
     * The height is defined along the width direction (`HeightDir()`).
     * It is defined as twice the half height from `HalfH()`:
     *  * Cylindrical detector: same as `HalfW()`
     *  * Spherical detector: same as `HalfW()`
     *  * Box: full side on height direction.
     */
    double Height() const { return fSpecs.sides.height; }
    
    /**
     * @brief Returns the direction the optical detector is most sensitive to.
     * @returns a unit vector
     * 
     * The normal direction can be thought as the direction that the optical
     * detector faces toward.
     * It may be explicitly provided on construction, or otherwise it is
     * assigned as the completion of a positive base with `WidthDir()` and
     * `HeightDir()`.
     */
    geo::Vector_t NormalDir() const { return fDirections.NormalDir(); }
    
    /// Returns the direction where length is measured (matches `NormalDir()`).
    geo::Vector_t LengthDir() const { return NormalDir(); }
    /**
     * @brief Returns the direction the width is measured on.
     * @returns a unit vector
     * 
     * The width direction may be explicitly provided on construction, or
     * otherwise it is assigned as the transformation in the world frame of the
     * local _x_ coordinate.
     */
    geo::Vector_t WidthDir() const { return fDirections.MainDir(); }
    
    /**
     * @brief Returns the direction the height is measured on.
     * @returns a unit vector
     * 
     * The height direction may be explicitly provided on construction, or
     * otherwise it is assigned as the transformation in the world frame of the
     * local _y_ coordinate.
     */
    geo::Vector_t HeightDir() const { return fDirections.SecondaryDir(); }

    /// @}
    // --- END ---- Directions and sizes ---------------------------------------
    
    /// Returns the vector from the detector center (`GetCenter()`) to `point`.
    geo::Vector_t fromCenter(geo::Point_t const& point) const
      { return point - GetCenter(); }
    
    double ThetaZ() const;  ///< returns angle of detector
                            ///< with respect to z axis
                            ///< in the Y-Z plane, in radians
    double ThetaZ(bool degrees) const; ///< returns angle of detector
                                       ///< with respect to z axis
                                       ///< in the Y-Z plane
    //@{
    /// Get cos(angle) to normal of this detector - used for solid angle calcs
    double CosThetaFromNormal(geo::Point_t const& point) const;
    double CosThetaFromNormal(double const* xyz) const;
    //@}
    //@{
    /// Returns the distance of the specified point from detector center [cm]
    double DistanceToPoint(geo::Point_t const& point) const;
    double DistanceToPoint(double const* xyz) const;
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
    void LocalToWorld(const double* opdet, double* world) const
      { fTrans.LocalToWorld(opdet, world); }

    /// Transform point from local optical detector frame to world frame.
    geo::Point_t toWorldCoords(LocalPoint_t const& local) const
      { return fTrans.toWorldCoords(local); }

    /// Transform direction vector from local to world.
    void LocalToWorldVect(const double* opdet, double* world) const
      { fTrans.LocalToWorldVect(opdet, world); }

    /// Transform direction vector from local to world.
    geo::Vector_t toWorldCoords(LocalVector_t const& local) const
      { return fTrans.toWorldCoords(local); }

    /// Transform point from world frame to local optical detector frame.
    void WorldToLocal(const double* world, double* opdet) const
      { fTrans.WorldToLocal(world, opdet); }

    /// Transform point from world frame to local optical detector frame.
    LocalPoint_t toLocalCoords(geo::Point_t const& world) const
      { return fTrans.toLocalCoords(world); }

    /// Transform direction vector from world to local.
    void WorldToLocalVect(const double* world, double* opdet) const
      { fTrans.WorldToLocalVect(world, opdet); }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(geo::Vector_t const& world) const
      { return fTrans.toLocalCoords(world); }

    /// @}

    /// Returns the ROOT object describing the detector geometry.
    const TGeoNode*     Node() const { return fOpDetNode; }

    
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

    /**
     * @brief Performs all updates after cryostat has sorted optical detectors.
     * @param opdetid the ID to be assigned to this optical detector
     * @param directions reference directions for the orientation of detector
     * 
     * The content of `directions` is directly used to "orient" the optical
     * detector: `directions.MainDir()` becomes the width direction, and
     * `directions.SecondaryDir()` becomes the height direction.
     * If `directions` is `nullptr`, directions from `standardDirections()` are
     * used.
     */
    void UpdateAfterSorting(
      geo::OpDetID opdetid,
      geo::AffinePlaneBase<geo::Vector_t, geo::Point_t> const* directions
      );


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
    void PrintOpDetInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 0) const;

    /**
     * @brief Returns a string with optical detector information
     * @see `PrintOpDetInfo()`
     *
     * Arguments and provided information are the same as in `PrintOpDetInfo()`.
     */
    std::string OpDetInfo
      (std::string indent = "", unsigned int verbosity = 0) const;

    /// Maximum verbosity supported by `PrintOpDetInfo()`.
    static constexpr unsigned int MaxVerbosity = 2;

  private:
    /// Type for local-to-world transformations.
    using LocalTransformation_t = details::OpDetGeoLocalTransformation;

    /// Type holding the shape of the detector.
    using Shape_t = std::variant<
//       TGeoShape const*,
      TGeoSphere const*,
      TGeoTube const*,
      TGeoBBox const*
      >;
    
    using Specs_t = details::OpDetGeoSpecs;
    
    LocalTransformation_t fTrans; ///< Optical-detector-to-world transformation.
    const TGeoNode* fOpDetNode;  ///< Pointer to theopdet node
    
    /// Base directions of the detector geometry.
    geo::PlaneBase<geo::Vector_t> fDirections;
    
    Shape_t fShape; ///< Shape of the detector as `TGeoShape` object.
    
    Specs_t fSpecs; ///< Aggregate geometry information from the detector.
    
    /// Geometric center of the solid the optical detector is modeled after.
    geo::Point_t fBoxCenter;

    geo::OpDetID fID; ///< Identifier of this optical detector.

    /// Returns the geometry object as `TGeoTube`, `nullptr` if not a tube.
    TGeoTube const* asTube() const
      { return dynamic_cast<TGeoTube const*>(Shape()); }

    /// Returns the geometry object as `TGeoSphere`, `nullptr` if not a sphere.
    TGeoSphere const* asSphere() const
      { return dynamic_cast<TGeoSphere const*>(Shape()); }

    /// Returns the geometry object as `TGeoBBox`, `nullptr` if not box-derived.
    TGeoBBox const* asBox() const
      { return dynamic_cast<TGeoBBox const*>(Shape()); }
    
    /**
     * @brief Returns reference directions for this optical detector.
     * @param reference the candidate reference directions, with as origin
     *                  the location the detector should point toward
     * @return the reference directions for this object
     */
    geo::PlaneBase<geo::Vector_t> directionsFromReference
      (geo::AffinePlaneBase<geo::Vector_t, geo::Point_t> const& reference)
      const;
    
    
    /// Builds a shape object variant out of the specified node.
    /// @throws cet::exception (category: `Geometry`) if unsupported shape.
    static Shape_t makeShape(TGeoNode const& node);
    
    /// Turns a `geo::TransformationMatrix` into a `LocalTransformation_t`.
    static LocalTransformation_t makeTrans(geo::TransformationMatrix&& trans)
      { return { std::move(trans) }; }
    
    /// Returns "standard" directions (local x/y turned into world frame).
    static geo::PlaneBase<geo::Vector_t> standardDirections
      (LocalTransformation_t const& trans);

  }; // class OpDetGeo

} // namespace geo


//------------------------------------------------------------------------------
//--- template implementation
//---

template <typename ShapeObj>
bool geo::OpDetGeo::isShape() const {
  static_assert(std::is_base_of_v<TGeoShape, std::decay_t<ShapeObj>>);
  
  // C++ understanding of the business instead of ROOT's (no strong reason)
  TGeoShape const* shape = Shape(); // needed to convince Clang 7 I really mean it
  return typeid(*shape) == typeid(std::decay_t<ShapeObj>);
  
} // geo::OpDetGeo::isShape()

//------------------------------------------------------------------------------
template <typename ShapeObj>
bool geo::OpDetGeo::isShapeLike() const {
  static_assert(std::is_base_of_v<TGeoShape, std::decay_t<ShapeObj>>);
  
  // C++ understanding of the business instead of ROOT's (no strong reason)
  return dynamic_cast<std::decay_t<ShapeObj> const*>(Shape()) != nullptr;
  
} // geo::OpDetGeo::isShapeLike()


//------------------------------------------------------------------------------
template <typename Stream>
void geo::OpDetGeo::PrintOpDetInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 0 */
) const {
  
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
    out << ", bar size " << Width() << " x " << Height() << " x " << Length()
      << " cm";
  }
  else if (TGeoSphere const* sphere = asSphere(); sphere) {
    assert(isSphere());
    auto const [ th1, th2 ]
      = std::minmax({ sphere->GetTheta1(), sphere->GetTheta2() });
    out << ", ";
    // some information out of the interface
    if (cmp.zero(th1) && cmp.equal(th2, 180.0)) out << "spherical";
    else if ((cmp.zero(th1) && cmp.equal(th2, 90.0))
      || (cmp.equal(th1, 90.0) && cmp.equal(th2, 180.0)))
    {
      out << "hemispherical";
    }
    else out << "spherical portion (" << th1 << " -> " << th2 << " degree)";
    out << " with external radius " << RMax() << " cm";
  }
  else out << ", shape: '" << Shape()->IsA()->GetName() << "'";

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------
  out << ", theta(z): " << ThetaZ() << " rad";

//  if (verbosity-- <= 0) return; // 2

  //----------------------------------------------------------------------------

} // geo::OpDetGeo::PrintOpDetInfo()

//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_OPDETGEO_H
