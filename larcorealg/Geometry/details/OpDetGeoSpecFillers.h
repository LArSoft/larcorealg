/**
 * @file   larcorealg/Geometry/details/OpDetGeoSpecFillers.h
 * @brief  Definition of internal data of `geo::OpDetGeo`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    larcorealg/Geometry/details/OpDetGeoSpecFillers.cxx
 * 
 */

#ifndef LARCOREALG_GEOMETRY_DETAILS_OPDETGEOSPECFILLERS_H
#define LARCOREALG_GEOMETRY_DETAILS_OPDETGEOSPECFILLERS_H

// LArSoft libraries
#include "larcorealg/Geometry/details/OpDetGeoSpecs.h"
#include "larcorealg/Geometry/Decomposer.h" // geo::PlaneBase<>
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t...

// Framework libraries
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace geo::details {
  
  struct OpDetGeoSpecFillerBase;
  struct OpDetGeoSpecFillerFromTube;
  struct OpDetGeoSpecFillerFromSphere;
  struct OpDetGeoSpecFillerFromBox;
  
  struct OpDetGeoSpecFiller;
  
} // namespace geo::details

// foreign forward declarations
class TGeoTube;
class TGeoSphere;
class TGeoBBox;


// -----------------------------------------------------------------------------
struct geo::details::OpDetGeoSpecFillerBase {
  
  // no serviceable parts here
  
    protected:
  friend class OpDetGeoSpecFiller;
  
  /// Type for local transformation. 
  using LocalTransformation_t = OpDetGeoLocalTransformation;
  using LocalVector_t = LocalTransformation_t::LocalVector_t;
  using LocalPoint_t = LocalTransformation_t::LocalPoint_t;
  
  using Specs_t = OpDetGeoSpecs; ///< Type of specifications to be filled.
  
  LocalTransformation_t const& trans; ///< Local-to-world transformation.
  geo::PlaneBase<geo::Vector_t> const& dirs; ///< Reference directions.
  
  OpDetGeoSpecFillerBase(
    LocalTransformation_t const& trans,
    geo::PlaneBase<geo::Vector_t> const& directions
    );
  
}; // geo::details::OpDetGeoSpecFillerBase


// -----------------------------------------------------------------------------
struct geo::details::OpDetGeoSpecFillerFromTube: OpDetGeoSpecFillerBase {
  
  using OpDetGeoSpecFillerBase::OpDetGeoSpecFillerBase; // inherit constructor
  
  Specs_t operator()(TGeoTube const& disc) const;
  
    private:
  
  void checkAlignment(TGeoTube const& disc) const;
  double extractRadius(TGeoTube const& disc) const noexcept;
  double extractLength(TGeoTube const& disc) const noexcept;
  geo::Point_t extractCenter(TGeoTube const&) const noexcept;
  
}; // geo::details::OpDetGeoSpecFillerFromTube


// -----------------------------------------------------------------------------
struct geo::details::OpDetGeoSpecFillerFromSphere: OpDetGeoSpecFillerBase {

  using OpDetGeoSpecFillerBase::OpDetGeoSpecFillerBase; // inherit constructor
  
  Specs_t operator()(TGeoSphere const& sphere) const;
  
    private:
  
  /// @throws cet::exception (category `Geometry`) on failure check
  /// @return alignment of local z axis and normal (+1.0/-1.0)
  double checkShape(TGeoSphere const& sphere) const;
  
  double extractLength
    (TGeoSphere const& sphere, double alignment) const noexcept;
  
  double extractRadius(TGeoSphere const& sphere, double length) const noexcept;
  
  double extractDistanceToTip
    (TGeoSphere const& sphere, double length) const noexcept;
  
  geo::Point_t solidCenter() const noexcept;
  
  geo::Point_t extractCenter
    (TGeoSphere const& sphere, double length) const noexcept;
  
}; // geo::details::OpDetGeoSpecFillerFromSphere


// -----------------------------------------------------------------------------
struct geo::details::OpDetGeoSpecFillerFromBox: OpDetGeoSpecFillerBase {
  
  using OpDetGeoSpecFillerBase::OpDetGeoSpecFillerBase; // inherit constructor
  
  Specs_t operator()(TGeoBBox const& box) const;
  
    private:
  
  struct SideInfo_t {
    LocalVector_t localDir;
    geo::Vector_t dir;
    double size;
  };
  
  Specs_t::BoxSize_t extractSides(TGeoBBox const& box) const;
  
  std::array<SideInfo_t, 3U> makeSideInfo(TGeoBBox const& box) const;
  
  SideInfo_t const* matchDir
    (std::array<SideInfo_t const*, 3U>& sides, geo::Vector_t const& dir) const;
  
  void validateAlignment(
    SideInfo_t const* widthInfo,
    SideInfo_t const* heightInfo,
    SideInfo_t const* lengthInfo,
    std::array<SideInfo_t, 3U> const& sideInfo
    ) const;
  
  geo::Point_t extractCenter(TGeoBBox const& box) const noexcept;
  
}; // geo::OpDetGeo::OpDetGeoSpecFillerFromBox


// -----------------------------------------------------------------------------
/// This is an "overloader": provides a `operator()` for each type in `Shape_t`.
struct geo::details::OpDetGeoSpecFiller: OpDetGeoSpecFillerBase {
  
  OpDetGeoSpecFiller(
    LocalTransformation_t const& trans,
    geo::PlaneBase<geo::Vector_t> const& directions
    );
  
  Specs_t operator()(TGeoTube const* node) const;
  Specs_t operator()(TGeoSphere const* node) const;
  Specs_t operator()(TGeoBBox const* node) const;
  
}; // geo::details::OpDetGeoSpecFiller()


// -----------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_DETAILS_OPDETGEOSPECFILLERS_H
