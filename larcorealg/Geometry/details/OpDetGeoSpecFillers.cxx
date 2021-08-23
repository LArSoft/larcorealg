/**
 * @file   larcorealg/Geometry/details/OpDetGeoSpecFillers.cxx
 * @brief  Algorithms for optical detector specifications (implementations).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    larcorealg/Geometry/details/OpDetGeoSpecFillers.h
 * 
 */


// class header
#include "larcorealg/Geometry/details/OpDetGeoSpecFillers.h"

// ROOT libraries
#include "TGeoTube.h"
#include "TGeoSphere.h"
#include "TGeoBBox.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::DegreesToRadians()

#if 0
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::makeFromCoords()

// Framework libraries
#include "cetlib_except/exception.h"

#endif // 0

// C/C++ standard libraries
#include <array>
#include <cmath> // std::abs(), std::sqrt(), ...
#include <cassert>


// -----------------------------------------------------------------------------
// ---  geo::details::OpDetGeoSpecFillerBase
// -----------------------------------------------------------------------------
geo::details::OpDetGeoSpecFillerBase::OpDetGeoSpecFillerBase(
  LocalTransformation_t const& trans,
  geo::PlaneBase<geo::Vector_t> const& directions
  )
  : trans{ trans }, dirs{ directions }
{}


// -----------------------------------------------------------------------------
// ---  geo::details::OpDetGeoSpecFillerFromTube
// -----------------------------------------------------------------------------
auto geo::details::OpDetGeoSpecFillerFromTube::operator()
  (TGeoTube const& disc) const -> Specs_t 
{
  checkAlignment(disc);
  double const radius = extractRadius(disc);
  double const length = extractLength(disc);
  return { // C++20: specify data member names explicitly
      extractCenter(disc)                   // center
    , dirs.NormalDir()                      // normalDir
    , {                                     // sides
          radius                            //  .width
        , radius                            //  .height
        , length                            //  .length
      }
    , length / 2.0                          // centerToTip
  };
} // geo::details::OpDetGeoSpecFillerFromTube::operator()


// -----------------------------------------------------------------------------
void geo::details::OpDetGeoSpecFillerFromTube::checkAlignment
  (TGeoTube const& disc) const
{
  constexpr lar::util::Vector3DComparison cmp { 1e-5 };
  
  geo::Vector_t const& normalDir { dirs.NormalDir() };
  geo::Vector_t const shapeZdir
    = trans.toWorldCoords(geo::Zaxis<LocalVector_t>());
  if (cmp.nonParallel(shapeZdir, normalDir)) {
    throw cet::exception{ "Geometry" }
      << "geo::details (disc) at " << extractCenter(disc)
      << " has axis " << shapeZdir
      << " while it should match the provided normal " << normalDir
      << ".\n";
  }
} // geo::details::OpDetGeoSpecFillerFromTube::checkAlignment()


// -----------------------------------------------------------------------------
double geo::details::OpDetGeoSpecFillerFromTube::extractRadius
  (TGeoTube const& disc) const noexcept
  { return disc.GetRmax(); }


// -----------------------------------------------------------------------------
double geo::details::OpDetGeoSpecFillerFromTube::extractLength
  (TGeoTube const& disc) const noexcept
  { return 2.0 * disc.GetDZ(); }


// -----------------------------------------------------------------------------
geo::Point_t geo::details::OpDetGeoSpecFillerFromTube::extractCenter
  (TGeoTube const&) const noexcept
  { return { trans.toWorldCoords(geo::origin<LocalPoint_t>()) }; }


// -----------------------------------------------------------------------------
// ---  geo::details::OpDetGeoSpecFillerFromSphere
// -----------------------------------------------------------------------------
auto geo::details::OpDetGeoSpecFillerFromSphere::operator()
  (TGeoSphere const& sphere) const -> Specs_t
{
  // alignment of local z axis and normal: +1 = parallel, -1 = antiparallel
  double const alignment = checkShape(sphere);
  double const length = extractLength(sphere, alignment);
  double const radius = extractRadius(sphere, length);
  return { // C++20: specify data member names explicitly
      extractCenter(sphere, length)         // center
    , dirs.NormalDir()                      // normalDir
    , {                                     // sides
          radius                            //  .width
        , radius                            //  .height
        , length                            //  .length
      }
    , extractDistanceToTip(sphere, length)  // centerToTip
  };
} // geo::details::OpDetGeoSpecFillerFromSphere::operator()


// -----------------------------------------------------------------------------
double geo::details::OpDetGeoSpecFillerFromSphere::checkShape
  (TGeoSphere const& sphere) const
{
  /*
    * 1. if full sphere, all good already
    * 2. if partial sphere, local z axis parallel to normal direction
    * 3. if one angle is closed, the verse of normal must point to it
    */
  constexpr lar::util::RealComparisons<double> cmp { 1e-5 };
  
  geo::Vector_t const& normalDir { dirs.NormalDir() };
  geo::Vector_t const shapeZdir
    = trans.toWorldCoords(geo::Zaxis<LocalVector_t>());
  double const dotProd = geo::vect::dot(shapeZdir, normalDir);
  
  // 1.
  // GetTheta1() is the smaller angle with respect to positive local z axis.
  if (cmp.zero(sphere.GetTheta1()) && cmp.zero(180.0 - sphere.GetTheta2()))
    return dotProd;
  
  // 2.
  if (!cmp.equal(std::abs(dotProd), 1.0)) {
    throw cet::exception{ "Geometry" }
      << "geo::details (incomplete sphere) at " << solidCenter()
      << " has axis " << shapeZdir
      << " while it should match the provided normal " << normalDir
      << " (dot product: " << dotProd << ").\n";
  }
  
  // 3.
  double const frontAngle
    { (dotProd > 0)? sphere.GetTheta1(): 180.0 - sphere.GetTheta2() };
  if (cmp.nonZero(frontAngle)) {
    throw cet::exception{ "Geometry" }
      << "geo::details (incomplete sphere) at " << solidCenter()
      << " has the open side facing the provided normal " << normalDir
      << " (detector axis \"z\": " << shapeZdir
      << "; theta: (+z)=" << sphere.GetTheta1()
      << ", (-z)=" << sphere.GetTheta2()
      << ").\n";
  }
  
  return dotProd;
} // geo::details::OpDetGeoSpecFillerFromSphere::checkShape()


// -----------------------------------------------------------------------------
double geo::details::OpDetGeoSpecFillerFromSphere::extractLength
  (TGeoSphere const& sphere, double alignment) const noexcept
{
  double const backAngle
    { (alignment < 0)? 180.0 - sphere.GetTheta1(): sphere.GetTheta2() };
  return
    sphere.GetRmax() * (1.0 - std::cos(util::DegreesToRadians(backAngle)));
} // geo::details::OpDetGeoSpecFillerFromSphere::extractLength()


// -----------------------------------------------------------------------------
double geo::details::OpDetGeoSpecFillerFromSphere::extractRadius
  (TGeoSphere const& sphere, double length) const noexcept
{
  double const R = sphere.GetRmax();
  return (length >= R)? R: std::sqrt( (2.0 * R - length) * length );
}


// -----------------------------------------------------------------------------
double geo::details::OpDetGeoSpecFillerFromSphere::extractDistanceToTip
  (TGeoSphere const& sphere, double length) const noexcept
  { return std::min(length, sphere.GetRmax()); }


// -----------------------------------------------------------------------------
geo::Point_t
geo::details::OpDetGeoSpecFillerFromSphere::solidCenter() const noexcept
  { return trans.toWorldCoords(geo::origin<LocalPoint_t>()); }


// -----------------------------------------------------------------------------
geo::Point_t geo::details::OpDetGeoSpecFillerFromSphere::extractCenter
  (TGeoSphere const& sphere, double length) const noexcept
{
  geo::Point_t center { solidCenter() };
  double const R = sphere.GetRmax();
  if (length < R) {
    geo::Vector_t const& normalDir { dirs.NormalDir() };
    center += (R - length) * normalDir;
  }
  return center;
} // geo::details::OpDetGeoSpecFillerFromSphere::extractCenter()


// -----------------------------------------------------------------------------
// ---  geo::details::OpDetGeoSpecFillerFromBox
// -----------------------------------------------------------------------------
auto geo::details::OpDetGeoSpecFillerFromBox::operator()
  (TGeoBBox const& box) const -> Specs_t
{
  Specs_t::BoxSize_t const sides = extractSides(box);
  return { // C++20: specify data member names explicitly
      extractCenter(box)   // center
    , dirs.NormalDir()     // normalDir
    , sides                // sides
    , (sides.length / 2.0) // centerToTip
  };
} // geo::details::OpDetGeoSpecFillerFromBox::operator()


// -----------------------------------------------------------------------------
auto geo::details::OpDetGeoSpecFillerFromBox::extractSides
  (TGeoBBox const& box) const -> Specs_t::BoxSize_t
{
  /*
    * 1. determine the local axes of the box in world coordinates
    * 2. match with the proposed directions
    * 3. validate the alignment
    * 4. fill the side sizes
    */
  
  // 1.
  std::array<SideInfo_t, 3U> const sideInfo { makeSideInfo(box) };
  
  // 2.
  std::array<SideInfo_t const*, 3U> availSideInfo
    { &(sideInfo[0]), &(sideInfo[1]), &(sideInfo[2]) };
  
  SideInfo_t const* widthInfo  = matchDir(availSideInfo, dirs.MainDir());
  SideInfo_t const* heightInfo = matchDir(availSideInfo, dirs.SecondaryDir());
  SideInfo_t const* lengthInfo = matchDir(availSideInfo, dirs.NormalDir());
  
  // 3.
  validateAlignment(widthInfo, heightInfo, lengthInfo, sideInfo);
  
  // 4.
  auto sizeOf = [](SideInfo_t const* info){ return info? info->size: 0.0; };
  return { // C++20: spell out data member names
      sizeOf(widthInfo)   // width
    , sizeOf(heightInfo)  // height
    , sizeOf(lengthInfo)  // length
    };
  
} // geo::details::OpDetGeoSpecFillerFromBox::extractSides()


// -----------------------------------------------------------------------------
auto geo::details::OpDetGeoSpecFillerFromBox::makeSideInfo
  (TGeoBBox const& box) const -> std::array<SideInfo_t, 3U> 
{
  return { // C++20: explicitly spell the data member names out
    SideInfo_t{
        geo::Xaxis<LocalVector_t>()                       // localDir
      , trans.toWorldCoords(geo::Xaxis<LocalVector_t>())  // dir
      , box.GetDX() * 2.0                                 // size
    },
    SideInfo_t{
        geo::Yaxis<LocalVector_t>()                       // localDir
      , trans.toWorldCoords(geo::Yaxis<LocalVector_t>())  // dir
      , box.GetDY() * 2.0                                 // size
    },
    SideInfo_t{
        geo::Zaxis<LocalVector_t>()                       // localDir
      , trans.toWorldCoords(geo::Zaxis<LocalVector_t>())  // dir
      , box.GetDZ() * 2.0                                 // size
    }
    };
} // geo::details::OpDetGeoSpecFillerFromBox::makeSideInfo()


// -----------------------------------------------------------------------------
auto geo::details::OpDetGeoSpecFillerFromBox::matchDir
  (std::array<SideInfo_t const*, 3U>& sides, geo::Vector_t const& dir) const
  -> SideInfo_t const* 
{
  constexpr lar::util::Vector3DComparison cmp
    { lar::util::RealComparisons<double>{ 1e-3 } };
  
  SideInfo_t const* matched = nullptr;
  for (SideInfo_t const*& side: sides) {
    if (!side) continue;
    if (cmp.nonParallel(side->dir, dir)) continue;
    std::swap(matched, side);
    break;
  }
  return matched;
} // geo::details::OpDetGeoSpecFillerFromBox::matchDir()


// -----------------------------------------------------------------------------
void geo::details::OpDetGeoSpecFillerFromBox::validateAlignment(
  SideInfo_t const* widthInfo,
  SideInfo_t const* heightInfo,
  SideInfo_t const* lengthInfo,
  std::array<SideInfo_t, 3U> const& sideInfo
) const {
  unsigned int nErrors { 0 };
  for (SideInfo_t const* info: { widthInfo, heightInfo, lengthInfo })
    if (!info) ++nErrors;
  if (nErrors == 0) return;
  
  cet::exception e { "Geometry" };
  auto printInfo = [&e](SideInfo_t const* info)
    {
      if (info) {
        e << info->dir << " ( <- " << info->localDir << ", "
          << info->size << " cm)";
      }
      else e << " unmatched!";
    };
  e << "geo::details (incomplete sphere): failed to match " << nErrors
    << " directions of optical detector to the proposed frame:";
  e << "\n   W: " << dirs.MainDir() << " <=> ";
  printInfo(widthInfo);
  e << "\n   H: " << dirs.SecondaryDir() << " <=> ";
  printInfo(heightInfo);
  e << "\n   L: " << dirs.NormalDir() << " <=> ";
  printInfo(lengthInfo);
  throw e;
} // geo::details::OpDetGeoSpecFillerFromBox::validateAlignment()


// -----------------------------------------------------------------------------
geo::Point_t geo::details::OpDetGeoSpecFillerFromBox::extractCenter
  (TGeoBBox const& box) const noexcept
  { return trans.toWorldCoords(geo::origin<LocalPoint_t>()); }


// -----------------------------------------------------------------------------
// ---  geo::details::OpDetGeoSpecFiller
// -----------------------------------------------------------------------------
geo::details::OpDetGeoSpecFiller::OpDetGeoSpecFiller(
  LocalTransformation_t const& trans,
  geo::PlaneBase<geo::Vector_t> const& directions
  )
  : OpDetGeoSpecFillerBase(trans, directions)
  {}


// -----------------------------------------------------------------------------
auto geo::details::OpDetGeoSpecFiller::operator()
  (TGeoTube const* node) const -> Specs_t
  { assert(node); return OpDetGeoSpecFillerFromTube{ trans, dirs }(*node); }
auto geo::details::OpDetGeoSpecFiller::operator()
  (TGeoSphere const* node) const -> Specs_t
  { assert(node); return OpDetGeoSpecFillerFromSphere{ trans, dirs }(*node); }
auto geo::details::OpDetGeoSpecFiller::operator()
  (TGeoBBox const* node) const -> Specs_t
  { assert(node); return OpDetGeoSpecFillerFromBox{ trans, dirs }(*node); }


// -----------------------------------------------------------------------------
