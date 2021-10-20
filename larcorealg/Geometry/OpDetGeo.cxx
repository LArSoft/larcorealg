////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/OpDetGeo.cxx
/// \brief Encapsulate the geometry of an OpDet
///
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/OpDetGeo.h"

// LArSoft libraries
#include "larcorealg/Geometry/details/OpDetGeoSpecFillers.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::makeFromCoords()
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// Framework libraries
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TGeoNode.h"

// C/C++ standard libraries
#include <cmath>

namespace geo{

  //-----------------------------------------
  OpDetGeo::OpDetGeo(TGeoNode const& node, geo::TransformationMatrix&& trans)
    : fTrans{ makeTrans(std::move(trans)) }
    , fOpDetNode{ &node }
    , fBoxCenter{ toWorldCoords(geo::origin<LocalPoint_t>()) }
    , fShape{ makeShape(*fOpDetNode) }
    , fSpecs{ GetShapeCenter() } // only center initialized so far (legacy)
    {}
  
  //......................................................................

  void OpDetGeo::GetCenter(double* xyz, double offset) const
  {
    geo::vect::fillCoords(xyz, GetCenter() + offset * LengthDir());
  }

  //......................................................................

  double OpDetGeo::RMax() const
  {
    if (TGeoSphere const* sphere = asSphere(); sphere) {
      return sphere->GetRmax();
    }
    else if (TGeoTube const* tube = asTube(); tube) {
      return tube->GetRmax();
    }
    else {
      throw std::bad_cast{};
    }
  }

  //......................................................................

  double OpDetGeo::RMin() const
  {
    if (TGeoSphere const* sphere = asSphere(); sphere) {
      return sphere->GetRmin();
    }
    else if (TGeoTube const* tube = asTube(); tube) {
      return tube->GetRmin();
    }
    else {
      throw std::bad_cast{};
    }
  }

  //......................................................................
  double OpDetGeo::ThetaZ() const
  {
    double angle = std::acos(geo::vect::dot(LengthDir(), geo::Zaxis()));
    if (angle < 0.0) angle += util::pi();
    return angle;
  }

  //......................................................................
  double OpDetGeo::ThetaZ(bool degree) const
    { return degree? util::RadiansToDegrees(ThetaZ()): ThetaZ(); }



  //......................................................................
  double OpDetGeo::DistanceToPoint(geo::Point_t const& point) const
    { return fromCenter(point).R(); }
  double OpDetGeo::DistanceToPoint(double const* xyz) const
    { return DistanceToPoint(geo::vect::makeFromCoords<geo::Point_t>(xyz)); }

  //......................................................................
  std::string OpDetGeo::OpDetInfo
    (std::string indent /* = "" */, unsigned int verbosity /* = 0 */) const
  {
    std::ostringstream sstr;
    PrintOpDetInfo(sstr, indent, verbosity);
    return sstr.str();
  } // OpDetGeo::OpDetInfo()

  //......................................................................
  double OpDetGeo::CosThetaFromNormal(geo::Point_t const& point) const
    { return geo::vect::dot(LengthDir(), fromCenter(point)); }
  double OpDetGeo::CosThetaFromNormal(double const* xyz) const
    { return CosThetaFromNormal(geo::vect::makeFromCoords<geo::Point_t>(xyz)); }

  
  //......................................................................
  void OpDetGeo::UpdateAfterSorting(
    geo::OpDetID opdetid,
    geo::AffinePlaneBase<geo::Vector_t, geo::Point_t> const* directions
  ) {

    fID = opdetid;
    fDirections = directions
      ? directionsFromReference(*directions)
      : standardDirections(fTrans)
      ;
    fSpecs
      = std::visit(details::OpDetGeoSpecFiller{ fTrans, fDirections }, fShape);

  } // OpDetGeo::UpdateAfterSorting()


  //......................................................................
  geo::PlaneBase<geo::Vector_t> OpDetGeo::directionsFromReference
    (geo::AffinePlaneBase<geo::Vector_t, geo::Point_t> const& reference) const
  {
    /*
     * We adopt `reference` directions as our local reference directions,
     * but we require our final reference to have a normal direction pointing
     * toward the target point provided as origin of the reference.
     * Id that is not the case, we flip the normal direction
     * (and the secondary direction to preserve "handness").
     */
    
    bool const flip =
      geo::vect::dot(reference.NormalDir(), fromCenter(reference.Origin()))
      < 0.0
      ;
    
    return geo::PlaneBase<geo::Vector_t>{
      reference.MainDir(),
      flip? -reference.SecondaryDir(): reference.SecondaryDir()
      };
    
  } // OpDetGeo::directionsFromReference()
  
  
  //......................................................................
  geo::PlaneBase<geo::Vector_t> OpDetGeo::standardDirections
    (LocalTransformation_t const& trans)
  {
    return {
      trans.toWorldCoords(geo::Xaxis<LocalVector_t>()), // main (width)
      trans.toWorldCoords(geo::Yaxis<LocalVector_t>())  // secondary (height)
      };
  } // OpDetGeo::standardDirections()
  
  
  //......................................................................
  auto OpDetGeo::makeShape(TGeoNode const& node) -> Shape_t {
    TGeoShape const* shape = node.GetVolume()->GetShape();
    assert(shape); // supposedly ROOT guarantees a shape...
    if (auto obj = dynamic_cast<TGeoSphere const*>(shape)) return { obj };
    if (auto obj = dynamic_cast<TGeoTube const*>(shape)) return { obj };
    if (auto obj = dynamic_cast<TGeoBBox const*>(shape)) return { obj };
    throw cet::exception("Geometry")
      << "geo::OpDetGeo does not support optical detectors of shape '"
      << shape->GetName() << "'\n"
      ;
  } // OpDetGeo::makeShape()
  
  
  //......................................................................

}
////////////////////////////////////////////////////////////////////////
