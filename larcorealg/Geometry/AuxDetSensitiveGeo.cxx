////////////////////////////////////////////////////////////////////////
/// @file  larcorealg/Geometry/AuxDetSensitiveGeo.cxx
/// @brief Encapsulate the geometry of the sensitive portion of an auxilary detector
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect namespace

/// ROOT libraries
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TGeoTrd2.h"
#include "TGeoVolume.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard libraries
#include <sstream> // std::ostringstream

namespace geo {

  //-----------------------------------------
  AuxDetSensitiveGeo::AuxDetSensitiveGeo(TGeoNode const* node, TransformationMatrix&& trans)
    : fTrans(std::move(trans)), fTotalVolume(node->GetVolume())
  {

    MF_LOG_DEBUG("Geometry") << "detector sensitive total  volume is " << fTotalVolume->GetName();

    InitShapeSize();
  }

  //......................................................................
  Point_t AuxDetSensitiveGeo::GetCenter(double localz /* = 0.0 */) const
  {
    return toWorldCoords(LocalPoint_t{0.0, 0.0, localz});
  }

  //......................................................................

  // Return the unit normal vector (0,0,1) in local coordinates to global coordinates
  Vector_t AuxDetSensitiveGeo::GetNormalVector() const
  {
    return toWorldCoords(Zaxis<LocalVector_t>());
  }

  //......................................................................
  Length_t AuxDetSensitiveGeo::DistanceToPoint(double const* point) const
  {
    return DistanceToPoint(vect::makePointFromCoords(point));
  }

  //......................................................................
  std::string AuxDetSensitiveGeo::AuxDetInfo(std::string indent /* = "" */,
                                             unsigned int verbosity /* = 1 */) const
  {
    std::ostringstream sstr;
    PrintAuxDetInfo(sstr, indent, verbosity);
    return sstr.str();
  }

  //......................................................................
  void AuxDetSensitiveGeo::InitShapeSize()
  {
    // set the ends depending on whether the shape is a box or trapezoid
    std::string volName(fTotalVolume->GetName());
    if (volName.find("Trap") != std::string::npos) {

      //       Small Width
      //          ____          Height is the thickness
      //         /    \     T     of the trapezoid
      //        /      \    |
      //       /        \   | Length
      //      /__________\  _
      //         Width
      fHalfHeight = ((TGeoTrd2*)fTotalVolume->GetShape())->GetDy1(); // same as Dy2()
      fLength = 2.0 * ((TGeoTrd2*)fTotalVolume->GetShape())->GetDz();
      fHalfWidth1 = ((TGeoTrd2*)fTotalVolume->GetShape())->GetDx1(); // at -Dz
      fHalfWidth2 = ((TGeoTrd2*)fTotalVolume->GetShape())->GetDx2(); // at +Dz
    }
    else {
      fHalfWidth1 = ((TGeoBBox*)fTotalVolume->GetShape())->GetDX();
      fHalfHeight = ((TGeoBBox*)fTotalVolume->GetShape())->GetDY();
      fLength = 2.0 * ((TGeoBBox*)fTotalVolume->GetShape())->GetDZ();
      fHalfWidth2 = fHalfWidth1;
    }
  }
}
////////////////////////////////////////////////////////////////////////
