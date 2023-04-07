////////////////////////////////////////////////////////////////////////
/// \file TPCGeo.cxx
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/TPCGeo.h"

// LArSoft includes
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::fillCoords()

// Framework includes
#include "cetlib/container_algorithms.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"

// C/C++ standard libraries
#include <algorithm> // std::max(), std::copy()
#include <cassert>
#include <cmath>
#include <functional> // std::mem_fn()
#include <iterator>   // std::inserter()
#include <map>
#include <sstream> // std::ostringstream

namespace {
  geo::Vector_t to_vector(geo::DriftAxis const driftAxis)
  {
    auto const [axis, sign] = driftAxis;
    switch (axis) {
    case geo::Coordinate::X: return geo::Xaxis() * to_int(sign);
    case geo::Coordinate::Y: return geo::Yaxis() * to_int(sign);
    case geo::Coordinate::Z: return geo::Zaxis() * to_int(sign);
    }
    return {}; // unreachable
  }
}

namespace geo {

  //......................................................................
  TPCGeo::TPCGeo(TGeoNode const* node,
                 TransformationMatrix&& trans,
                 DriftAxis const driftAxis,
                 double const driftDistance)
    : fTrans(std::move(trans))
    , fDriftAxis{driftAxis}
    , fDriftDir{to_vector(driftAxis)}
    , fDriftDistance{driftDistance}
  {
    // all planes are going to be contained in the volume named volTPC
    // now get the total volume of the TPC
    TGeoVolume* vc = node->GetVolume();
    if (!vc) {
      throw cet::exception("Geometry")
        << "cannot find detector outline volume - bail ungracefully\n";
    }

    fTotalVolume = vc;

    auto* active_node = NodeForActiveVolume(node);
    // compute the active volume transformation too
    auto ActiveHMatrix(fTrans.Matrix());
    if (active_node) { ActiveHMatrix *= makeTransformationMatrix(*active_node->GetMatrix()); }

    fActiveCenter = vect::toPoint(ActiveHMatrix(ROOT::Math::Transform3D::Point{}));
    fActiveVolume = active_node->GetVolume();

    MF_LOG_DEBUG("Geometry") << "detector total  volume is " << fTotalVolume->GetName()
                             << "\ndetector active volume is " << fActiveVolume->GetName();

    // set the width, height, and lengths
    auto const* active_volume_box = static_cast<TGeoBBox const*>(fActiveVolume->GetShape());
    fActiveHalfWidth = active_volume_box->GetDX();
    fActiveHalfHeight = active_volume_box->GetDY();
    fActiveLength = 2.0 * active_volume_box->GetDZ();

    fHalfWidth = active_volume_box->GetDX();
    fHalfHeight = active_volume_box->GetDY();
    fLength = 2.0 * active_volume_box->GetDZ();

    // Check that the rotation matrix to the world is the identity, if not we need to
    // change the width, height and length values; the correspondence of these to x, y and
    // z are not guaranteed to be trivial, so we store the two independently (cartesian
    // dimensions in the bounding boxes, the sizes in data members directly).

    // TODO: there must be a more general way to do this...
    double Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz;
    fTrans.Matrix().Rotation().GetComponents(Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz);
    auto const* total_volume_box = static_cast<TGeoBBox const*>(fTotalVolume->GetShape());
    if (Rxx != 1) {
      if (std::abs(Rxz) == 1) {
        fActiveHalfWidth = active_volume_box->GetDZ();
        fHalfWidth = total_volume_box->GetDZ();
        fWidthDir = Zaxis();
      }
      if (std::abs(Rxy) == 1) {
        fActiveHalfWidth = active_volume_box->GetDY();
        fHalfWidth = total_volume_box->GetDY();
        fWidthDir = Yaxis();
      }
    }
    if (Ryy != 1) {
      if (std::abs(Rxy) == 1) {
        fActiveHalfHeight = active_volume_box->GetDX();
        fHalfHeight = total_volume_box->GetDX();
        fHeightDir = Xaxis();
      }
      if (std::abs(Rzy) == 1) {
        fActiveHalfHeight = active_volume_box->GetDZ();
        fHalfHeight = total_volume_box->GetDZ();
        fHeightDir = Zaxis();
      }
    }
    if (Rzz != 1) {
      if (std::abs(Rzx) == 1) {
        fActiveLength = 2. * active_volume_box->GetDX();
        fLength = 2. * total_volume_box->GetDX();
        fLengthDir = Xaxis();
      }
      if (std::abs(Ryz) == 1) {
        fActiveLength = 2. * active_volume_box->GetDY();
        fLength = 2. * total_volume_box->GetDY();
        fLengthDir = Yaxis();
      }
    }

    InitTPCBoundaries();
  }

  //......................................................................
  Point_t TPCGeo::GetCathodeCenter() const
  {
    // 1. find the center of the face of the TPC opposite to the anode
    // 2. compute the distance of it from the last wire plane

    //
    // find the cathode center
    //
    Point_t cathodeCenter = GetActiveVolumeCenter();
    auto const [axis, sign] = DriftAxisWithSign();
    switch (axis) {
    case Coordinate::X:
      cathodeCenter.SetX(cathodeCenter.X() - to_int(sign) * ActiveHalfWidth());
      break;
    case Coordinate::Y:
      cathodeCenter.SetY(cathodeCenter.Y() - to_int(sign) * ActiveHalfHeight());
      break;
    case Coordinate::Z:
      cathodeCenter.SetZ(cathodeCenter.Z() - to_int(sign) * ActiveLength() / 2.0);
      break;
    }

    return cathodeCenter;
  }

  //......................................................................
  double TPCGeo::DriftDistance() const { return fDriftDistance; }

  //......................................................................
  std::string TPCGeo::TPCInfo(std::string indent /* = "" */, unsigned int verbosity /* = 1 */) const
  {
    std::ostringstream sstr;
    PrintTPCInfo(sstr, indent, verbosity);
    return sstr.str();
  }

  //......................................................................
  void TPCGeo::InitTPCBoundaries()
  {
    // note that this assumes no rotations of the TPC (except for rotations of a flat
    // angle around one of the three main axes); to avoid this, we should transform the
    // six vertices rather than just the centre

    // we rely on the asumption that the center of TPC is at the local origin
    SetBoundaries(toWorldCoords(LocalPoint_t(-HalfWidth(), -HalfHeight(), -HalfLength())),
                  toWorldCoords(LocalPoint_t(+HalfWidth(), +HalfHeight(), +HalfLength())));

    // the center of the active volume may be elsewhere than the local origin:
    auto const& activeCenter = GetActiveVolumeCenter();
    fActiveBox.SetBoundaries(activeCenter.X() - ActiveHalfWidth(),
                             activeCenter.X() + ActiveHalfWidth(),
                             activeCenter.Y() - ActiveHalfHeight(),
                             activeCenter.Y() + ActiveHalfHeight(),
                             activeCenter.Z() - ActiveHalfLength(),
                             activeCenter.Z() + ActiveHalfLength());
  }

  //......................................................................
  geo::Point_t TPCGeo::GetFrontFaceCenter() const
  {
    auto const& activeBox = ActiveBoundingBox();
    return {activeBox.CenterX(), activeBox.CenterY(), activeBox.MinZ()};
  }

  //......................................................................
  void TPCGeo::UpdateAfterSorting(TPCID tpcid) { fID = tpcid; }

  //......................................................................
  TGeoNode* TPCGeo::NodeForActiveVolume(TGeoNode const* tpc)
  {
    for (int i = 0, nd = tpc->GetNdaughters(); i < nd; ++i) {
      auto daughter = tpc->GetDaughter(i);
      if (daughter && strncmp(daughter->GetName(), "volTPCActive", 12) == 0) return daughter;
    }
    return nullptr;
  }
}
