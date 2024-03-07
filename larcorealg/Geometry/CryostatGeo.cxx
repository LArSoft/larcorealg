////////////////////////////////////////////////////////////////////////
/// @file larcorealg/Geometry/CryostatGeo.cxx
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/CryostatGeo.h"

// LArSoft includes
#include "larcorealg/Geometry/GeoObjectSorter.h" // for GeoObjectSo...

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TClass.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TGeoShape.h" // for TGeoShape

// C++ standard libraries
#include <algorithm>
#include <limits>  // std::numeric_limits<>
#include <sstream> // std::ostringstream
#include <utility> // std::move()
#include <vector>

namespace geo {

  //......................................................................
  CryostatGeo::CryostatGeo(TGeoNode const* node,
                           TransformationMatrix&& trans,
                           TPCList_t&& TPCs,
                           OpDetList_t&& OpDets)
    : fTrans(std::move(trans))
    , fTPCs(std::move(TPCs))
    , fOpDets(std::move(OpDets))
    , fVolume(node->GetVolume())
  {
    // all planes are going to be contained in the volume named volCryostat
    if (!fVolume) throw cet::exception("CryostatGeo") << "cannot find cryostat outline volume\n";

    MF_LOG_DEBUG("Geometry") << "cryostat  volume is " << fVolume->GetName();

    // set the bounding box
    InitCryoBoundaries();
  }

  //......................................................................
  // sort the TPCGeo objects, and the PlaneGeo objects inside
  void CryostatGeo::SortSubVolumes(GeoObjectSorter const& sorter)
  {
    sorter.sort(fTPCs);
    sorter.sort(fOpDets);
  }

  //......................................................................
  void CryostatGeo::UpdateAfterSorting(CryostatID cryoid)
  {
    // update the cryostat ID
    fID = cryoid;

    // trigger all the optical detectors to update
    for (unsigned int opdet = 0; opdet < NOpDet(); ++opdet)
      fOpDets[opdet].UpdateAfterSorting(OpDetID(fID, opdet));

    // trigger all the TPCs to update as well
    for (unsigned int tpc = 0; tpc < NTPC(); ++tpc)
      fTPCs[tpc].UpdateAfterSorting(TPCID(fID, tpc));
  }

  //......................................................................
  TPCGeo const& CryostatGeo::TPC(unsigned int itpc) const
  {
    TPCGeo const* pTPC = TPCPtr(itpc);
    if (!pTPC) {
      throw cet::exception("TPCOutOfRange") << "Request for non-existant TPC " << itpc << "\n";
    }
    return *pTPC;
  }

  //......................................................................
  OpDetGeo const& CryostatGeo::OpDet(unsigned int iopdet) const
  {
    if (iopdet >= fOpDets.size()) {
      throw cet::exception("OpDetOutOfRange") << "Request for non-existant OpDet " << iopdet;
    }
    return fOpDets[iopdet];
  }

  //......................................................................
  auto CryostatGeo::IterateElements() const -> ElementIteratorBox
  {
    return fTPCs;
  }

  //......................................................................
  // wiggle is 1+a small number to allow for rounding errors on the
  // passed in world loc relative to the boundaries.
  TPCID CryostatGeo::PositionToTPCID(Point_t const& point, double wiggle) const
  {
    TPCGeo const* tpc = PositionToTPCptr(point, wiggle);
    return tpc ? tpc->ID() : TPCID{};
  }

  //......................................................................
  // wiggle is 1+a small number to allow for rounding errors on the
  // passed in world loc relative to the boundaries.
  TPCGeo const& CryostatGeo::PositionToTPC(Point_t const& point, double wiggle) const
  {
    TPCGeo const* tpc = PositionToTPCptr(point, wiggle);
    if (!tpc) {
      throw cet::exception("CryostatGeo")
        << "Can't find any TPC for position " << point << " within " << ID() << "\n";
    }
    return *tpc;
  }

  //......................................................................
  TPCGeo const* CryostatGeo::PositionToTPCptr(Point_t const& point, double wiggle) const
  {
    for (auto const& tpc : IterateTPCs())
      if (tpc.ContainsPosition(point, wiggle)) return &tpc;
    return nullptr;
  }

  //......................................................................
  double CryostatGeo::HalfWidth() const
  {
    return static_cast<TGeoBBox const*>(fVolume->GetShape())->GetDX();
  }

  //......................................................................
  double CryostatGeo::HalfHeight() const
  {
    return static_cast<TGeoBBox const*>(fVolume->GetShape())->GetDY();
  }

  //......................................................................
  double CryostatGeo::HalfLength() const
  {
    return static_cast<TGeoBBox const*>(fVolume->GetShape())->GetDZ();
  }

  //......................................................................
  void CryostatGeo::Boundaries(double* boundaries) const
  {
    boundaries[0] = MinX();
    boundaries[1] = MaxX();
    boundaries[2] = MinY();
    boundaries[3] = MaxY();
    boundaries[4] = MinZ();
    boundaries[5] = MaxZ();
  }

  //......................................................................
  std::string CryostatGeo::CryostatInfo(std::string indent /* = "" */,
                                        unsigned int verbosity /* = 1 */) const
  {
    std::ostringstream sstr;
    PrintCryostatInfo(sstr, indent, verbosity);
    return sstr.str();
  }

  //......................................................................
  // Find the nearest opdet to point in this cryostat

  OpDetGeo const* CryostatGeo::GetClosestOpDetPtr(Point_t const& point) const
  {
    unsigned int iOpDet = GetClosestOpDet(point);
    return (iOpDet == std::numeric_limits<double>::max()) ? nullptr : &OpDet(iOpDet);
  }

  //......................................................................
  unsigned int CryostatGeo::GetClosestOpDet(Point_t const& point) const
  {
    unsigned int ClosestDet = std::numeric_limits<unsigned int>::max();
    double ClosestDist = std::numeric_limits<double>::max();

    for (unsigned int o = 0U; o < NOpDet(); ++o) {
      double const ThisDist = OpDet(o).DistanceToPoint(point);
      if (ThisDist < ClosestDist) {
        ClosestDist = ThisDist;
        ClosestDet = o;
      }
    } // for
    return ClosestDet;
  }

  //......................................................................
  unsigned int CryostatGeo::GetClosestOpDet(double const* point) const
  {
    return GetClosestOpDet(vect::makePointFromCoords(point));
  }

  //......................................................................
  void CryostatGeo::InitCryoBoundaries()
  {
    // check that this is indeed a box
    if (!dynamic_cast<TGeoBBox*>(Volume()->GetShape())) {
      // at initialisation time we don't know yet our real ID
      throw cet::exception("CryostatGeo")
        << "Cryostat is not a box! (it is a " << Volume()->GetShape()->IsA()->GetName() << ")\n";
    }

    // get the half width, height, etc of the cryostat
    double const halflength = HalfLength();
    double const halfwidth = HalfWidth();
    double const halfheight = HalfHeight();

    SetBoundaries(toWorldCoords(LocalPoint_t{-halfwidth, -halfheight, -halflength}),
                  toWorldCoords(LocalPoint_t{+halfwidth, +halfheight, +halflength}));
  }

  //......................................................................

} // namespace geo
////////////////////////////////////////////////////////////////////////
