////////////////////////////////////////////////////////////////////////
/// \file larcorealg/Geometry/CryostatGeo.cxx
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////


// class header
#include "larcorealg/Geometry/CryostatGeo.h"

// LArSoft includes
#include "larcorealg/Geometry/GeoObjectSorter.h"          // for GeoObjectSo...

#include "larcorealg/CoreUtils/counter.h"

// Framework includes
#include "cetlib/pow.h" // cet::square()
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TGeoShape.h"                                    // for TGeoShape
#include "TClass.h"

// C++ standard libraries
#include <sstream> // std::ostringstream
#include <limits> // std::numeric_limits<>
#include <vector>
#include <utility> // std::move()


namespace {
  
  //----------------------------------------------------------------------------
  double computeMinDistSq
    (geo::Point_t const& point, geo::BoxBoundedGeo const& box)
  {
    // (squared) distance of the point from the closest face of the box
    using cet::square;
    return std::min({
      square(point.X() - box.MinX()), square(point.X() - box.MaxX()),
      square(point.Y() - box.MinY()), square(point.Y() - box.MaxY()),
      square(point.Z() - box.MinZ()), square(point.Z() - box.MaxZ())
      });
  } // computeMinDistSq()
  
  //----------------------------------------------------------------------------
  
} // local namespace

namespace geo{

  //......................................................................
  CryostatGeo::CryostatGeo(
    TGeoNode const& node, geo::TransformationMatrix&& trans,
    TPCList_t&& TPCs,
    OpDetList_t&& OpDets
    )
    : fTrans(std::move(trans))
    , fTPCs(std::move(TPCs))
    , fOpDets(std::move(OpDets))
    , fVolume(nullptr)
  {

    // all planes are going to be contained in the volume named volCryostat
    // now get the total volume of the Cryostat
    fVolume = node.GetVolume();
    if(!fVolume)
      throw cet::exception("CryostatGeo") << "cannot find cryostat outline volume\n";

    MF_LOG_DEBUG("Geometry") << "cryostat  volume is " << fVolume->GetName();

    // set the bounding box
    InitCryoBoundaries();

    // Set OpDetName;
    fOpDetGeoName = "volOpDetSensitive";

  }


  //......................................................................
  // sort the TPCGeo objects, and the PlaneGeo objects inside
  void CryostatGeo::SortSubVolumes(geo::GeoObjectSorter const& sorter)
  {
    sorter.SortTPCs(fTPCs);

    for (geo::TPCGeo& TPC: fTPCs) {
      TPC.SortSubVolumes(sorter);
    }

    sorter.SortOpDets(fOpDets);
  } // CryostatGeo::SortSubVolumes()


  //......................................................................
  void CryostatGeo::UpdateAfterSorting(geo::CryostatID cryoid) {

    // update the cryostat ID
    fID = cryoid;

    // trigger all the TPCs to update as well
    for (unsigned int tpc = 0; tpc < NTPC(); ++tpc)
      fTPCs[tpc].UpdateAfterSorting(geo::TPCID(fID, tpc));

    // trigger update of all the optical detectors, providing them a reference
    for (auto const iOpDet: util::counter(fOpDets.size())) {
      geo::OpDetGeo& opDet = fOpDets[iOpDet];
      
      std::optional<RefFrame_t const> const refFrame
        { referenceFromClosestTPC(opDet.GetShapeCenter(), fTPCs) };
      
      opDet.UpdateAfterSorting
        (geo::OpDetID(fID, iOpDet), refFrame? &*refFrame: nullptr);
    }

  } // CryostatGeo::UpdateAfterSorting()


  //......................................................................
  const TPCGeo& CryostatGeo::TPC(unsigned int itpc) const
  {
    TPCGeo const* pTPC = TPCPtr(itpc);
    if(!pTPC){
      throw cet::exception("TPCOutOfRange") << "Request for non-existant TPC "
                                            << itpc << "\n";
    }

    return *pTPC;
  }



  //......................................................................
  const OpDetGeo& CryostatGeo::OpDet(unsigned int iopdet) const
  {
    if(iopdet >= fOpDets.size()){
      throw cet::exception("OpDetOutOfRange") << "Request for non-existant OpDet "
                                              << iopdet;
    }

    return fOpDets[iopdet];
  }


  //......................................................................
  auto CryostatGeo::IterateElements() const -> ElementIteratorBox
    { return fTPCs; }


  //......................................................................
  // wiggle is 1+a small number to allow for rounding errors on the
  // passed in world loc relative to the boundaries.
  geo::TPCID::TPCID_t CryostatGeo::FindTPCAtPosition
    (double const worldLoc[3], double wiggle) const
  {
    geo::TPCID tpcid
      = PositionToTPCID(geo::vect::makePointFromCoords(worldLoc), wiggle);
    return tpcid? tpcid.TPC: geo::TPCID::InvalidID;
  } // CryostatGeo::FindTPCAtPosition()

  //......................................................................
  // wiggle is 1+a small number to allow for rounding errors on the
  // passed in world loc relative to the boundaries.
  geo::TPCID CryostatGeo::PositionToTPCID
    (geo::Point_t const& point, double wiggle) const
  {
    geo::TPCGeo const* tpc = PositionToTPCptr(point, wiggle);
    return tpc? tpc->ID(): geo::TPCID{};
  }

  //......................................................................
  // wiggle is 1+a small number to allow for rounding errors on the
  // passed in world loc relative to the boundaries.
  TPCGeo const& CryostatGeo::PositionToTPC
    (geo::Point_t const& point, double wiggle) const
  {
    geo::TPCGeo const* tpc = PositionToTPCptr(point, wiggle);
    if (!tpc) {
      throw cet::exception("CryostatGeo")
        << "Can't find any TPC for position " << point << " within " << ID()
        << "\n";
    }
    return *tpc;
  }

  //......................................................................
  geo::TPCGeo const* CryostatGeo::PositionToTPCptr
    (geo::Point_t const& point, double wiggle) const
  {
    for (auto const& tpc: IterateTPCs())
      if (tpc.ContainsPosition(point, wiggle)) return &tpc;
    return nullptr;
  } // CryostatGeo::PositionToTPCptr()


  //......................................................................
  unsigned int CryostatGeo::MaxPlanes() const {
    unsigned int maxPlanes = 0;
    for (geo::TPCGeo const& TPC: fTPCs) {
      unsigned int maxPlanesInTPC = TPC.Nplanes();
      if (maxPlanesInTPC > maxPlanes) maxPlanes = maxPlanesInTPC;
    } // for
    return maxPlanes;
  } // CryostatGeo::MaxPlanes()

  //......................................................................
  unsigned int CryostatGeo::MaxWires() const {
    unsigned int maxWires = 0;
    for (geo::TPCGeo const& TPC: fTPCs) {
      unsigned int maxWiresInTPC = TPC.MaxWires();
      if (maxWiresInTPC > maxWires) maxWires = maxWiresInTPC;
    } // for
    return maxWires;
  } // CryostatGeo::MaxWires()

  //......................................................................
  double CryostatGeo::HalfWidth()  const
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
  void CryostatGeo::Boundaries(double* boundaries) const {
    boundaries[0] = MinX();
    boundaries[1] = MaxX();
    boundaries[2] = MinY();
    boundaries[3] = MaxY();
    boundaries[4] = MinZ();
    boundaries[5] = MaxZ();
  } // CryostatGeo::CryostatBoundaries(double*)


  //......................................................................
  std::string CryostatGeo::CryostatInfo
    (std::string indent /* = "" */, unsigned int verbosity /* = 1 */) const
  {
    std::ostringstream sstr;
    PrintCryostatInfo(sstr, indent, verbosity);
    return sstr.str();
  } // CryostatGeo::CryostatInfo()

  //......................................................................
  // Find the nearest opdet to point in this cryostat

  geo::OpDetGeo const* CryostatGeo::GetClosestOpDetPtr
    (geo::Point_t const& point) const
  {
    unsigned int iOpDet = GetClosestOpDet(point);
    return
      (iOpDet == std::numeric_limits<double>::max())? nullptr: &OpDet(iOpDet);
  }

  //......................................................................
  unsigned int CryostatGeo::GetClosestOpDet(geo::Point_t const& point) const {
    unsigned int ClosestDet = std::numeric_limits<unsigned int>::max();
    double ClosestDist = std::numeric_limits<double>::max();

    for(unsigned int o = 0U; o < NOpDet(); ++o) {
      double const ThisDist = OpDet(o).DistanceToPoint(point);
      if(ThisDist < ClosestDist) {
        ClosestDist = ThisDist;
        ClosestDet  = o;
      }
    } // for
    return ClosestDet;
  } // CryostatGeo::GetClosestOpDet(geo::Point_t)

  //......................................................................
  unsigned int CryostatGeo::GetClosestOpDet(double const* point) const
    { return GetClosestOpDet(geo::vect::makePointFromCoords(point)); }

  //......................................................................
  void CryostatGeo::InitCryoBoundaries() {

    // check that this is indeed a box
    if (!dynamic_cast<TGeoBBox*>(Volume()->GetShape())) {
      // at initialisation time we don't know yet our real ID
      throw cet::exception("CryostatGeo") << "Cryostat is not a box! (it is a "
        << Volume()->GetShape()->IsA()->GetName() << ")\n";
    }

    // get the half width, height, etc of the cryostat
    const double halflength = HalfLength();
    const double halfwidth  = HalfWidth();
    const double halfheight = HalfHeight();

    SetBoundaries(
      toWorldCoords(LocalPoint_t{ -halfwidth, -halfheight, -halflength }),
      toWorldCoords(LocalPoint_t{ +halfwidth, +halfheight, +halflength })
      );

  } // CryostatGeo::InitCryoBoundaries()


  //......................................................................
  geo::TPCGeo const* CryostatGeo::findClosestTPC
    (geo::Point_t const& point, TPCList_t const& TPCs)
  {
    geo::TPCGeo const* closestTPC = nullptr;
    double minDist2 = std::numeric_limits<double>::max();
    for (geo::TPCGeo const& TPC: TPCs) {
      double const dist2 = computeMinDistSq(point, TPC);
      if (dist2 > minDist2) continue;
      closestTPC = &TPC;
      minDist2 = dist2;
    } // for
    return closestTPC;
  } // CryostatGeo::findClosestTPC()


  //......................................................................
  geo::PlaneGeo const* CryostatGeo::getReferencePlane(geo::TPCGeo const* TPC) {
    return (!TPC || (TPC->Nplanes() == 0))? nullptr: &(TPC->FirstPlane());
  } // CryostatGeo::getReferencePlane()


  //......................................................................
  auto CryostatGeo::referenceFromClosestTPC
    (geo::Point_t const& point, TPCList_t const& TPCs)
    -> std::optional<RefFrame_t>
  {
    geo::TPCGeo const* closestTPC = findClosestTPC(point, TPCs);
    geo::PlaneGeo const* refPlane = getReferencePlane(closestTPC);
    return refPlane
      ? std::make_optional<RefFrame_t>
        (refPlane->WidthDir<geo::Vector_t>(), refPlane->DepthDir<geo::Vector_t>())
      : std::nullopt
      ;
    
  } // CryostatGeo::referenceFromClosestTPC()


  //......................................................................

} // namespace geo
////////////////////////////////////////////////////////////////////////
