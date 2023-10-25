////////////////////////////////////////////////////////////////////////
/// @file  larcorealg/Geometry/WireGeo.cxx
/// @brief Encapsulate the geometry of a wire
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/WireGeo.h"

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h"                // geo::vect
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util ns

// framework
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT
#include "TGeoNode.h"
#include "TGeoTube.h"

// C/C++ libraries
#include <algorithm> // std::clamp()...
#include <cassert>
#include <cmath> // std::cos(), std::isfinite()...
#include <sstream>

namespace geo {

  //-----------------------------------------
  WireGeo::WireGeo(TGeoNode const* node, TransformationMatrix&& trans)
    : fWireNode(node), fTrans(std::move(trans)), flipped(false)
  {
    fHalfL = ((TGeoTube*)fWireNode->GetVolume()->GetShape())->GetDZ();

    // uncomment the following to check the paths to the wires
    //   std::string p(base);
    //   for(int i = 0; i <= depth; ++i){
    //     p += "/";
    //     p += path[i]->GetName();
    //   }
    //   std::cout << p.c_str() << std::endl;

    // determine the orientation of the wire
    auto lp = origin<LocalPoint_t>();
    fCenter = toWorldCoords(lp);

    lp.SetZ(fHalfL);
    auto end = toWorldCoords(lp);

    fThetaZ = std::acos(std::clamp((end.Z() - fCenter.Z()) / fHalfL, -1.0, +1.0));

    // check to see if it runs "forward" or "backwards" in z
    // check is made looking at the y position of the end point
    // relative to the center point because you want to know if
    // the end point is above or below the center of the wire in
    // the yz plane
    if (end.Y() < fCenter.Y()) fThetaZ *= -1.;

    //This ensures we are looking at the angle between 0 and Pi
    //as if the wire runs at one angle it also runs at that angle +-Pi
    if (fThetaZ < 0) fThetaZ += util::pi();

    assert(std::isfinite(fThetaZ));

  } // geo::WireGeo::WireGeo()

  //......................................................................
  double geo::WireGeo::RMax() const
  {
    return ((TGeoTube*)fWireNode->GetVolume()->GetShape())->GetRmax();
  }

  //......................................................................
  double geo::WireGeo::RMin() const
  {
    return ((TGeoTube*)fWireNode->GetVolume()->GetShape())->GetRmin();
  }

  //......................................................................
  double WireGeo::ThetaZ(bool degrees) const
  {
    return degrees ? util::RadiansToDegrees(fThetaZ) : fThetaZ;
  }

  //......................................................................
  std::string WireGeo::WireInfo(std::string indent /* = "" */,
                                unsigned int verbosity /* = 1 */) const
  {
    std::ostringstream sstr;
    PrintWireInfo(sstr, indent, verbosity);
    return sstr.str();
  } // WireGeo::WireInfo()

  //......................................................................
  double WireGeo::DistanceFrom(WireGeo const& wire) const
  {
    //
    // The algorithm assumes that picking any point on the wire will do,
    // that is, that the wires are parallel.
    //

    if (!isParallelTo(wire)) return 0;

    // a vector connecting to the other wire
    auto const toWire = wire.GetCenter() - GetCenter();

    // the distance is that vector, times the sine of the angle with the wire
    // direction; we get that in a generic way with a cross product.
    // We don't even care about the sign here
    // (if we did, we would do a dot-product with the normal to the plane,
    // and we should get a positive distance if the other wire has larger wire
    // coordinate than this one).
    return toWire.Cross(Direction()).R();

  } // WireGeo::DistanceFrom()

  //......................................................................
  void WireGeo::UpdateAfterSorting(WireID const&, bool flip)
  {

    // flip, if asked
    if (flip) Flip();

  } // WireGeo::UpdateAfterSorting()

  //......................................................................
  void WireGeo::Flip()
  {
    // we don't need to do much to flip so far:
    // - ThetaZ() is defined in [0, pi], invariant to flipping
    // - we don't change the transformation matrices, that we want to be
    //   untouched and coherent with the original geometry source
    // - center is invariant for flipping
    // - start and end are computed on the fly (taking flipping into account)
    // - ... and we chose to leave half length unsigned and independent

    // change the flipping bit
    flipped = !flipped;

  } // WireGeo::Flip()

  //......................................................................

}
////////////////////////////////////////////////////////////////////////
