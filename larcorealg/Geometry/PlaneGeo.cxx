////////////////////////////////////////////////////////////////////////
/// \file larcorealg/Geometry/PlaneGeo.cxx
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/PlaneGeo.h"

// LArSoft includes
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireError
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"                // geo::vect::convertTo()
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// Framework includes
#include "cetlib/pow.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TClass.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TMath.h"

// C/C++ standard library
#include <algorithm>
#include <array>
#include <cassert>
#include <functional>  // std::less<>, std::greater<>, std::transform()
#include <iterator>    // std::back_inserter()
#include <sstream>     // std::ostringstream
#include <type_traits> // std::is_same<>, std::decay_t<>

using namespace geo;

namespace {
  /// Returns the offset to apply to value to move it inside [ -limit, +limit ].
  template <typename T>
  T symmetricCapDelta(T value, T limit)
  {
    return (value < -limit) ? -limit - value : (value > +limit) ? +limit - value : 0.0;
  }

  //......................................................................
  /**
   * @brief Class computing the active area of the plane.
   *
   * The active area is defined in the width/depth space which include approximatively all
   * wires.
   *
   * This algorithm requires the frame reference and the wire pitch to be already defined.
   *
   * That area is tuned so that all its points are closer than half a wire pitch from a
   * wire.
   */

  struct ActiveAreaCalculator {

    ActiveAreaCalculator(PlaneGeo const& plane, double margin)
      : plane(plane), wMargin(margin), dMargin(margin)
    {}

    operator PlaneGeo::Rect() { return recomputeArea(); }

  private:
    using Projection_t = ROOT::Math::PositionVector2D<ROOT::Math::Cartesian2D<double>,
                                                      PlaneGeo::WidthDepthReferenceTag>;
    using Vector_t = PlaneGeo::WidthDepthDisplacement_t;

    static_assert(!std::is_same<Projection_t, PlaneGeo::WidthDepthProjection_t>::value,
                  "Necessary maintenance: remove the now optional conversions");

    static constexpr std::size_t kFirstWireStart = 0;
    static constexpr std::size_t kFirstWireEnd = 1;
    static constexpr std::size_t kLastWireStart = 2;
    static constexpr std::size_t kLastWireEnd = 3;

    PlaneGeo const& plane;     ///< Plane to work on.
    double const wMargin;      ///< Margin subtracted from each side of width.
    double const dMargin;      ///< Margin subtracted from each side of depth.
    PlaneGeo::Rect activeArea; ///< Result.

    /// Cache: wire end projections.
    Projection_t wireEnds[4];

    void initializeWireEnds()
    {
      // Collect the projections of the relevant points.
      //
      // Sorted so that start points have width not larger than end points.
      //
      // PointWidthDepthProjection() erroneously returns a vector rather
      // than a point, so a conversion is required
      auto project = [](auto v) -> Projection_t { return {v.X(), v.Y()}; };

      wireEnds[kFirstWireStart] =
        project(plane.PointWidthDepthProjection(plane.FirstWire().GetStart()));
      wireEnds[kFirstWireEnd] =
        project(plane.PointWidthDepthProjection(plane.FirstWire().GetEnd()));
      if (wireEnds[kFirstWireStart].X() > wireEnds[kFirstWireEnd].X())
        std::swap(wireEnds[kFirstWireStart], wireEnds[kFirstWireEnd]);
      wireEnds[kLastWireStart] =
        project(plane.PointWidthDepthProjection(plane.LastWire().GetStart()));
      wireEnds[kLastWireEnd] = project(plane.PointWidthDepthProjection(plane.LastWire().GetEnd()));
      if (wireEnds[kLastWireStart].X() > wireEnds[kLastWireEnd].X())
        std::swap(wireEnds[kLastWireStart], wireEnds[kLastWireEnd]);
    } // initializeWireEnds()

    void includeAllWireEnds()
    {
      //
      // Find the basic area containing all the coordinates.
      //

      // include all the coordinates of the first and last wire
      for (auto const& aWireEnd : wireEnds) {
        activeArea.width.extendToInclude(aWireEnd.X());
        activeArea.depth.extendToInclude(aWireEnd.Y());
      }

    } // includeAllWireEnds()

    void adjustCorners()
    {
      // Modify the corners so that none is father than half a pitch from all wires.
      //
      // directions in wire/depth plane
      Vector_t const widthDir = {1.0, 0.0};
      Vector_t const depthDir = {0.0, 1.0};
      Vector_t wireCoordDir = plane.VectorWidthDepthProjection(plane.GetIncreasingWireDirection());
      double const hp = plane.WirePitch() / 2.0; // half pitch

      // The plan: identify if wires are across or corner, and then:
      // - across:
      //   - identify which sides
      //   - set the farther end of the wire from the side to be p/2 from its corner

      // - corner:
      //   - identify which corners
      //   - move the corners to p/2 from the wires

      // are the wires crossing side to side, as opposed to cutting corners?

      // these are the angles of the original wire coordinate direction
      double const cosAngleWidth = vect::dot(wireCoordDir, widthDir);
      double const cosAngleDepth = vect::dot(wireCoordDir, depthDir);
      // if the wire coordinate direction is on first or third quadrant:
      bool const bPositiveAngle = none_or_both((wireCoordDir.X() >= 0), (wireCoordDir.Y() >= 0));

      // now we readjust the wire coordinate direction to always point toward positive
      // width; this breaks the relation between wireCoordDir and which is the first/last
      // wire
      if (cosAngleWidth < 0) wireCoordDir = -wireCoordDir;

      // let's study the first wire (ends are sorted by width)
      assert(wireEnds[kFirstWireEnd].X() >= wireEnds[kFirstWireStart].X());
      bool const bAlongWidth // horizontal
        = equal(wireEnds[kFirstWireEnd].X(), activeArea.width.upper) &&
          equal(wireEnds[kFirstWireStart].X(), activeArea.width.lower);
      bool const bAlongDepth =
        !bAlongWidth && // vertical
        equal(std::max(wireEnds[kFirstWireStart].Y(), wireEnds[kFirstWireEnd].Y()),
              activeArea.depth.upper) &&
        equal(std::min(wireEnds[kFirstWireStart].Y(), wireEnds[kFirstWireEnd].Y()),
              activeArea.depth.lower);
      assert(!(bAlongWidth && bAlongDepth));

      if (bAlongWidth) { // horizontal

        // +---------+
        // |   ___,--|  upper width bound
        // |--'      |

        // find which is the wire with higher width:
        // the last wire is highest if the wire coordinate direction (which is defined by
        // what is first and what is last) is parallel to the width direction
        std::size_t const iUpperWire = (cosAngleDepth > 0) ? kLastWireStart : kFirstWireStart;
        // largest distance from upper depth bound of the two ends of wire
        double const maxUpperDistance =
          activeArea.depth.upper -
          std::min(wireEnds[iUpperWire].Y(), wireEnds[iUpperWire ^ 0x1].Y());
        // set the upper side so that the maximum distance is p/2 (it may be actually less
        // if the wire is not perfectly horizontal)
        activeArea.depth.upper += (hp - maxUpperDistance);

        // |--.___   |
        // |      `--|  deal with the lower bound now
        // +---------+

        std::size_t const iLowerWire = (cosAngleDepth > 0) ? kFirstWireStart : kLastWireStart;
        // largest distance from lower depth bound of the two ends of wire
        double const maxLowerDistance =
          std::max(wireEnds[iLowerWire].Y(), wireEnds[iLowerWire ^ 0x1].Y()) -
          activeArea.depth.lower;
        // set the upper side so that the minimum distance is p/2
        activeArea.depth.lower -= (hp - maxLowerDistance);
      }                       // horizontal wires
      else if (bAlongDepth) { // vertical
        // --,---+
        //   |   |
        //    \  |
        //    |  |  upper depth bound
        //     \ |
        //     | |
        // ------+

        // find which is the wire with higher depth:
        // the last wire is highest if the wire coordinate direction (which is defined by
        // what is first and what is last) is parallel to the depth direction
        std::size_t const iUpperWire = (cosAngleWidth > 0) ? kLastWireStart : kFirstWireStart;
        // largest distance from upper depth bound of the two ends of wire
        double const maxUpperDistance =
          activeArea.width.upper -
          std::min(wireEnds[iUpperWire].X(), wireEnds[iUpperWire ^ 0x1].X());
        // set the upper side so that the minimum distance is p/2
        activeArea.width.upper += (hp - maxUpperDistance);

        // +-,----
        // | |
        // |  \     .
        // |  |     deal with the lower bound now
        // |   \    .
        // |   |
        // +------
        std::size_t const iLowerWire = (cosAngleWidth > 0) ? kFirstWireStart : kLastWireStart;
        // largest distance from lower width bound of the two ends of wire
        double const maxLowerDistance =
          std::max(wireEnds[iLowerWire].X(), wireEnds[iLowerWire ^ 0x1].X()) -
          activeArea.width.lower;
        // set the upper side so that the minimum distance is p/2
        activeArea.width.lower -= (hp - maxLowerDistance);

      }                          // vertical wires
      else if (bPositiveAngle) { // wires are not going across: corners!
        // corners at (lower width, lower depth), (upper width, upper depth)

        // -._------+
        //    `-._  |  upper width corner (upper depth)
        //        `-|

        // start of the wire on the upper corner (width coordinate is lower for start than
        // for end)
        std::size_t const iUpperWire = (cosAngleWidth > 0) ? kLastWireStart : kFirstWireStart;

        double const upperDistance =
          vect::dot(Vector_t(activeArea.width.upper - wireEnds[iUpperWire].X(), 0.0), wireCoordDir);
        // make the upper distance become p/2
        auto const upperDelta = (hp - upperDistance) * wireCoordDir;
        activeArea.width.upper += upperDelta.X();
        activeArea.depth.upper += upperDelta.Y();

        // |-._
        // |   `-._    lower width corner (lower depth)
        // +-------`-

        // end of the wire on the lower corner (width coordinate is lower than the end)
        std::size_t const iLowerWire = (cosAngleWidth > 0) ? kFirstWireEnd : kLastWireEnd;

        double const lowerDistance =
          vect::dot(Vector_t(wireEnds[iLowerWire].X() - activeArea.width.lower, 0.0), wireCoordDir);
        // make the lower distance become p/2 (note direction of wire coord)
        auto const lowerDelta = (hp - lowerDistance) * wireCoordDir;
        activeArea.width.lower -= lowerDelta.X();
        activeArea.depth.lower -= lowerDelta.Y();
      }
      else { // !bPositiveAngle
        // corners at (lower width, upper depth), (upper width, lower depth)

        //       _,-|
        //   _,-'   |  upper width corner (lower depth)
        // -'-------+

        // start of the wire on the upper corner (width coordinate is lower than the end)
        std::size_t const iUpperWire = (cosAngleWidth > 0) ? kLastWireStart : kFirstWireStart;

        double const upperDistance =
          vect::dot(Vector_t(activeArea.width.upper - wireEnds[iUpperWire].X(), 0.0), wireCoordDir);
        // make the upper distance become p/2
        auto const upperDelta = (hp - upperDistance) * wireCoordDir;
        activeArea.width.upper += upperDelta.X();
        activeArea.depth.lower += upperDelta.Y();

        // +------_,-
        // |  _,-'     lower width corner (upper depth)
        // |-'

        // end of the wire on the lower corner (width coordinate is lower than the end)
        std::size_t const iLowerWire = (cosAngleWidth > 0) ? kFirstWireEnd : kLastWireEnd;

        double const lowerDistance =
          vect::dot(Vector_t(wireEnds[iLowerWire].X() - activeArea.width.lower, 0.0), wireCoordDir);
        // make the lower distance become p/2 (note direction of wire coord)
        auto const lowerDelta = (hp - lowerDistance) * wireCoordDir;
        activeArea.width.lower -= lowerDelta.X();
        activeArea.depth.upper -= lowerDelta.Y();

      } // if ...

    } // adjustCorners()

    void applyMargin()
    {
      if (wMargin != 0.0) {
        activeArea.width.lower += wMargin;
        activeArea.width.upper -= wMargin;
      }
      if (dMargin != 0.0) {
        activeArea.depth.lower += dMargin;
        activeArea.depth.upper -= dMargin;
      }
    } // applyMargin()

    PlaneGeo::Rect recomputeArea()
    {
      activeArea = {};

      // 0. collect the projections of the relevant point
      initializeWireEnds();

      // 1. find the basic area containing all the coordinates
      includeAllWireEnds();

      // 2. adjust area so that no corner is father than half a wire pitch
      adjustCorners();

      // 3. apply an absolute margin
      applyMargin();

      return activeArea;
    } // computeArea()

    /// Returns true if a and b are both true or both false (exclusive nor).
    static bool none_or_both(bool a, bool b) { return a == b; }

    /// Returns whether the two numbers are the same, lest a tolerance.
    template <typename T>
    static bool equal(T a, T b, T tol = T(1e-5))
    {
      return std::abs(a - b) <= tol;
    }

  }; // struct ActiveAreaCalculator
}

namespace geo {

  //----------------------------------------------------------------------------
  //---  geo::PlaneGeo
  //---
  PlaneGeo::PlaneGeo(TGeoNode const* node, TransformationMatrix&& trans, WireCollection_t&& wires)
    : fNode(node), fTrans(std::move(trans)), fWire(std::move(wires))
  {
    if (!fNode->GetVolume()) {
      throw cet::exception("PlaneGeo")
        << "Plane geometry node " << fNode->IsA()->GetName() << "[" << fNode->GetName() << ", #"
        << fNode->GetNumber() << "] has no volume!\n";
    }

    // view is now set at TPC level with SetView

    DetectGeometryDirections();
    UpdateWirePitchSlow();
  }

  //......................................................................

  BoxBoundedGeo PlaneGeo::BoundingBox() const
  {
    //
    // The algorithm is not very refined...
    //

    TGeoBBox const* pShape = dynamic_cast<TGeoBBox const*>(Volume().GetShape());
    if (!pShape) {
      throw cet::exception("PlaneGeo")
        << "BoundingBox(): volume " << Volume().IsA()->GetName() << " is not a TGeoBBox!";
    }

    BoxBoundedGeo box;
    unsigned int points = 0;
    for (double dx : {-(pShape->GetDX()), +(pShape->GetDX())}) {
      for (double dy : {-(pShape->GetDY()), +(pShape->GetDY())}) {
        for (double dz : {-(pShape->GetDZ()), +(pShape->GetDZ())}) {

          auto const p = toWorldCoords(LocalPoint_t{dx, dy, dz});

          if (points++ == 0)
            box.SetBoundaries(p, p);
          else
            box.ExtendToInclude(p);

        } // for z
      }   // for y
    }     // for x
    return box;
  }

  //......................................................................

  WireGeo const& PlaneGeo::Wire(unsigned int iwire) const
  {
    if (WireGeo const* pWire = WirePtr(iwire)) { return *pWire; }
    throw cet::exception("WireOutOfRange") << "Request for non-existant wire " << iwire << "\n";
  }

  //......................................................................
  void PlaneGeo::SortWires(Compare<WireGeo> cmp) { std::sort(fWire.begin(), fWire.end(), cmp); }

  //......................................................................
  bool PlaneGeo::WireIDincreasesWithZ() const
  {
    return lar::util::RealComparisons(1e-3).nonNegative(GetIncreasingWireDirection().Z());
  }

  //......................................................................
  std::string PlaneGeo::PlaneInfo(std::string indent /* = "" */,
                                  unsigned int verbosity /* = 1 */) const
  {
    std::ostringstream sstr;
    PrintPlaneInfo(sstr, indent, verbosity);
    return sstr.str();
  }

  //......................................................................
  PlaneGeo::WidthDepthProjection_t PlaneGeo::DeltaFromPlane(WidthDepthProjection_t const& proj,
                                                            double wMargin,
                                                            double dMargin) const
  {
    return {symmetricCapDelta(proj.X(), fFrameSize.HalfWidth() - wMargin),
            symmetricCapDelta(proj.Y(), fFrameSize.HalfDepth() - dMargin)};
  }

  //......................................................................
  PlaneGeo::WidthDepthProjection_t PlaneGeo::DeltaFromActivePlane(
    WidthDepthProjection_t const& proj,
    double wMargin,
    double dMargin) const
  {
    return {fActiveArea.width.delta(proj.X(), wMargin), fActiveArea.depth.delta(proj.Y(), dMargin)};
  }

  //......................................................................
  bool PlaneGeo::isProjectionOnPlane(Point_t const& point) const
  {
    auto const deltaProj = DeltaFromPlane(PointWidthDepthProjection(point));
    return (deltaProj.X() == 0.) && (deltaProj.Y() == 0.);
  }

  //......................................................................
  PlaneGeo::WidthDepthProjection_t PlaneGeo::MoveProjectionToPlane(
    WidthDepthProjection_t const& proj) const
  {
    // We have a more complicated implementation to avoid rounding errors.  In this way,
    // the result is really guaranteed to be exactly on the border.
    auto const delta = DeltaFromPlane(proj);
    return {(delta.X() == 0.0) ?
              proj.X() :
              ((delta.X() > 0) ? -fFrameSize.HalfWidth() // delta positive -> proj on negative side
                                 :
                                 fFrameSize.HalfWidth()),
            (delta.Y() == 0.0) ?
              proj.Y() :
              ((delta.Y() > 0) ? -fFrameSize.HalfDepth() // delta positive -> proj on negative side
                                 :
                                 fFrameSize.HalfDepth())};
  }

  //......................................................................
  Point_t PlaneGeo::MovePointOverPlane(Point_t const& point) const
  {
    // This implementation is subject to rounding errors, since the result of the addition
    // might jitter above or below the border.

    auto const deltaProj = DeltaFromPlane(PointWidthDepthProjection(point));
    return point + deltaProj.X() * WidthDir() + deltaProj.Y() * DepthDir();
  }

  //......................................................................
  WireID PlaneGeo::NearestWireID(Point_t const& pos) const
  {
    // 1) compute the wire coordinate of the point
    // 2) get the closest wire number
    // 3) check if the wire does exist
    // 4) build and return the wire ID

    // this line merges parts (1) and (2); add 0.5 to have the correct rounding:
    int nearestWireNo = int(0.5 + WireCoordinate(pos));

    // if we are outside of the wireplane range, throw an exception
    if ((nearestWireNo < 0) || ((unsigned int)nearestWireNo >= Nwires())) {

      auto wireNo = nearestWireNo; // save for the output

      if (nearestWireNo < 0)
        wireNo = 0;
      else
        wireNo = Nwires() - 1;

      throw InvalidWireError("Geometry", ID(), nearestWireNo, wireNo)
        << "Can't find nearest wire for position " << pos << " in plane " << std::string(ID())
        << " approx wire number # " << wireNo << " (capped from " << nearestWireNo << ")\n";
    } // if invalid

    return {ID(), (WireID::WireID_t)nearestWireNo};
  }

  //......................................................................
  WireGeo const& PlaneGeo::NearestWire(Point_t const& point) const
  {
    // Note that this code is ready for when NearestWireID() will be changed to return an
    // invalid ID instead of throwing.  As things are now, `NearestWireID()` will never
    // return an invalid ID, but it will throw an exception similar to this one.

    WireID const wireID = NearestWireID(point);
    if (wireID) return Wire(wireID); // we have that wire, so we return it

    // wire ID is invalid, meaning it's out of range. Throw an exception!
    WireID const closestID = ClosestWireID(wireID);
    throw InvalidWireError("Geometry", ID(), closestID.Wire, wireID.Wire)
      << "Can't find nearest wire for position " << point << " in plane " << std::string(ID())
      << " approx wire number # " << closestID.Wire << " (capped from " << wireID.Wire << ")\n";
  }

  //......................................................................
  double PlaneGeo::InterWireProjectedDistance(WireCoordProjection_t const& projDir) const
  {
    assert(lar::util::Vector2DComparison{1e-6}.nonZero(projDir));
    return std::sqrt(cet::square(projDir.X() / projDir.Y()) + 1.0) * fWirePitch;
  }

  //......................................................................
  double PlaneGeo::InterWireDistance(Vector_t const& dir) const
  {
    // the secondary component of the wire decomposition basis is wire coord.
    double const r = dir.R();
    assert(r >= 1.e-6);

    double const absWireCoordProj = std::abs(fDecompWire.VectorSecondaryComponent(dir));
    return r / absWireCoordProj * fWirePitch;
  }

  //......................................................................
  double PlaneGeo::ThetaZ() const { return FirstWire().ThetaZ(); }

  //......................................................................
  void PlaneGeo::UpdateAfterSorting(PlaneID planeid, BoxBoundedGeo const& TPCbox)
  {
    // the order here matters

    // reset our ID
    fID = planeid;

    UpdatePlaneNormal(TPCbox);
    UpdateWidthDepthDir();
    UpdateIncreasingWireDir();

    // update wires
    WireID::WireID_t wireNo = 0;
    for (auto& wire : fWire) {
      wire.UpdateAfterSorting(WireID(fID, wireNo), shouldFlipWire(wire));
      ++wireNo;
    } // for wires

    UpdateDecompWireOrigin();
    UpdateWireDir();
    UpdateWirePlaneCenter();
    UpdateOrientation();
    UpdateWirePitch();
    UpdateActiveArea();
    UpdatePhiZ();
    UpdateView();
  }

  //......................................................................
  std::string PlaneGeo::ViewName(View_t view)
  {
    switch (view) {
    case kU: return "U";
    case kV: return "V";
    case kZ: return "Z";
    case kY: return "Y";
    case kX: return "X";
    case k3D: return "3D";
    case kUnknown: return "?";
    default: return "<UNSUPPORTED (" + std::to_string((int)view) + ")>";
    }
  }

  //......................................................................
  std::string PlaneGeo::OrientationName(Orient_t orientation)
  {
    switch (orientation) {
    case kHorizontal: return "horizontal"; break;
    case kVertical: return "vertical"; break;
    default: return "unexpected"; break;
    }
  }

  //......................................................................
  void PlaneGeo::DetectGeometryDirections()
  {
    // We need to identify which are the "long" directions of the plane.  We assume it is
    // a box, and the shortest side is excluded.  The first direction ("width") is given
    // by preference to z.  If z is the direction of the normal to the plane... oh well.
    // Let's say privilege to the one which comes from local z, then y.  That means:
    // undefined.
    //
    // Requirements:
    //  - ROOT geometry information (shapes and transformations)
    //  - the shape must be a box (an error is PRINTED if not)
    //  - center of the wire plane (not just the center of the plane box)

    // how do they look like in the world?
    TGeoBBox const* pShape = dynamic_cast<TGeoBBox const*>(Volume().GetShape());
    if (!pShape) {
      mf::LogError("BoxInfo") << "Volume " << Volume().IsA()->GetName()
                              << " is not a TGeoBBox! Dimensions won't be available.";
      // set it invalid
      fDecompFrame.SetOrigin(origin());
      fDecompFrame.SetMainDir({0., 0., 0.});
      fDecompFrame.SetSecondaryDir({0., 0., 0.});
      fFrameSize = {0.0, 0.0};
      return;
    }

    std::array<Vector_t, 3U> sides;
    size_t iSmallest = 3;
    {
      size_t iSide = 0;

      sides[iSide] = toWorldCoords(LocalVector_t{pShape->GetDX(), 0.0, 0.0});
      iSmallest = iSide;
      ++iSide;

      sides[iSide] = toWorldCoords(LocalVector_t{0.0, pShape->GetDY(), 0.0});
      if (sides[iSide].Mag2() < sides[iSmallest].Mag2()) iSmallest = iSide;
      ++iSide;

      sides[iSide] = toWorldCoords(LocalVector_t{0.0, 0.0, pShape->GetDZ()});
      if (sides[iSide].Mag2() < sides[iSmallest].Mag2()) iSmallest = iSide;
      ++iSide;
    }

    // which are the largest ones?
    size_t kept[2];
    {
      size_t iKept = 0;
      for (size_t i = 0; i < 3; ++i)
        if (i != iSmallest) kept[iKept++] = i;
    }

    // which is which?
    //
    // Pick width as the most z-like.
    size_t const iiWidth =
      std::abs(sides[kept[0]].Unit().Z()) > std::abs(sides[kept[1]].Unit().Z()) ? 0 : 1;
    size_t const iWidth = kept[iiWidth];
    size_t const iDepth = kept[1 - iiWidth]; // the other

    fDecompFrame.SetMainDir(vect::rounded01(sides[iWidth].Unit(), 1e-4));
    fDecompFrame.SetSecondaryDir(vect::rounded01(sides[iDepth].Unit(), 1e-4));
    fFrameSize.halfWidth = sides[iWidth].R();
    fFrameSize.halfDepth = sides[iDepth].R();
  }

  //......................................................................
  Vector_t PlaneGeo::GetNormalAxis() const
  {
    unsigned const int NWires = Nwires();
    if (NWires < 2) return {}; // why are we even here?

    // 1) get the direction of the middle wire
    auto const WireDir = Wire(NWires / 2).Direction();

    // 2) get the direction between the middle wire and the next one
    auto const ToNextWire = Wire(NWires / 2 + 1).GetCenter() - Wire(NWires / 2).GetCenter();

    // 3) get the direction perpendicular to the plane
    // 4) round it
    // 5) return its norm
    return vect::rounded01(WireDir.Cross(ToNextWire).Unit(), 1e-4);
  }

  //......................................................................
  void PlaneGeo::UpdateOrientation()
  {
    // this algorithm needs to know about the axis; the normal is expected to be already
    // updated.

    // sanity check
    if (fWire.size() < 2) {
      // this likely means construction is not complete yet
      throw cet::exception("NoWireInPlane")
        << "PlaneGeo::UpdateOrientation(): only " << fWire.size() << " wires!\n";
    } // if

    auto normal = GetNormalDirection();

    if (std::abs(std::abs(normal.X()) - 1.) < 1e-3)
      fOrientation = kVertical;
    else if (std::abs(std::abs(normal.Y()) - 1.) < 1e-3)
      fOrientation = kHorizontal;
    else {
      // at this point, the only problem is the lack of a label for this orientation;
      // probably introducing a geo::kOtherOrientation would suffice
      throw cet::exception("Geometry")
        << "Plane with unsupported orientation (normal: " << normal << ")\n";
    }
  }

  //......................................................................
  void PlaneGeo::UpdateWirePitch()
  {
    // pick long wires around the center of the detector, so that their coordinates are
    // defined with better precision
    assert(Nwires() > 1);

    auto const iWire = Nwires() / 2;

    fWirePitch = WireGeo::WirePitch(Wire(iWire - 1), Wire(iWire));
  }

  //......................................................................
  void PlaneGeo::UpdatePhiZ()
  {
    auto const& wire_coord_dir = GetIncreasingWireDirection();
    fCosPhiZ = wire_coord_dir.Z();
    fSinPhiZ = wire_coord_dir.Y();
  }

  void PlaneGeo::UpdateView()
  {
    /**
     * This algorithm assigns views according to the angle the wire axis cuts with y axis
     * ("thetaY"), but from the point of view of the center of the TPC.  A special case is
     * when the drift axis is on y axis.
     *
     * In the normal case, the discrimination happens on the the arctangent of the point {
     * (y,w), (y x n,w) }, where w is the wire direction, y is the coordinate axis and n
     * the normal to the wire plane. This definition gives the same value regardless of
     * the direction of w on its axis.
     *
     * If thetaY is 0, wires are parallel to the y axis, the view is assigned as kX or kZ
     * depending on whether the plane normal is closer to the z axis or the x axis,
     * respectively (the normal describes a direction _not_ measured by the wires).
     *
     * If thetaY is a right angle, the wires are orthogonal to y axis and view kY view is
     * assigned.
     * If thetaY is smaller than 0, the view is called "U".
     * If thetaY is larger than 0, the view is called "V".
     *
     * The special case where the drift axis is on y axis is treated separately.  In that
     * case, the role of y axis is replaced by the z axis and the discriminating figure is
     * equivalent to the usual ThetaZ().
     *
     * If thetaZ is 0, the wires are measuring x and kX view is chosen.
     * If thetaZ is a right angle, the wires are measuring z and kZ view is chosen.
     * If thetaZ is smaller than 0, the view is called "U".
     * If thetaZ is larger than 0, the view is called "V".
     */

    auto const& normalDir = GetNormalDirection();
    auto const& wireDir = GetWireDirection();

    // normal direction has been rounded, so exact comparison can work
    if (std::abs(normalDir.Y()) != 1.0) {
      // normal case: drift direction is not along y (vertical)

      // yw is pretty much GetWireDirection().Y()...  thetaY is related to atan2(ynw, yw)
      double const yw = vect::dot(wireDir, Yaxis());
      double const ynw = vect::mixedProduct(Yaxis(), normalDir, wireDir);

      if (std::abs(yw) < 1.0e-4) { // wires orthogonal to y axis
        double const closeToX = std::abs(vect::dot(normalDir, Xaxis()));
        double const closeToZ = std::abs(vect::dot(normalDir, Zaxis()));
        SetView((closeToZ > closeToX) ? kX : kY);
      }
      else if (std::abs(ynw) < 1.0e-4) { // wires parallel to y axis
        SetView(kZ);
      }
      else if ((ynw * yw) < 0)
        SetView(kU); // different sign => thetaY > 0
      else if ((ynw * yw) > 0)
        SetView(kV); // same sign => thetaY < 0
      else
        assert(false); // logic error?!
    }
    else { // if drift is vertical
      // special case: drift direction is along y (vertical)

      // zw is pretty much GetWireDirection().Z()...
      double const zw = vect::dot(wireDir, Zaxis());
      // while GetNormalDirection() axis is on y, its direction is not fixed:
      double const znw = vect::mixedProduct(Zaxis(), normalDir, wireDir);

      // thetaZ is std::atan(znw/zw)

      if (std::abs(zw) < 1.0e-4) { // orthogonal to z, orthogonal to y...
        // this is equivalent to thetaZ = +/- pi/2
        SetView(kZ);
      }
      else if (std::abs(znw) < 1.0e-4) { // parallel to z, orthogonal to y...
        // this is equivalent to thetaZ = 0
        SetView(kX);
      }
      else if ((znw * zw) < 0)
        SetView(kU); // different sign => thetaZ > 0
      else if ((znw * zw) > 0)
        SetView(kV); // same sign => thetaZ < 0
      else
        assert(false); // logic error?!

    } // if drift direction... else
  }

  //......................................................................
  void PlaneGeo::UpdatePlaneNormal(BoxBoundedGeo const& TPCbox)
  {
    // direction normal to the wire plane, points toward the center of TPC

    // start from the axis
    fNormal = GetNormalAxis();

    // now evaluate where we are pointing
    auto const towardCenter = TPCbox.Center() - GetBoxCenter();

    // if they are pointing in opposite directions, flip the normal
    if (fNormal.Dot(towardCenter) < 0) fNormal = -fNormal;
    vect::round01(fNormal, 1e-3);
  }

  //......................................................................
  void PlaneGeo::UpdateWidthDepthDir()
  {
    // fix the positiveness of the width/depth/normal frame

    // The basis is already set and orthonormal, with only the width and depth directions
    // arbitrary.  We choose the direction of the secondary axis ("depth") so that the
    // frame normal is oriented in the general direction of the plane normal (the latter
    // is computed independently).
    if (WidthDir().Cross(DepthDir()).Dot(GetNormalDirection()) < 0.0) {
      fDecompFrame.SetSecondaryDir(vect::rounded01(-fDecompFrame.SecondaryDir(), 1e-4));
    }
  }

  //......................................................................
  void PlaneGeo::UpdateIncreasingWireDir()
  {
    // Direction measured by the wires, pointing toward increasing wire number; requires:
    // - the normal to the plane to be correct
    // - wires to be sorted

    // 1) get the direction of the middle wire
    auto refWireNo = Nwires() / 2;
    if (refWireNo == Nwires() - 1) --refWireNo;
    auto const& refWire = Wire(refWireNo);
    auto const& WireDir = refWire.Direction(); // we only rely on the axis

    // 2) get the axis perpendicular to it on the wire plane (arbitrary direction)
    auto wireCoordDir = GetNormalDirection().Cross(WireDir).Unit();

    // 3) where is the next wire?
    auto toNextWire = Wire(refWireNo + 1).GetCenter() - refWire.GetCenter();

    // 4) if wireCoordDir is pointing away from the next wire, flip it
    if (wireCoordDir.Dot(toNextWire) < 0) { wireCoordDir = -wireCoordDir; }
    fDecompWire.SetSecondaryDir(vect::rounded01(wireCoordDir, 1e-4));
  }

  //......................................................................
  void PlaneGeo::UpdateWireDir()
  {
    fDecompWire.SetMainDir(vect::rounded01(FirstWire().Direction(), 1e-4));

    // check that the resulting normal matches the plane one
    assert(
      lar::util::makeVector3DComparison(1e-5).equal(fDecompWire.NormalDir(), GetNormalDirection()));
  }

  //......................................................................
  void PlaneGeo::UpdateWirePitchSlow()
  {
    if (fWire.empty()) { return; }

    // Compare one wire (the first one, for convenience) with all other wires; the wire
    // pitch is the smallest distance we find.
    //
    // This algorithm assumes wire pitch is constant, but it does not assume wire ordering
    auto firstWire = fWire.cbegin(), wire = firstWire, wend = fWire.cend();
    fWirePitch = WireGeo::WirePitch(*firstWire, *(++wire));

    while (++wire != wend) {
      auto wirePitch = WireGeo::WirePitch(*firstWire, *wire);
      if (wirePitch < 1e-4) continue; // it's 0!
      if (wirePitch < fWirePitch) fWirePitch = wirePitch;
    }
  }

  //......................................................................
  void PlaneGeo::UpdateDecompWireOrigin()
  {
    // update the origin of the reference frame (the middle of the first wire)
    fDecompWire.SetOrigin(vect::toPoint(FirstWire().GetCenter()));
  }

  //......................................................................
  void PlaneGeo::UpdateActiveArea()
  {
    // The active area is defined in the width/depth space which include approximatively
    // all wires.
    //
    // See `ActiveAreaCalculator` for details of the algorithm.

    // we scratch 1 um from each side to avoid rounding errors later
    fActiveArea = ActiveAreaCalculator(*this, 0.0001);
  }

  //......................................................................
  void PlaneGeo::UpdateWirePlaneCenter()
  {
    // The center of the wire plane is defined as the center of the plane box, translated
    // to the plane the wires lie on.  This assumes that the thickness direction of the
    // box is aligned with the drift direction, so that the translated point is still in
    // the middle of width and depth dimensions.  It is possible to remove that assumption
    // by translating the center of the box along the thickness direction enough to bring
    // it to the wire plane.  The math is just a bit less straightforward, so we don't
    // bother yet.
    //
    // Requirements:
    //  * the wire decomposition frame must be set up (at least its origin and normal
    //    direction)

    fCenter = GetBoxCenter();

    DriftPoint(fCenter, DistanceFromPlane(fCenter));

    vect::round0(fCenter, 1e-7); // round dimensions less than 1 nm to 0

    fDecompFrame.SetOrigin(fCenter); // equivalent to GetCenter() now

  } // PlaneGeo::UpdateWirePlaneCenter()

  //......................................................................
  bool PlaneGeo::shouldFlipWire(WireGeo const& wire) const
  {
    // The correct orientation is so that:
    //
    // (direction) x (wire coordinate direction) . (plane normal)
    //
    // is positive; it it's negative, then we should flip the wire.
    //
    // Note that the increasing wire direction comes from the wire frame, while the normal
    // direction is computed independently by geometry.  The resulting normal in the wire
    // frame is expected to be the same as the plane normal from GetNormalDirection(); if
    // this is not the case, flipping the wire direction should restore it.

    return wire.Direction().Cross(GetIncreasingWireDirection()).Dot(GetNormalDirection()) <
           +0.5; // should be in fact exactly +1
  }

} // namespace geo
