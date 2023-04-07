/**
 * @file   larcorealg/Geometry/BoxBoundedGeo.cxx
 * @brief  Provides a additional hit information of a line through the box
 * @author Christoph Rudolf von Rohr (crohr@fnal.gov)
 * @date   July 27th, 2015
 */

#include "larcorealg/Geometry/BoxBoundedGeo.h"

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h"

// ROOT libraries
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"

#include <array>
#include <type_traits>
#include <vector>

namespace geo {
  //----------------------------------------------------------------------------
  bool BoxBoundedGeo::ContainsPosition(TVector3 const& point, double wiggle /* = 1.0 */) const
  {
    return ContainsPosition(vect::toPoint(point), wiggle);
  }

  //----------------------------------------------------------------------------
  bool BoxBoundedGeo::ContainsPosition(double const* point, double wiggle /* = 1.0 */) const
  {
    return ContainsPosition(vect::makePointFromCoords(point), wiggle);
  }

  //----------------------------------------------------------------------------
  std::vector<Point_t> BoxBoundedGeo::GetIntersections(Point_t const& TrajectoryStart,
                                                       Vector_t const& TrajectoryDirect) const
  {
    std::vector<Point_t> IntersectionPoints;
    std::vector<double> LineParameters;

    // Generate normal vectors and offsets for every plane of the box
    // All normal vectors are headed outwards
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    static std::array<Vector_t, 6U> const NormalVectors = {{
      -Xaxis<Vector_t>(),
      Xaxis<Vector_t>(), // anode, cathode,
      -Yaxis<Vector_t>(),
      Yaxis<Vector_t>(), // bottom, top,
      -Zaxis<Vector_t>(),
      Zaxis<Vector_t>() // upstream, downstream
    }};
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    std::array<Point_t, 6U> const NormalVectorOffsets = {{
      Point_t{Min().X(), Min().Y(), Min().Z()}, // Anode offset
      Point_t{Max().X(), Min().Y(), Min().Z()}, // Cathode
      Point_t{Min().X(), Min().Y(), Min().Z()}, // Bottom
      Point_t{Min().X(), Max().Y(), Min().Z()}, // Top
      Point_t{Min().X(), Min().Y(), Min().Z()}, // upstream
      Point_t{Min().X(), Min().Y(), Max().Z()}  // downstream
    }};

    // Loop over all surfaces of the box
    for (unsigned int face_no = 0; face_no < NormalVectors.size(); face_no++) {
      // Check if trajectory and surface are not parallel
      if (NormalVectors[face_no].Dot(TrajectoryDirect)) {
        // Calculate the line parameter for the intersection points
        LineParameters.push_back(
          NormalVectors[face_no].Dot(NormalVectorOffsets.at(face_no) - TrajectoryStart) /
          NormalVectors[face_no].Dot(TrajectoryDirect));
      }
      else
        continue;

      // Calculate intersection point using the line parameter
      IntersectionPoints.push_back(TrajectoryStart + LineParameters.back() * TrajectoryDirect);

      // Coordinate which should be ignored when checking for limits added by Christoph
      // Rudolf von Rohr 05/21/2016
      unsigned int NoCheckCoord;

      // Calculate NoCheckCoord out of the face_no
      if (face_no % 2) {
        // Convert odd face number to coordinate
        NoCheckCoord = (face_no - 1) / 2;
      }
      else {
        // Convert even face number to coordinate
        NoCheckCoord = face_no / 2;
      }

      // Loop over all three space coordinates
      unsigned int coord = 0;
      for (auto extractCoord : vect::coordReaders<Point_t>()) {
        auto const lastPointCoord = vect::bindCoord(IntersectionPoints.back(), extractCoord);
        auto const minCoord = vect::bindCoord(c_min, extractCoord);
        auto const maxCoord = vect::bindCoord(c_max, extractCoord);

        // Changed by Christoph Rudolf von Rohr 05/21/2016 Then check if point is not
        // within the surface limits at this coordinate, without looking at the plane
        // normal vector coordinate. We can assume, that our algorithm already found this
        // coordinate correctily.  In rare cases, looking for boundaries in this
        // coordinate causes unexpected behavior due to floating point inaccuracies.
        if (coord++ != NoCheckCoord &&
            ((lastPointCoord() > maxCoord()) || (lastPointCoord() < minCoord))) {
          // if off limits, get rid of the useless data and break the coordinate loop
          LineParameters.pop_back();
          IntersectionPoints.pop_back();
          break;
        }
      } // coordinate loop
    }   // Surcaces loop

    // sort points according to their parameter value (first is entry, second is exit)
    if (LineParameters.size() == 2 && LineParameters.front() > LineParameters.back()) {
      std::swap(IntersectionPoints.front(), IntersectionPoints.back());
    }

    return IntersectionPoints;
  } // GetIntersections()

  //----------------------------------------------------------------------------
  std::vector<TVector3> BoxBoundedGeo::GetIntersections(TVector3 const& TrajectoryStart,
                                                        TVector3 const& TrajectoryDirect) const
  {
    std::vector<TVector3> intersections;
    for (auto const& point : GetIntersections(Point_t(TrajectoryStart), Vector_t(TrajectoryDirect)))
      intersections.emplace_back(point.X(), point.Y(), point.Z());
    return intersections;
  }

  //----------------------------------------------------------------------------
  void BoxBoundedGeo::SortCoordinates()
  {
    for (auto coordMan : vect::coordManagers<Point_t>()) {
      auto min = vect::bindCoord(c_min, coordMan);
      auto max = vect::bindCoord(c_max, coordMan);
      if (min() > max()) {
        auto temp = min();
        min = max();
        max = temp;
      }
    } // for
  }   // BoxBoundedGeo::SortCoordinates()

  //----------------------------------------------------------------------------

} // namespace geo
