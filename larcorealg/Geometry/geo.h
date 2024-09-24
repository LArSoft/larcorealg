/// @file    geo.h
/// @brief   Collect all the geometry header files together

#ifndef LARCOREALG_GEOMETRY_GEO_H
#define LARCOREALG_GEOMETRY_GEO_H

#include "cetlib_except/exception.h"

#include <stdexcept>

/// Detector geometry definition and interface
namespace geo {
  void ProjectToBoxEdge(double const xyz[],
                        double const dxyz[],
                        double xlo,
                        double xhi,
                        double ylo,
                        double yhi,
                        double zlo,
                        double zhi,
                        double xyzout[]);
  double ClosestApproach(double const point[],
                         double const intercept[],
                         double const slopes[],
                         double closest[]);
  bool CrossesBoundary(double x0[],
                       double gradient[],
                       double x_lo,
                       double x_hi,
                       double y_lo,
                       double y_hi,
                       double z_lo,
                       double z_hi,
                       double point[]);
}

///
/// Project along a direction from a particular starting point to the
/// edge of a box
///
/// \param xyz    - The starting x,y,z location. Must be inside box.
/// \param dxyz   - Direction vector
/// \param xlo    - Low edge of box in x
/// \param xhi    - Low edge of box in x
/// \param ylo    - Low edge of box in y
/// \param yhi    - Low edge of box in y
/// \param zlo    - Low edge of box in z
/// \param zhi    - Low edge of box in z
/// \param xyzout - On output, the position at the box edge
///
/// Note: It should be safe to use the same array for input and
/// output.
///
inline void geo::ProjectToBoxEdge(double const xyz[],
                                  double const dxyz[],
                                  double xlo,
                                  double xhi,
                                  double ylo,
                                  double yhi,
                                  double zlo,
                                  double zhi,
                                  double xyzout[])
{
  // Make sure we're inside the box!
  if (!(xyz[0] >= xlo && xyz[0] <= xhi) || !(xyz[1] >= ylo && xyz[1] <= yhi) ||
      !(xyz[2] >= zlo && xyz[2] <= zhi))
    throw cet::exception("ProjectToBoxEdge") << "desired point is not"
                                             << " in the specififed box\n";

  // Compute the distances to the x/y/z walls
  double dx = 99.E99;
  double dy = 99.E99;
  double dz = 99.E99;
  if (dxyz[0] > 0.0) { dx = (xhi - xyz[0]) / dxyz[0]; }
  else if (dxyz[0] < 0.0) {
    dx = (xlo - xyz[0]) / dxyz[0];
  }
  if (dxyz[1] > 0.0) { dy = (yhi - xyz[1]) / dxyz[1]; }
  else if (dxyz[1] < 0.0) {
    dy = (ylo - xyz[1]) / dxyz[1];
  }
  if (dxyz[2] > 0.0) { dz = (zhi - xyz[2]) / dxyz[2]; }
  else if (dxyz[2] < 0.0) {
    dz = (zlo - xyz[2]) / dxyz[2];
  }

  // Choose the shortest distance
  double d = 0.0;
  if (dx < dy && dx < dz)
    d = dx;
  else if (dy < dz && dy < dx)
    d = dy;
  else if (dz < dx && dz < dy)
    d = dz;

  // Make the step
  for (int i = 0; i < 3; ++i) {
    xyzout[i] = xyz[i] + dxyz[i] * d;
  }
}

///
/// Find the distance of closest approach between point and line
///
/// \param point - xyz coordinates of point
/// \param intercept - xyz coodinates of point on line
/// \param slopes - unit vector direction (need not be normalized)
/// \param closest - on output, point on line that is closest
///
/// \returns distance from point to line
///
inline double geo::ClosestApproach(double const point[],
                                   double const intercept[],
                                   double const slopes[],
                                   double closest[])
{
  double s = (slopes[0] * (point[0] - intercept[0]) + slopes[1] * (point[1] - intercept[1]) +
              slopes[2] * (point[2] - intercept[2]));
  double sd = (slopes[0] * slopes[0] + slopes[1] * slopes[1] + slopes[2] * slopes[2]);
  if (sd > 0.0) {
    s /= sd;
    closest[0] = intercept[0] + s * slopes[0];
    closest[1] = intercept[1] + s * slopes[1];
    closest[2] = intercept[2] + s * slopes[2];
  }
  else {
    // How to handle this zero gracefully? Assume that the intercept
    // is a particle vertex and "slopes" are momenta. In that case,
    // the particle goes nowhere and the closest approach is the
    // distance from the intercept to point
    closest[0] = intercept[0];
    closest[1] = intercept[1];
    closest[2] = intercept[2];
  }
  return std::sqrt(pow((point[0] - closest[0]), 2) + pow((point[1] - closest[1]), 2) +
                   pow((point[2] - closest[2]), 2));
}

///
/// Determine whether or not track intersects box of volume:
///   ( x_hi - x_lo ) x ( y_hi - y_lo ) x ( z_hi - z_lo )
///
/// \param x_hi - x box coordinates in space w.r.t. the origin
/// \param x_lo - x box coordinates in space w.r.t. the origin
/// \param y_hi - y box coordinates in space w.r.t. the origin
/// \param y_lo - y box coordinates in space w.r.t. the origin
/// \param z_hi - z box coordinates in space w.r.t. the origin
/// \param z_lo - z box coordinates in space w.r.t. the origin
/// \param x0[] - initial position of the particle
/// \param gradient[] - initial gradient of particle position
/// \param point[] (_output_) - position of the particle at intersection
///
/// @see `geo::BoxBoundedGeo::GetIntersection()`
///
/// *** assumes particle's track is linear
///
inline bool geo::CrossesBoundary(double x0[],       // initial particle position
                                 double gradient[], // initial particle gradient
                                 double x_lo,       // -
                                 double x_hi,       //  |
                                 double y_lo,       //  |- box coordinates
                                 double y_hi,       //  |  (make into vectors?)
                                 double z_lo,       //  |
                                 double z_hi,       // -
                                 double point[])
{

  double distance[3]; // distance to plane

  // puts box coordinates into more useful vectors (for loop later)
  double lo[3] = {x_lo, y_lo, z_lo};
  double hi[3] = {x_hi, y_hi, z_hi};

  // puts box coordinates into more useful vector (for loop later)
  double facecoord[6] = {lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]};

  int intersect[6] = {0, 0, 0, 0, 0, 0}; // initialize intersection tally vector
  int count = 0;
  // iterates through spatial axes (0,1,2) = (x,y,z)
  for (int i = 0; i < 3; i++) {
    // try both planes with normal parallel to axis
    for (int p = 0; p < 2; p++) {

      point[i] = facecoord[count];    // point on face
      distance[i] = point[i] - x0[i]; // calculate x-coordinate of track distance

      double C_dg = distance[i] / gradient[i]; // calculate correlation b/w gradient and distance

      for (int m = 0; m < 3; m++) {
        distance[m] = C_dg * gradient[m];
      }
      for (int n = 0; n < 3; n++) {
        point[n] = x0[n] + distance[n];
      }

      int j, k;
      switch (i) {
      case 0:
        j = 1;
        k = 2;
        break;
      case 1:
        j = 2;
        k = 0;
        break;
      case 2:
        j = 0;
        k = 1;
        break;
      default: throw std::logic_error("Big trouble");
      } // switch

      // now want to check to see if the point is in the right plane
      if (lo[j] < point[j] && point[j] < hi[j] && lo[k] < point[k] && point[k] < hi[k]) {

        //           double length = std::sqrt( distance[0]*distance[0]
        //                               + distance[1]*distance[1]
        //                               + distance[2]*distance[2] );

        // direction of motion w.r.t. start point
        int direction = distance[i] * gradient[i] /
                        std::sqrt((distance[i] * distance[i]) * (gradient[i] * gradient[i]));
        bool directed = (direction + 1) / 2;

        // checks if particle passes through face
        // and also checks to see whether it passes inward or outward
        //if ( track_length > length && directed ) {
        if (directed) {
          int normal = pow(-1, count + 1);
          int thru = normal * gradient[i] / std::sqrt(gradient[i] * gradient[i]);
          intersect[count] = thru;
        }
      }
      count++;
    }
  }

  // count faces it passes through,
  // ... not necessary now, but maybe useful in the future
  int passes = 0;
  for (int face = 0; face < 6; ++face) {
    passes += (intersect[face] * intersect[face]);
  }

  if (passes == 0) { return 0; }
  else {
    return 1;
  }
}

#endif // LARCOREALG_GEOMETRY_GEO_H
