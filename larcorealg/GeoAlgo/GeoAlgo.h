/**
 * \file GeoAlgo.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class GeoAlgo
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOALGO_H
#define BASICTOOL_GEOALGO_H

#include "larcorealg/GeoAlgo/GeoAABox.h"
#include "larcorealg/GeoAlgo/GeoHalfLine.h"
#include "larcorealg/GeoAlgo/GeoLine.h"
#include "larcorealg/GeoAlgo/GeoLineSegment.h"
#include "larcorealg/GeoAlgo/GeoSphere.h"
#include "larcorealg/GeoAlgo/GeoTrajectory.h"
#include "larcorealg/GeoAlgo/GeoVector.h"

#include <vector>

namespace geoalgo {

  /**
     \class GeoAlgo
     @brief Algorithm to compute various geometrical relation among geometrical objects.
     In particular functions to inspect following relations are implemented: \n
     0) Distance between geometrical objects \n
     1) Closest point of approach            \n
     2) Intersection points                  \n
     3) Containment/Overlap of objects       \n
     4) Common Origin functions              \n
     5) Bounding Sphere functions            \n

     Most functions are taken from the reference Real-Time-Collision-Detection (RTCD):
     Ref: http://realtimecollisiondetection.net
  */
  class GeoAlgo {

  public:
    //
    // Intersections
    //

    /// Intersection between a HalfLine and an AABox
    std::vector<Point_t> Intersection(AABox_t const& box,
                                      HalfLine_t const& line,
                                      bool back = false) const;
    /// Intersection between a HalfLine and an AABox
    std::vector<Point_t> Intersection(HalfLine_t const& line,
                                      AABox_t const& box,
                                      bool back = false) const
    {
      return Intersection(box, line, back);
    }

    /// Intersection between LineSegment and an AABox
    std::vector<Point_t> Intersection(AABox_t const& box, LineSegment_t const& l) const;
    /// Intersection between LineSegment and an AABox
    std::vector<Point_t> Intersection(LineSegment_t const& l, AABox_t const& box) const
    {
      return Intersection(box, l);
    }

    /// Intersection between Trajectory and an AABox
    std::vector<Point_t> Intersection(AABox_t const& box, Trajectory_t const& trj) const;
    /// Intersection between Trajectory and an AABox
    std::vector<Point_t> Intersection(Trajectory_t const& trj, AABox_t const& box) const
    {
      return Intersection(box, trj);
    }

    /// LineSegment sub-segment of HalfLine inside an AABox
    LineSegment_t BoxOverlap(AABox_t const& box, HalfLine_t const& line) const;
    /// LineSegment sub-segment of HalfLine inside an AABox
    LineSegment_t BoxOverlap(HalfLine_t const& line, AABox_t const& box) const
    {
      return BoxOverlap(box, line);
    }

    /// Get Trajectory inside box given some input trajectory -> now assumes trajectory cannot exit and re-enter box
    Trajectory_t BoxOverlap(AABox_t const& box, Trajectory_t const& trj) const;
    /// Get Trajectory inside box given some input trajectory -> now assumes trajectory cannot exit and re-enter box
    Trajectory_t BoxOverlap(Trajectory_t const& trj, AABox_t const& box) const
    {
      return BoxOverlap(box, trj);
    }

    //************************************************
    //CLOSEST APPROACH BETWEEN POINT AND INFINITE LINE
    //************************************************
    // min distance between point and infinite line
    double SqDist(Line_t const& line, Point_t const& pt) const
    {
      pt.compat(line.Pt1());
      return _SqDist_(line, pt);
    }
    // min distance between point and infinite line
    double SqDist(Point_t const& pt, Line_t const& line) const { return SqDist(line, pt); }
    // closest point (on infinite line) from point
    Point_t ClosestPt(Line_t const& line, Point_t const& pt) const
    {
      pt.compat(line.Pt1());
      return _ClosestPt_(pt, line);
    }
    // closest point (on infinite line) from point
    Point_t ClosestPt(Point_t const& pt, Line_t const& line) const { return _ClosestPt_(pt, line); }

    //*******************************************
    //CLOSEST APPROACH BETWEEN TWO INFINITE LINES
    //*******************************************
    // Closest approach between two infinite line segments - keep track of closest approach points
    double SqDist(Line_t const& l1, Line_t const& l2, Point_t& L1, Point_t& L2) const
    {
      l1.Pt1().compat(l2.Pt1());
      return _SqDist_(l1, l2, L1, L2);
    }
    // Closest approach between two infinite line segments - don't keep track of closest approach points
    double SqDist(Line_t const& l1, Line_t const& l2) const
    {
      Point_t L1;
      Point_t L2;
      return SqDist(l1, l2, L1, L2);
    }

    //************************************************
    //CLOSEST APPROACH BETWEEN TWO HALF-INFINITE LINES
    //************************************************
    // Closest approach between two infinite line segments - keep track of closest approach points
    double SqDist(HalfLine_t const& l1, HalfLine_t const& l2, Point_t& L1, Point_t& L2) const
    {
      l1.Start().compat(l2.Start());
      return _SqDist_(l1, l2, L1, L2);
    }
    // Closest approach between two infinite line segments - don't keep track of closest approach points
    double SqDist(HalfLine_t const& l1, HalfLine_t const& l2) const
    {
      Point_t L1;
      Point_t L2;
      return SqDist(l1, l2, L1, L2);
    }

    //******************************************
    //CLOSEST APPROACH BETWEEN TWO LINE SEGMENTS
    //******************************************
    /// LineSegment_t & LineSegment_t distance - keep track of points
    double SqDist(LineSegment_t const& seg1,
                  LineSegment_t const& seg2,
                  Point_t& c1,
                  Point_t& c2) const
    {
      seg1.Start().compat(seg2.Start());
      return _SqDist_(seg1, seg2, c1, c2);
    }
    /// LineSegment & LineSegment, don't keep track of points
    double SqDist(LineSegment_t const& seg1, LineSegment_t const& seg2) const
    {
      Point_t c1;
      Point_t c2;
      return SqDist(seg1, seg2, c1, c2);
    }

    //******************************************
    //CLOSEST APPROACH BETWEEN SEGMENT AND TRACK
    //******************************************
    /// LineSegment & Trajectory, keep track of points
    double SqDist(LineSegment_t const& seg,
                  Trajectory_t const& trj,
                  Point_t& c1,
                  Point_t& c2) const;
    /// LineSegment & Trajectory, keep track of points
    double SqDist(Trajectory_t const& trj, LineSegment_t const& seg, Point_t& c1, Point_t& c2) const
    {
      return SqDist(seg, trj, c1, c2);
    }
    /// LineSegment & Trajectory, don't keep track of points
    double SqDist(Trajectory_t const& trj, LineSegment_t const& seg) const
    {
      Point_t c1;
      Point_t c2;
      return SqDist(seg, trj, c1, c2);
    }
    /// LineSegment & Trajectory, don't keep track of points
    double SqDist(LineSegment_t const& seg, Trajectory_t const& trj) const
    {
      Point_t c1;
      Point_t c2;
      return SqDist(seg, trj, c1, c2);
    }

    //****************************************
    //CLOSEST APPROACH BETWEEN TRACK AND TRACK
    //****************************************
    /// Trajectory & Trajectory, keep track of points
    double SqDist(Trajectory_t const& trj1,
                  Trajectory_t const& trj2,
                  Point_t& c1,
                  Point_t& c2) const;
    /// Trajectory & Trajectory, don't keep track of points
    double SqDist(Trajectory_t const& trj1, Trajectory_t const& trj2) const
    {
      Point_t c1;
      Point_t c2;
      return SqDist(trj1, trj2, c1, c2);
    }

    //*****************************************************
    //CLOSEST APPROACH BETWEEN SEGMENT AND VECTOR OF TRACKS
    //*****************************************************
    /// LineSegment & vector of Trajectories, keep track of points
    double SqDist(LineSegment_t const& seg,
                  const std::vector<geoalgo::Trajectory_t>& trj,
                  Point_t& c1,
                  Point_t& c2,
                  int& trackIdx) const;
    /// LineSegment & vector of Trajectories, keep track of points
    double SqDist(const std::vector<geoalgo::Trajectory_t>& trj,
                  LineSegment_t const& seg,
                  Point_t& c1,
                  Point_t& c2,
                  int& trackIdx) const
    {
      return SqDist(seg, trj, c2, c1, trackIdx);
    }
    /// LineSegment & vector of Trajectories, don't keep track of points
    double SqDist(const std::vector<geoalgo::Trajectory_t>& trj, LineSegment_t const& seg) const
    {
      Point_t c1;
      Point_t c2;
      int trackIdx;
      return SqDist(seg, trj, c1, c2, trackIdx);
    }
    /// LineSegment & vector of Trajectories, don't keep track of points
    double SqDist(LineSegment_t const& seg, const std::vector<geoalgo::Trajectory_t>& trj) const
    {
      Point_t c1;
      Point_t c2;
      int trackIdx;
      return SqDist(seg, trj, c1, c2, trackIdx);
    }

    //*******************************************
    //CLOSEST APPROACH BETWEEN HALFLINE AND TRACK
    //*******************************************
    /// HalfLine & Trajectory, keep track of points
    double SqDist(HalfLine_t const& hline, Trajectory_t const& trj, Point_t& c1, Point_t& c2) const;
    /// HalfLine & Trajectory, keep track of points
    double SqDist(Trajectory_t const& trj, HalfLine_t const& hline, Point_t& c1, Point_t& c2) const
    {
      return SqDist(hline, trj, c2, c1);
    }
    /// HalfLine & Trajectory, don't keep track of points
    double SqDist(Trajectory_t const& trj, HalfLine_t const& hline) const
    {
      Point_t c1;
      Point_t c2;
      return SqDist(hline, trj, c1, c2);
    }
    /// HalfLine & Trajectory, don't keep track of points
    double SqDist(HalfLine_t const& hline, Trajectory_t const& trj) const
    {
      Point_t c1;
      Point_t c2;
      return SqDist(hline, trj, c1, c2);
    }

    //****************************************
    //CLOSEST APPROACH BETWEEN POINT AND TRACK
    //****************************************
    /// Point_t & Trajectory_t distance
    double SqDist(Point_t const& pt, Trajectory_t const& trj) const;
    /// Point_t & Trajectory_t distance
    double SqDist(Trajectory_t const& trj, Point_t const& pt) const { return SqDist(pt, trj); }
    /// Point_t & Trajectory_t closest point
    Point_t ClosestPt(Point_t const& pt, Trajectory_t const& trj) const
    {
      int idx = 0;
      return ClosestPt(pt, trj, idx);
    }
    /// Point_t & Trajectory_t closest point
    Point_t ClosestPt(Trajectory_t const& trj, Point_t const& pt) const
    {
      int idx = 0;
      return ClosestPt(pt, trj, idx);
    }
    /// Point_t & Trajectory_t closest point. Keep track of index of segment
    Point_t ClosestPt(Point_t const& pt, Trajectory_t const& trj, int& idx) const;
    /// Point_t & Trajectory_t closest point. Keep track of index of segment
    Point_t ClosestPt(Trajectory_t const& trj, Point_t const& pt, int& idx) const
    {
      return ClosestPt(pt, trj, idx);
    }

    //***************************************************
    //CLOSEST APPROACH BETWEEN POINT AND VECTOR OF TRACKS
    //***************************************************
    /// Point_t & Trajectory_t distance - keep track of which track
    double SqDist(Point_t const& pt,
                  const std::vector<geoalgo::Trajectory_t>& trj,
                  int& trackIdx) const;
    /// Point_t & Trajectory_t distance - keep track of which track
    double SqDist(const std::vector<geoalgo::Trajectory_t>& trj,
                  Point_t const& pt,
                  int& trackIdx) const
    {
      return SqDist(pt, trj, trackIdx);
    }
    /// Point_t & Trajectory_t distance - don't keep track
    double SqDist(Point_t const& pt, const std::vector<geoalgo::Trajectory_t>& trj) const
    {
      int trackIdx;
      return SqDist(pt, trj, trackIdx);
    }
    /// Point_t & Trajectory_t distance - don't keep track
    double SqDist(const std::vector<geoalgo::Trajectory_t>& trj, Point_t const& pt) const
    {
      int trackIdx;
      return SqDist(pt, trj, trackIdx);
    }
    /// Point_t & Trajectory_t closest point - keep track of which track is closest
    Point_t ClosestPt(Point_t const& pt,
                      const std::vector<geoalgo::Trajectory_t>& trj,
                      int& trackIdx) const;
    /// Point_t & Trajectory_t closest point - keep track of which track is closest
    Point_t ClosestPt(const std::vector<geoalgo::Trajectory_t>& trj,
                      Point_t const& pt,
                      int& trackIdx) const
    {
      return ClosestPt(pt, trj, trackIdx);
    }
    /// Point_t & Trajectory_t closest point - don't keep track of which track is closest
    Point_t ClosestPt(Point_t const& pt, const std::vector<geoalgo::Trajectory_t>& trj) const
    {
      int trackIdx;
      return ClosestPt(pt, trj, trackIdx);
    }
    /// Point_t & Trajectory_t closest point - don't keep track of which track is closest
    Point_t ClosestPt(const std::vector<geoalgo::Trajectory_t>& trj, Point_t const& pt) const
    {
      int trackIdx;
      return ClosestPt(pt, trj, trackIdx);
    }

    //***********************************************
    //CLOSEST APPROACH BETWEEN POINT AND LINE SEGMENT
    //***********************************************
    /// Point & LineSegment_t distance
    double SqDist(Point_t const& pt, LineSegment_t const& line) const
    {
      pt.compat(line.Start());
      return _SqDist_(pt, line);
    }
    /// Point & LineSegment distance
    double SqDist(LineSegment_t const& line, Point_t const& pt) const { return SqDist(pt, line); }
    /// Point & LineSegment closest point
    Point_t ClosestPt(Point_t const& pt, LineSegment_t const& line) const
    {
      pt.compat(line.Start());
      return _ClosestPt_(pt, line);
    }
    /// Point & LineSegment closest point
    Point_t ClosestPt(LineSegment_t const& line, Point_t const& pt) const
    {
      return ClosestPt(pt, line);
    }

    //********************************************
    //CLOSEST APPROACH BETWEEN POINT AND HALF LINE
    //********************************************
    /// Point & HalfLine distance
    double SqDist(Point_t const& pt, HalfLine_t const& line) const
    {
      pt.compat(line.Start());
      return _SqDist_(pt, line);
    }
    /// Point & HalfLine distance
    double SqDist(HalfLine_t const& line, Point_t const& pt) const { return SqDist(pt, line); }
    /// Point & HalfLine closest point
    Point_t ClosestPt(Point_t const& pt, HalfLine_t const& line) const
    {
      pt.compat(line.Start());
      return _ClosestPt_(pt, line);
    }
    /// Point & HalfLine closest point
    Point_t ClosestPt(HalfLine_t const& line, Point_t const& pt) const
    {
      return ClosestPt(pt, line);
    }

    //***************************************************
    //CLOSEST APPROACH BETWEEN HALF LINE AND LINE SEGMENT
    //***************************************************
    // half-line and line-segment. keep track of closest approach points
    double SqDist(HalfLine_t const& hline, LineSegment_t const& seg, Point_t& L1, Point_t& L2) const
    {
      hline.Start().compat(seg.Start());
      return _SqDist_(hline, seg, L1, L2);
    }
    // half-line and line-segment. keep track of closest approach points
    double SqDist(LineSegment_t const& seg, HalfLine_t const& hline, Point_t& L1, Point_t& L2) const
    {
      return SqDist(hline, seg, L2, L1);
    }
    // half-line and line-segment. Do not keep track of closest approach points
    double SqDist(HalfLine_t const& hline, LineSegment_t const& seg) const
    {
      Point_t L1;
      Point_t L2;
      return SqDist(hline, seg, L1, L2);
    }
    // half-line and line-segment. Do not keep track of closest approach points
    double SqDist(LineSegment_t const& seg, HalfLine_t const& hline) const
    {
      return SqDist(hline, seg);
    }

    //***************************************************
    //CLOSEST APPROACH BETWEEN POINT AND AXIS ALIGNED BOX
    //***************************************************
    /// Point & AABox distance
    double SqDist(Point_t const& pt, AABox_t const& box) const
    {
      pt.compat(box.Min());
      return _SqDist_(pt, box);
    }
    /// Point & AABox distance
    double SqDist(AABox_t const& box, Point_t const& pt) { return SqDist(pt, box); }
    /// Point & AABox closest point
    Point_t ClosestPt(Point_t const& pt, AABox_t const& box) const
    {
      pt.compat(box.Min());
      return _ClosestPt_(pt, box);
    }
    /// Point & AABox closest point
    Point_t ClosestPt(AABox_t const& box, Point_t const& pt) const { return ClosestPt(pt, box); }

    //***************************************************************************
    //COMMON ORIGIN ALGORITHMS: DETERMINE IF TWO GEO-OBJECTS HAVE A COMMON ORIGIN
    //***************************************************************************
    /// Common origin: Line Segment & Line Segment. Do not keep track of origin
    double commonOrigin(Line_t const& lin1, Line_t const& lin2) const
    {
      Point_t origin(lin1.Pt1().size());
      return commonOrigin(lin1, lin2, origin);
    }
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double commonOrigin(Line_t const& lin1, Line_t const& lin2, Point_t& origin) const
    {
      lin1.Pt1().compat(lin2.Pt1());
      return _commonOrigin_(lin1, lin2, origin);
    }
    /// Common origin: Line Segment & Line Segment. Do not keep track of origin
    double commonOrigin(LineSegment_t const& seg1,
                        LineSegment_t const& seg2,
                        bool backwards = false) const
    {
      Point_t origin(seg1.Start().size());
      return commonOrigin(seg1, seg2, origin, backwards);
    }
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double commonOrigin(LineSegment_t const& seg1,
                        LineSegment_t const& seg2,
                        Point_t& origin,
                        bool backwards = false) const
    {
      seg1.Start().compat(seg2.Start());
      return _commonOrigin_(seg1, seg2, origin, backwards);
    }
    /// Common origin: Line Segment & Half Line. Do not keep track of origin
    double commonOrigin(HalfLine_t const& lin,
                        LineSegment_t const& seg,
                        bool backwards = false) const
    {
      Point_t origin(lin.Start().size());
      return commonOrigin(lin, seg, origin, backwards);
    }
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double commonOrigin(HalfLine_t const& lin,
                        LineSegment_t const& seg,
                        Point_t& origin,
                        bool backwards = false) const
    {
      lin.Start().compat(seg.Start());
      return _commonOrigin_(lin, seg, origin, backwards);
    }
    /// Common origin: Line Segment & Half Line. Do not keep track of origin
    double commonOrigin(LineSegment_t const& seg,
                        HalfLine_t const& lin,
                        bool backwards = false) const
    {
      Point_t origin(lin.Start().size());
      return commonOrigin(lin, seg, origin, backwards);
    }
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double commonOrigin(LineSegment_t const& seg,
                        HalfLine_t const& lin,
                        Point_t& origin,
                        bool backwards = false) const
    {
      lin.Start().compat(seg.Start());
      return _commonOrigin_(lin, seg, origin, backwards);
    }
    /// Common origin: Half Line & Half Line. Do not keep track of origin
    double commonOrigin(HalfLine_t const& lin1,
                        HalfLine_t const& lin2,
                        bool backwards = false) const
    {
      Point_t origin(lin1.Start().size());
      return commonOrigin(lin1, lin2, origin, backwards);
    }
    /// Common origin: Half Line & Half Line. Keep track of origin
    double commonOrigin(HalfLine_t const& lin1,
                        HalfLine_t const& lin2,
                        Point_t& origin,
                        bool backwards = false) const
    {
      lin1.Start().compat(lin2.Start());
      return _commonOrigin_(lin1, lin2, origin, backwards);
    }
    /// Common origin: Trajectory & Trajectory. Do not keep track of origin
    double commonOrigin(Trajectory_t const& trj1,
                        Trajectory_t const& trj2,
                        bool backwards = false) const
    {
      Point_t origin(trj1.front().size());
      return commonOrigin(trj1, trj2, origin, backwards);
    }
    /// Common origin: Trajectory & Trajectory. Keep track of origin
    double commonOrigin(Trajectory_t const& trj1,
                        Trajectory_t const& trj2,
                        Point_t& origin,
                        bool backwards = false) const
    {
      trj1.front().compat(trj2.front());
      return _commonOrigin_(trj1, trj2, origin, backwards);
    }
    /// Common origin: Trajectory & Half Line. Do not keep track of origin
    double commonOrigin(Trajectory_t const& trj,
                        HalfLine_t const& lin,
                        bool backwards = false) const
    {
      Point_t origin(trj.front().size());
      return commonOrigin(trj, lin, origin, backwards);
    }
    /// Common origin: Trajectory & Half Line. Keep track of origin
    double commonOrigin(Trajectory_t const& trj,
                        HalfLine_t const& lin,
                        Point_t& origin,
                        bool backwards = false) const
    {
      trj.front().compat(lin.Start());
      return _commonOrigin_(trj, lin, origin, backwards);
    }
    /// Common origin: Trajectory & Half Line. Do not keep track of origin
    double commonOrigin(HalfLine_t const& lin,
                        Trajectory_t const& trj,
                        bool backwards = false) const
    {
      Point_t origin(trj.front().size());
      return commonOrigin(trj, lin, origin, backwards);
    }
    /// Common origin: Trajectory & Half Line. Keep track of origin
    double commonOrigin(HalfLine_t const& lin,
                        Trajectory_t const& trj,
                        Point_t& origin,
                        bool backwards = false) const
    {
      trj.front().compat(lin.Start());
      return _commonOrigin_(trj, lin, origin, backwards);
    }
    /// Common origin: Trajectory & Line Segment. Do not keep track of origin
    double commonOrigin(Trajectory_t const& trj,
                        LineSegment_t const& seg,
                        bool backwards = false) const
    {
      Point_t origin(trj.front().size());
      return commonOrigin(trj, seg, origin, backwards);
    }
    /// Common origin: Trajectory & Line Segment. Keep track of origin
    double commonOrigin(Trajectory_t const& trj,
                        LineSegment_t const& seg,
                        Point_t& origin,
                        bool backwards = false) const
    {
      trj.front().compat(seg.Start());
      return _commonOrigin_(trj, seg, origin, backwards);
    }
    /// Common origin: Trajectory & Line Segment. Do not keep track of origin
    double commonOrigin(LineSegment_t const& seg,
                        Trajectory_t const& trj,
                        bool backwards = false) const
    {
      Point_t origin(trj.front().size());
      return commonOrigin(trj, seg, origin, backwards);
    }
    /// Common origin: Trajectory & Line Segment. Keep track of origin
    double commonOrigin(LineSegment_t const& seg,
                        Trajectory_t const& trj,
                        Point_t& origin,
                        bool backwards = false) const
    {
      trj.front().compat(seg.Start());
      return _commonOrigin_(trj, seg, origin, backwards);
    }

    //************************************************************************
    //BOUNDING SPHERE ALGORITHM: RETURN SMALLEST SPHERE THAT BOUNDS ALL POINTS
    //************************************************************************
    // Bounding Sphere problem given a vector of 3D points
    Sphere_t boundingSphere(const std::vector<Point_t>& pts) const
    {
      for (auto& p : pts) {
        pts.front().compat(p);
      }
      return _boundingSphere_(pts);
    }

  protected:
    /// Line & Line distance w/o dimensionality check
    double _SqDist_(Line_t const& l1, Line_t const& l2, Point_t& L1, Point_t& L2) const;

    /// HalfLine & HalfLine distance w/o dimensionality check
    double _SqDist_(HalfLine_t const& l1, HalfLine_t const& l2, Point_t& L1, Point_t& L2) const;

    /// Point & LineSegment distance w/o dimensionality check
    double _SqDist_(Point_t const& pt, LineSegment_t const& line) const
    {
      return _SqDist_(pt, line.Start(), line.End());
    }

    /// Point & LineSegment distance w/o dimensionality check
    double _SqDist_(Point_t const& pt, Point_t const& line_s, Point_t const& line_e) const;

    /// Point & LineSegment distance w/o dimensionality check
    double _SqDist_(LineSegment_t const& line, Point_t const& pt) const
    {
      return _SqDist_(pt, line);
    }

    /// HalfLine & LineSegment distance w/o dimensionality check
    double _SqDist_(HalfLine_t const& hline,
                    LineSegment_t const& seg,
                    Point_t& L1,
                    Point_t& L2) const;

    /// LineSegment & LineSegment distance w/o dimensionality check
    double _SqDist_(LineSegment_t const& seg1,
                    LineSegment_t const& seg2,
                    Point_t& c1,
                    Point_t& c2) const;

    // Point & LineSegment closest point w/o dimensionality check
    Point_t _ClosestPt_(Point_t const& pt, LineSegment_t const& line) const;
    // Point & LineSegment closest point w/o dimensionality check
    Point_t _ClosestPt_(LineSegment_t const& line, Point_t const& pt) const
    {
      return _ClosestPt_(pt, line);
    }

    /// Point & HalfLine distance w/o dimensionality check
    double _SqDist_(Point_t const& pt, HalfLine_t const& line) const;
    /// Point & HalfLine distance w/o dimensionality check
    double _SqDist_(HalfLine_t const& line, Point_t const& pt) const { return _SqDist_(pt, line); }
    // Point & HalfLine closest point w/o dimensionality check
    Point_t _ClosestPt_(Point_t const& pt, HalfLine_t const& line) const;
    // Point & HalfLine closest point w/o dimensionality check
    Point_t _ClosestPt_(HalfLine_t const& line, Point_t const& pt) const
    {
      return _ClosestPt_(pt, line);
    }

    // Point & InfLine closest point w/o dimensionality check
    Point_t _ClosestPt_(Line_t const& line, Point_t const& pt) const;
    // Point & InfLine closest point w/o dimensionality check
    Point_t _ClosestPt_(Point_t const& pt, Line_t const& line) const
    {
      return _ClosestPt_(line, pt);
    }
    // Point & InfLine  distance w/o dimensionality check
    double _SqDist_(Line_t const& line, Point_t const& pt) const;
    // Point & InfLine  distance w/o dimensionality check
    double _SqDist_(Point_t const& pt, Line_t const& line) const { return _SqDist_(line, pt); }

    /// Point & AABox distance w/o dimensionality check
    double _SqDist_(Point_t const& pt, AABox_t const& box) const;
    /// Point & AABox distance w/o dimensionality check
    double _SqDist_(AABox_t const& box, Point_t const& pt) const { return _SqDist_(pt, box); }

    /// Point & AABox closest point w/o dimensionality check
    Point_t _ClosestPt_(Point_t const& pt, AABox_t const& box) const;
    /// Point & AABox closest point w/o dimensionality check
    Point_t _ClosestPt_(AABox_t const& box, Point_t const& pt) const
    {
      return _ClosestPt_(pt, box);
    }

    /// Common origin: Line & Line. Keep track of origin
    double _commonOrigin_(Line_t const& lin1, Line_t const& lin2, Point_t& origin) const;
    /// Common origin: Half Line & Half Line. Keep track of origin
    double _commonOrigin_(HalfLine_t const& lin1,
                          HalfLine_t const& lin2,
                          Point_t& origin,
                          bool backwards) const;
    /// Common origin: Line Segment & Half Line. Keep track of origin
    double _commonOrigin_(HalfLine_t const& lin,
                          LineSegment_t const& seg,
                          Point_t& origin,
                          bool backwards) const;
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double _commonOrigin_(LineSegment_t const& seg1,
                          LineSegment_t const& seg2,
                          Point_t& origin,
                          bool backwards) const;
    /// Common origin: Trajectory & Trajectory. Keep track of origin
    double _commonOrigin_(Trajectory_t const& trj1,
                          Trajectory_t const& trj2,
                          Point_t& origin,
                          bool backwards) const;
    /// Common origin: Trajectory & Line Segment. Keep track of origin
    double _commonOrigin_(Trajectory_t const& trj,
                          LineSegment_t const& seg,
                          Point_t& origin,
                          bool backwards) const;
    /// Common origin: Trajectory & Half Line. Keep track of origin
    double _commonOrigin_(Trajectory_t const& trj,
                          HalfLine_t const& lin,
                          Point_t& origin,
                          bool backwards) const;

    // Bounding Sphere given a vector of points
    Sphere_t _boundingSphere_(const std::vector<Point_t>& pts) const;
    Sphere_t _RemainingPoints_(std::vector<Point_t>& remaining, Sphere_t const& thisSphere) const;
    Sphere_t _WelzlSphere_(const std::vector<Point_t>& pts,
                           int numPts,
                           std::vector<Point_t> sosPts) const;

    /// Clamp function: checks if value out of bounds
    double _Clamp_(double const n, double const min, double const max) const;

    /// Swap two points if min & max are inverted
    inline void _Swap_(double& tmin, double& tmax) const
    {
      if (tmin > tmax) std::swap(tmin, tmax);
    }
  };
}

#endif
/** @} */ // end of doxygen group
