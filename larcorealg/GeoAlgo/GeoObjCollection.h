/**
 * \file GeoObjCollection.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class GeoObjCollection
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOOBJCOLLECTION_H
#define BASICTOOL_GEOOBJCOLLECTION_H

#include "larcorealg/GeoAlgo/GeoAABox.h"
#include "larcorealg/GeoAlgo/GeoCone.h"
#include "larcorealg/GeoAlgo/GeoHalfLine.h"
#include "larcorealg/GeoAlgo/GeoLineSegment.h"
#include "larcorealg/GeoAlgo/GeoSphere.h"
#include "larcorealg/GeoAlgo/GeoTrajectory.h"
#include "larcorealg/GeoAlgo/GeoVector.h"

#include <map>
#include <stddef.h>
#include <string>
#include <vector>

namespace geoalgo {

  /**
     \class GeoObjCollection
  */
  class GeoObjCollection {

  public:
    void Clear();

    void Add(Point_t const& pt, std::string label = "", std::string c = "");

    void Add(AABox_t const& box, std::string label = "", std::string c = "");

    void Add(LineSegment_t const& seg, std::string label = "", std::string c = "");

    void Add(HalfLine_t const& seg, std::string label = "", std::string c = "");

    void Add(Trajectory_t const& trj, std::string label = "", std::string c = "");

    void Add(Cone_t const& cone, std::string label = "", std::string c = "");

    void Add(Sphere_t const& sphere, std::string label = "", std::string c = "");

    const std::vector<geoalgo::Point_t>& Point() const { return _pt_v; }
    const std::vector<std::string>& PointColor() const { return _pt_col; }

    const std::vector<geoalgo::AABox_t>& AABox() const { return _box_v; }
    const std::vector<std::string>& AABoxColor() const { return _box_col; }

    const std::vector<geoalgo::LineSegment_t>& LineSegment() const { return _seg_v; }
    const std::vector<std::string>& LineSegmentColor() const { return _seg_col; }

    const std::vector<geoalgo::HalfLine_t>& HalfLine() const { return _lin_v; }
    const std::vector<std::string>& HalfLineColor() const { return _lin_col; }

    const std::vector<geoalgo::Trajectory_t>& Trajectory() const { return _trj_v; }
    const std::vector<std::string>& TrajectoryColor() const { return _trj_col; }

    const std::vector<geoalgo::Cone_t>& Cone() const { return _cone_v; }
    const std::vector<std::string>& ConeColor() const { return _cone_col; }

    const std::vector<geoalgo::Sphere_t>& Sphere() const { return _sphere_v; }
    const std::vector<std::string>& SphereColor() const { return _sphere_col; }

    const std::map<geoalgo::Point_t, std::string>& Labels() const { return _labels; }

  protected:
    Point_t const& _Point_(size_t i) const { return _pt_v[i]; }

    AABox_t const& _AABox_(size_t i) const { return _box_v[i]; }

    LineSegment_t const& _LineSegment_(size_t i) const { return _seg_v[i]; }

    Trajectory_t const& _Trajectory_(size_t i) const { return _trj_v[i]; }

    Cone_t const& _Cone_(size_t i) const { return _cone_v[i]; }

    Sphere_t const& _Sphere_(size_t i) const { return _sphere_v[i]; }

    void _AddLabel_(Point_t const& pt, std::string label);

    std::vector<geoalgo::Point_t> _pt_v;
    std::vector<std::string> _pt_col;
    std::vector<geoalgo::AABox_t> _box_v;
    std::vector<std::string> _box_col;
    std::vector<geoalgo::LineSegment_t> _seg_v;
    std::vector<std::string> _seg_col;
    std::vector<geoalgo::HalfLine_t> _lin_v;
    std::vector<std::string> _lin_col;
    std::vector<geoalgo::Trajectory_t> _trj_v;
    std::vector<std::string> _trj_col;
    std::vector<geoalgo::Cone_t> _cone_v;
    std::vector<std::string> _cone_col;
    std::vector<geoalgo::Sphere> _sphere_v;
    std::vector<std::string> _sphere_col;
    std::map<geoalgo::Point_t, std::string> _labels;
  };

}
#endif
/** @} */ // end of doxygen group
