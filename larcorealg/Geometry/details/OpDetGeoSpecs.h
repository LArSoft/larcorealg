/**
 * @file   larcorealg/Geometry/details/OpDetGeoSpecs.h
 * @brief  Definition of internal data types of `geo::OpDetGeo`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * 
 * This library is header-only.
 */

#ifndef LARCOREALG_GEOMETRY_DETAILS_OPDETGEOSPECS_H
#define LARCOREALG_GEOMETRY_DETAILS_OPDETGEOSPECS_H


// LArSoft libraries
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_optical_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t, ...


// -----------------------------------------------------------------------------
namespace geo::details {
  
  using OpDetGeoLocalTransformation = geo::LocalTransformationGeo
    <ROOT::Math::Transform3D, geo::OpticalPoint_t, geo::OpticalVector_t>;
  
  struct OpDetGeoSpecs;
  
} // namespace geo::details

// -----------------------------------------------------------------------------
/// Internal data for `geo::OpDetGeo` implementation.
struct geo::details::OpDetGeoSpecs {
  
  /// All dimensions of the optical detector.
  struct BoxSize_t {
    double width  { 0.0 };
    double height { 0.0 };
    double length { 0.0 };
  };
  
  geo::Point_t center; ///< Effective center of the detector.
  geo::Vector_t normalDir; ///< Normal to the sensitive surface.
  
  BoxSize_t sides; ///< Dimensions of a box around the sensitive detector.
  double centerToTip; ///< Distance between center and tip [cm]
  
}; // geo::details::OpDetGeoSpecs


// -----------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_DETAILS_OPDETGEOSPECS_H
