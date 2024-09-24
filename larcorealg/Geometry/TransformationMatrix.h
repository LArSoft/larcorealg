/**
 * @file   larcorealg/Geometry/TransformationMatrix.h
 * @brief  Selection of the type of transformation matrix used in geometry.
 * @see    `larcorealg/Geometry/GeoVectorLocalTransformation.h`
 * @ingroup Geometry
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_GEOMETRY_TRANSFORMATIONMATRIX_H
#define LARCOREALG_GEOMETRY_TRANSFORMATIONMATRIX_H

// LArSoft libraries
#include "larcorealg/Geometry/GeoVectorLocalTransformation.h"

// ROOT libraries
#include "Math/GenVector/Transform3D.h"

namespace geo {

  /**
   * @brief Type of transformation matrix used in geometry.
   *
   * This type is used for storing the transformation matrices of the various
   * geometry description objects (e.g. `geo::WireGeo`, `geo::OpDetGeo`, ...).
   */
  using TransformationMatrix = ROOT::Math::Transform3D;

  /// Converts a transformation matrix into a `geo::TransformationMatrix`.
  template <typename Trans>
  decltype(auto) makeTransformationMatrix(Trans&& trans)
  {
    return convertTransformationMatrix<TransformationMatrix>(trans);
  }

} // namespace geo

#endif // LARCOREALG_GEOMETRY_TRANSFORMATIONMATRIX_H
