/**
 * @file   larcorealg/Geometry/GeoVectorLocalTransformation.cxx
 * @brief  Specialization of local-to-world transformations for ROOT GenVector.
 * @see    `larcorealg/Geometry/GeoVectorLocalTransformation.h`
 * @ingroup Geometry
 *
 * This is a header-only library.
 */

// LArSoft libraries
#include "larcorealg/Geometry/GeoVectorLocalTransformation.h"

// ROOT libraries
#include "TGeoMatrix.h"

// C++ standard library
#include <stdexcept> // std::runtime_error

//------------------------------------------------------------------------------
template <>
void geo::LocalTransformation<ROOT::Math::Transform3D>::LocalToWorld(double const* local,
                                                                     double* world) const
{
  details::checkVectorBufferOverlap(local, world);

  // need direct transformation
  auto const local_v = vect::makeFromCoords<typename TransformationMatrix_t::Point>(local);
  auto const world_v = fGeoMatrix(local_v);
  vect::fillCoords(world, world_v);
}

//------------------------------------------------------------------------------
template <>
void geo::LocalTransformation<ROOT::Math::Transform3D>::LocalToWorldVect(double const* local,
                                                                         double* world) const
{
  details::checkVectorBufferOverlap(local, world);

  // need direct transformation
  auto const local_v = vect::makeFromCoords<typename TransformationMatrix_t::Vector>(local);
  auto const world_v = fGeoMatrix(local_v);
  vect::fillCoords(world, world_v);
}

//------------------------------------------------------------------------------
template <>
void geo::LocalTransformation<ROOT::Math::Transform3D>::WorldToLocal(double const* world,
                                                                     double* local) const
{
  details::checkVectorBufferOverlap(local, world);

  // need inverse transformation
  auto const world_v = vect::makeFromCoords<typename TransformationMatrix_t::Point>(world);
  auto const local_v = fGeoMatrix.ApplyInverse(world_v);
  vect::fillCoords(local, local_v);
}

//------------------------------------------------------------------------------
template <>
void geo::LocalTransformation<ROOT::Math::Transform3D>::WorldToLocalVect(const double* world,
                                                                         double* local) const
{
  details::checkVectorBufferOverlap(local, world);

  // need inverse transformation
  auto const world_v = vect::makeFromCoords<typename TransformationMatrix_t::Vector>(world);
  auto const local_v = fGeoMatrix.ApplyInverse(world_v);
  vect::fillCoords(local, local_v);
}

//------------------------------------------------------------------------------
ROOT::Math::Transform3D
geo::details::TransformationMatrixConverter<ROOT::Math::Transform3D, TGeoMatrix>::convert(
  TGeoMatrix const& trans)
{
  double const* rot = trans.GetRotationMatrix();
  double const* transl = trans.GetTranslation();
  double const* scale = trans.GetScale();

  for (auto ptr = scale; ptr != scale + 3; ++ptr)
    if (*ptr != 1.0)
      throw std::runtime_error("Matrix with scaling can't be converted to Transform3D");

  return {rot[0],
          rot[1],
          rot[2],
          transl[0],
          rot[3],
          rot[4],
          rot[5],
          transl[1],
          rot[6],
          rot[7],
          rot[8],
          transl[2]};
}

//------------------------------------------------------------------------------
