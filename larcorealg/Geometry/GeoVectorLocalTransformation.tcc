/**
 * @file   larcorealg/Geometry/GeoVectorLocalTransformation.tcc
 * @brief  Specialization of local-to-world transformations for ROOT GenVector
 *         (template implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 31, 2019
 * @see    `GeoVectorLocalTransformation.h`
 *
 * This file is expected to be included directly in the header.
 *
 */

#ifndef LARCOREALG_GEOMETRY_GEOVECTORLOCALTRANSFORMATION_TCC
#define LARCOREALG_GEOMETRY_GEOVECTORLOCALTRANSFORMATION_TCC

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::fillCoords()...

// framework libraries
#include "cetlib_except/exception.h"

// ROOT
#include "TGeoMatrix.h"
#include "TGeoNode.h"

// CLHEP
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/Rotation.h"           // CLHEP::HepRotation
#include "CLHEP/Vector/RotationInterfaces.h" // CLHEP::HepRep3x3
#include "CLHEP/Vector/ThreeVector.h"        // CLHEP::Hep3Vector

// C standard library
#include <cassert>

//------------------------------------------------------------------------------
template <>
struct geo::details::TransformationMatrixConverter<ROOT::Math::Transform3D, TGeoMatrix> {
  static ROOT::Math::Transform3D convert(TGeoMatrix const& trans);
};

//------------------------------------------------------------------------------
template <>
template <typename DestPoint, typename SrcPoint>
DestPoint geo::LocalTransformation<ROOT::Math::Transform3D>::WorldToLocalImpl(
  SrcPoint const& world) const
{
  // need inverse transformation
  using geo::vect::convertTo;
  return convertTo<DestPoint>(
    fGeoMatrix.ApplyInverse(convertTo<typename TransformationMatrix_t::Point>(world)));
}

//......................................................................
template <>
template <typename DestVector, typename SrcVector>
DestVector geo::LocalTransformation<ROOT::Math::Transform3D>::WorldToLocalVectImpl(
  SrcVector const& world) const
{
  // need inverse transformation
  using geo::vect::convertTo;
  return convertTo<DestVector>(
    fGeoMatrix.ApplyInverse(convertTo<typename TransformationMatrix_t::Vector>(world)));
}

//......................................................................
template <>
template <typename DestPoint, typename SrcPoint>
DestPoint geo::LocalTransformation<ROOT::Math::Transform3D>::LocalToWorldImpl(
  SrcPoint const& local) const
{
  // need direct transformation
  using geo::vect::convertTo;
  return convertTo<DestPoint>(fGeoMatrix(convertTo<typename TransformationMatrix_t::Point>(local)));
}

//......................................................................
template <>
template <typename DestVector, typename SrcVector>
DestVector geo::LocalTransformation<ROOT::Math::Transform3D>::LocalToWorldVectImpl(
  SrcVector const& local) const
{
  // need direct transformation
  using geo::vect::convertTo;
  return convertTo<DestVector>(
    fGeoMatrix(convertTo<typename TransformationMatrix_t::Vector>(local)));
}

//------------------------------------------------------------------------------
template <>
inline ROOT::Math::Transform3D geo::transformationFromPath<ROOT::Math::Transform3D>(
  GeoNodeIterator_t begin,
  GeoNodeIterator_t end)
{
  if (begin == end) return {}; // identity by default construction
  auto iNode = begin;
  ROOT::Math::Transform3D matrix =
    convertTransformationMatrix<ROOT::Math::Transform3D>(*((*iNode)->GetMatrix()));
  while (++iNode != end) {
    matrix *= convertTransformationMatrix<ROOT::Math::Transform3D>(*(*iNode)->GetMatrix());
  }
  return matrix;
}

//------------------------------------------------------------------------------
template <>
inline ROOT::Math::Transform3D geo::transformationFromPath<ROOT::Math::Transform3D>(
  std::vector<TGeoNode const*> const& path,
  size_t depth)
{
  return transformationFromPath<ROOT::Math::Transform3D>(path.begin(), path.begin() + depth + 1);
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_GEOVECTORLOCALTRANSFORMATION_TCC

// Local variables:
// mode: c++
// End:
