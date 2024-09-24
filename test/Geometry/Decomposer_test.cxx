/**
 * @file    Decomposer_test.cxx
 * @brief   Unit test for Decomposer class
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    December 8th, 2016
 *
 * Linking information: currently, both Decomposer.h and pow.h are
 * header-only libraries; the former does require ROOT vectors. Therefore the
 * only needed external library (beside Boost) is ROOT for the vector data type.
 *
 * Run the test without arguments.
 *
 */

// Boost libraries
#define BOOST_TEST_MODULE (decomposer_test)
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/Geometry/Decomposer.h"
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h"

// ROOT libraries
#include "Math/GenVector/Cartesian2D.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector2D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "TVector2.h"
#include "TVector3.h"

// C/C++ standard libraries
#include <ostream>
#include <utility> // std::move()

using boost::test_tools::tolerance;

//------------------------------------------------------------------------------
template <typename Vector, typename Point, typename Proj>
void StandardDecomposerTest()
{
  //
  // Test including all methods of geo::Decomposer<>
  //
  using Decomposer_t = geo::Decomposer<Vector, Point, Proj>;

  using Vector3D_t = typename Decomposer_t::Vector_t;
  using Point3D_t = typename Decomposer_t::Point_t;
  using Projection_t = typename Decomposer_t::Projection_t;
  using DecomposedVector_t = typename Decomposer_t::DecomposedVector_t;
  using AffinePlaneBase_t = typename Decomposer_t::AffinePlaneBase_t;

  //
  // preparation
  //
  static Point3D_t const Origin{0.0, 0.0, 0.0};
  static Point3D_t const ReferencePoint{-5.0, 10.0, 15.0};
  static Vector3D_t const NullVector{0.0, 0.0, 0.0};
  static Vector3D_t const Xaxis{1.0, 0.0, 0.0};
  static Vector3D_t const Yaxis{0.0, 1.0, 0.0};
  static Vector3D_t const Zaxis{0.0, 0.0, 1.0};

  //
  // constructors
  //

  //
  //  Default constructor: projection on (x,y) with origin (0, 0, 0)
  //
  //    Decomposer()
  //
  //  Constructor: specifies a base (an origin and two direction vectors)
  //
  //    Decomposer(AffinePlaneBase_t&& base)
  //    Decomposer(AffinePlaneBase_t const& base)
  //

  Decomposer_t defaultBase;
  BOOST_TEST(defaultBase.MainDir() == Xaxis);
  BOOST_TEST(defaultBase.SecondaryDir() == Yaxis);
  BOOST_TEST(defaultBase.NormalDir() == Zaxis);
  BOOST_TEST(defaultBase.ReferencePoint() == Origin);

  AffinePlaneBase_t rotatedBase(ReferencePoint, Yaxis, Zaxis);

  Decomposer_t rotatedBase1_1(rotatedBase);
  BOOST_TEST(rotatedBase1_1.MainDir() == Yaxis);
  BOOST_TEST(rotatedBase1_1.SecondaryDir() == Zaxis);
  BOOST_TEST(rotatedBase1_1.NormalDir() == Xaxis);
  BOOST_TEST(rotatedBase1_1.ReferencePoint() == ReferencePoint);

  Decomposer_t rotatedBase1_2(std::move(rotatedBase));
  BOOST_TEST(rotatedBase1_2.MainDir() == Yaxis);
  BOOST_TEST(rotatedBase1_2.SecondaryDir() == Zaxis);
  BOOST_TEST(rotatedBase1_2.NormalDir() == Xaxis);
  BOOST_TEST(rotatedBase1_2.ReferencePoint() == ReferencePoint);

  //
  // setters
  //
  AffinePlaneBase_t negativeBase(ReferencePoint, Yaxis, Xaxis);

  //
  //  Change projection base
  //
  //    void SetBase(AffinePlaneBase_t&& base)
  //
  //  Change projection base
  //
  //    void SetBase(AffinePlaneBase_t const& base)
  //
  //  Change the 3D point of the reference frame origin
  //
  //    void SetOrigin(Point_t const& point)
  //
  //  Change the main direction of the projection base
  //
  //    void SetMainDir(Vector_t const& dir)
  //
  //  Change the secondary direction of the projection base
  //
  //    void SetSecondaryDir(Vector_t const& dir)
  //
  //
  defaultBase.SetBase(negativeBase);
  BOOST_TEST(defaultBase.MainDir() == Yaxis);
  BOOST_TEST(defaultBase.SecondaryDir() == Xaxis);
  BOOST_TEST(defaultBase.NormalDir() == -Zaxis);
  BOOST_TEST(defaultBase.ReferencePoint() == ReferencePoint);

  defaultBase.SetOrigin(Origin);
  BOOST_TEST(defaultBase.MainDir() == Yaxis);
  BOOST_TEST(defaultBase.SecondaryDir() == Xaxis);
  BOOST_TEST(defaultBase.NormalDir() == -Zaxis);
  BOOST_TEST(defaultBase.ReferencePoint() == Origin);

  defaultBase.SetMainDir(Zaxis);
  BOOST_TEST(defaultBase.MainDir() == Zaxis);
  BOOST_TEST(defaultBase.SecondaryDir() == Xaxis);
  BOOST_TEST(defaultBase.NormalDir() == Yaxis);
  BOOST_TEST(defaultBase.ReferencePoint() == Origin);

  defaultBase.SetSecondaryDir(Yaxis);
  BOOST_TEST(defaultBase.MainDir() == Zaxis);
  BOOST_TEST(defaultBase.SecondaryDir() == Yaxis);
  BOOST_TEST(defaultBase.NormalDir() == -Xaxis);
  BOOST_TEST(defaultBase.ReferencePoint() == Origin);

  defaultBase.SetMainDir(Xaxis);
  BOOST_TEST(defaultBase.MainDir() == Xaxis);
  BOOST_TEST(defaultBase.SecondaryDir() == Yaxis);
  BOOST_TEST(defaultBase.NormalDir() == Zaxis);
  BOOST_TEST(defaultBase.ReferencePoint() == Origin);

  //
  // getters
  //

  AffinePlaneBase_t decompBase(ReferencePoint, Yaxis, Zaxis);
  Decomposer_t decomp(decompBase);

  //  Returns the reference point for the plane coordinate, as a 3D point
  //
  //    Point_t ReferencePoint() const
  //
  //  Returns the base of the decomposition
  //
  //    AffinePlaneBase_t const& Base() const
  //
  //  Returns the plane main axis direction
  //
  //    Vector_t const& MainDir() const
  //
  //  Returns the plane secondary axis direction
  //
  //    Vector_t const& SecondaryDir() const
  //
  //  Returns the plane normal axis direction
  //
  //    Vector_t const& NormalDir() const
  //

  BOOST_TEST(decomp.ReferencePoint() == ReferencePoint);
  BOOST_TEST(decomp.MainDir() == Yaxis);
  BOOST_TEST(decomp.SecondaryDir() == Zaxis);
  BOOST_TEST(decomp.NormalDir() == Xaxis);
  BOOST_TEST(decomp.Base().Origin() == ReferencePoint);
  BOOST_TEST(decomp.Base().MainDir() == Yaxis);
  BOOST_TEST(decomp.Base().SecondaryDir() == Zaxis);
  BOOST_TEST(decomp.Base().NormalDir() == Xaxis);

  //
  // projections: of 3D points
  //

  // decompBase: ( y, z, x ); ReferencePoint: { -5.0, 10.0, 15.0 }
  Point3D_t point(ReferencePoint + Vector3D_t(3.0, 4.0, 5.0)); // { -2, 14, 20 }

  //
  //  Returns the main component of a point
  //
  //    auto PointMainComponent(Point_t const& point) const
  //
  //  Returns the secondary component of a point
  //
  //    auto PointSecondaryComponent(Point_t const& point) const
  //
  //  Returns the secondary component of a point
  //
  //    auto PointNormalComponent(Point_t const& point) const
  //
  //  Returns the projection of the specified point on the plane
  //
  //    Projection_t ProjectPointOnPlane(Point_t const& point) const
  //
  //  Decomposes a 3D point in two components
  //
  //    DecomposedVector_t DecomposePoint(Point_t const& point) const
  //

  auto const tol = 0.01 % tolerance();
  BOOST_TEST(decomp.PointMainComponent(point) == 4.0, tol);
  BOOST_TEST(decomp.PointSecondaryComponent(point) == 5.0, tol);
  BOOST_TEST(decomp.PointNormalComponent(point) == 3.0, tol);

  {
    auto proj = decomp.ProjectPointOnPlane(point);
    BOOST_TEST(proj.X() == 4.0, tol);
    BOOST_TEST(proj.Y() == 5.0, tol);
  }

  {
    auto decomposed = decomp.DecomposePoint(point);
    BOOST_TEST(decomposed.distance == 3.0, tol);
    BOOST_TEST(decomposed.projection.X() == 4.0, tol);
    BOOST_TEST(decomposed.projection.Y() == 5.0, tol);
  }

  //
  // projections: of 3D vectors
  //

  // decompBase: ( y, z, x ); ReferencePoint: { -5.0, 10.0, 15.0 }
  Vector3D_t vector(2.0, 3.0, 4.0);

  //
  //  Returns the main component of a vector
  //
  //    auto VectorMainComponent(Vector_t const& v) const
  //
  //  Returns the secondary component of a vector
  //
  //    auto VectorSecondaryComponent(Vector_t const& v) const
  //
  //  Returns the secondary component of a vector
  //
  //    auto VectorNormalComponent(Vector_t const& v) const
  //
  //  Returns the projection of the specified vector on the plane
  //
  //    Projection_t ProjectVectorOnPlane(Point_t const& v) const
  //
  //  Decomposes a 3D vector in two components
  //
  //    DecomposedVector_t DecomposeVector(Vector_t const& v) const
  //
  //  Returns the main component of a projection vector
  //
  //    auto MainComponent(Projection_t const& v) const
  //
  //  Returns the secondary component of a projection vector
  //
  //    auto SecondaryComponent(Projection_t const& v) const
  //

  BOOST_TEST(decomp.VectorMainComponent(vector) == 3.0, tol);
  BOOST_TEST(decomp.VectorSecondaryComponent(vector) == 4.0, tol);
  BOOST_TEST(decomp.VectorNormalComponent(vector) == 2.0, tol);

  {
    auto proj = decomp.ProjectVectorOnPlane(vector);
    BOOST_TEST(proj.X() == 3.0, tol);
    BOOST_TEST(proj.Y() == 4.0, tol);
  }

  {
    auto decomposed = decomp.DecomposeVector(vector);
    BOOST_TEST(decomposed.distance == 2.0, tol);
    BOOST_TEST(decomposed.projection.X() == 3.0, tol);
    BOOST_TEST(decomposed.projection.Y() == 4.0, tol);
  }

  //
  // projection: projected vector on the plane
  //
  Projection_t proj(3.0, 4.0);

  BOOST_TEST(decomp.MainComponent(proj) == 3.0, tol);
  BOOST_TEST(decomp.SecondaryComponent(proj) == 4.0, tol);

  //
  // composition: point
  //
  // decompBase: ( y, z, x ); ReferencePoint: { -5.0, 10.0, 15.0 }
  DecomposedVector_t decomposed(3.0,       // distance from plane
                                {4.0, 5.0} // projection on plane in local plane coordinates
  );

  //  Returns the 3D point from composition of projection and distance
  //
  //    Point_t ComposePoint(DecomposedVector_t const& decomp) const
  //
  //  Returns the 3D point from composition of projection and distance
  //
  //    Point_t ComposePoint(double distance, Projection_t const& proj) const
  //
  //  Returns the 3D vector from composition of projection and distance
  //
  //    Vector_t ComposeVector(DecomposedVector_t const& decomp) const
  //
  //  Returns the 3D vector from composition of projection and distance
  //
  //    Vector_t ComposeVector(double distance, Projection_t const& proj) const
  //

  {
    auto comp = decomp.ComposePoint(decomposed.distance, decomposed.projection);
    BOOST_TEST(comp.X() == -2.0, tol);
    BOOST_TEST(comp.Y() == 14.0, tol);
    BOOST_TEST(comp.Z() == 20.0, tol);
  }

  {
    auto comp = decomp.ComposePoint(decomposed);
    BOOST_TEST(comp.X() == -2.0, tol);
    BOOST_TEST(comp.Y() == 14.0, tol);
    BOOST_TEST(comp.Z() == 20.0, tol);
  }

  //
  // composition: vector
  //

  {
    auto comp = decomp.ComposeVector(decomposed.distance, decomposed.projection);
    BOOST_TEST(comp.X() == 3.0, tol);
    BOOST_TEST(comp.Y() == 4.0, tol);
    BOOST_TEST(comp.Z() == 5.0, tol);
  }

  {
    auto comp = decomp.ComposeVector(decomposed);
    BOOST_TEST(comp.X() == 3.0, tol);
    BOOST_TEST(comp.Y() == 4.0, tol);
    BOOST_TEST(comp.Z() == 5.0, tol);
  }

} // StandardDecomposerTest<>()

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(TVectorDecomposerTestCase)
{

  StandardDecomposerTest<TVector3, TVector3, TVector2>();

} // BOOST_AUTO_TEST_CASE(TVectorDecomposerTestCase)

BOOST_AUTO_TEST_CASE(GenVectorDecomposerTestCase)
{

  using Point_t = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>>;
  using Vector_t = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>>;
  using Projection_t = ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<double>>;

  StandardDecomposerTest<Vector_t, Point_t, Projection_t>();

} // BOOST_AUTO_TEST_CASE(GenVectorDecomposerTestCase)
