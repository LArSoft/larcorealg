/**
 * @file   topology_test.cxx
 * @brief  Test for neighbourhood discovery in simple geometries.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 12, 2017
 *
 * Usage: just run the executable.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate
// a main() function; Boost is pulled in by boost_unit_test_base.h
#define BOOST_TEST_MODULE topology_test

// LArSoft libraries
#include "larcorealg/Geometry/BoxBoundedGeo.h"

// utility libraries
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp> // BOOST_CHECK_CLOSE()

// C/C++ standard libraries
#include <vector>
#include <iterator> // std::back_inserter()
#include <algorithm> // std::transform()
#include <memory> // std::address_of()
#include <cstdlib> // std::size_t

//------------------------------------------------------------------------------


struct Neighbors_t {

  static constexpr std::size_t Width = 0;
  static constexpr std::size_t Depth = 1;
  static constexpr std::size_t Drift = 2;
  static constexpr std::size_t NDims = 3;
  static constexpr std::size_t Pos = 0;
  static constexpr std::size_t Neg = 1;
  static constexpr std::size_t NDirs = 2;

  std::vector<geo::BoxBoundedGeo const*> neighbors[NDims][NDirs];

}; // Neighbors_t


//------------------------------------------------------------------------------
void BoxedTopologyTest() {

  //
  // definition of the topology; something simple here: a 3x3x3 cube
  //
  double const side = 2.0;
  std::vector<geo::BoxBoundedGeo> boxes;
  boxes.reserve(27);
  for (int x = -1; x <= 1; ++x)
    for (int y = -1; y <= 1; ++y)
      for (int z = -1; z <= 1; ++z)
        boxes.emplace_back(
          side * (x - 0.5), side * (x + 0.5),
          side * (y - 0.5), side * (y + 0.5),
          side * (z - 0.5), side * (z + 0.5)
          );
  std::vector<geo::BoxBoundedGeo const*> volumes;
  std::transform(boxes.cbegin(), boxes.cend(), std::back_inserter(volumes),
    [](auto& obj){ return std::addressof(obj); });

  //
  // discovery of the neighbourhood
  //

  //
  // verification of the result
  //



} // BoxedTopologyTest()







//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(BoxedTopologyTestCase) {
  BoxedTopologyTest();
} // BoxedTopologyTestCase
