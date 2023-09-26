/**
 * @file   StandaloneGeometrySetup.cxx
 * @brief  Utilities for one-line geometry initialization.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 22, 2017
 */

#include "larcorealg/Geometry/StandaloneGeometrySetup.h"

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderStandard.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

// Framework libraries
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory>  // std::make_unique(), std::make_shared()
#include <utility> // std::move()

//------------------------------------------------------------------------------
std::unique_ptr<geo::GeometryCore> lar::standalone::GeometryFor(
  fhicl::ParameterSet const& pset,
  std::unique_ptr<geo::GeoObjectSorter> sorter)
{
  return std::make_unique<geo::GeometryCore>(
    pset,
    std::make_unique<geo::GeometryBuilderStandard>(pset.get<fhicl::ParameterSet>("Builder", {})),
    std::move(sorter));
}

//------------------------------------------------------------------------------
std::unique_ptr<geo::AuxDetGeometryCore> lar::standalone::AuxDetGeometryFor(
  fhicl::ParameterSet const& pset,
  std::unique_ptr<geo::AuxDetGeoObjectSorter> sorter,
  std::unique_ptr<geo::AuxDetInitializer> initializer)
{
  return std::make_unique<geo::AuxDetGeometryCore>(pset, std::move(sorter), std::move(initializer));
}

//------------------------------------------------------------------------------
