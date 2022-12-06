/**
 * @file   StandaloneGeometrySetup.cxx
 * @brief  Utilities for one-line geometry initialization.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 22, 2017
 */

#include "larcorealg/Geometry/StandaloneGeometrySetup.h"

// LArSoft libraries
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"

// CET libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT libraries
#include "TGeoManager.h"

// C/C++ standard libraries
#include <memory> // std::make_unique(), std::make_shared()
#include <string>
#include <utility> // std::forward()

//------------------------------------------------------------------------------
std::unique_ptr<geo::GeometryCore> lar::standalone::GeometryFor(
  fhicl::ParameterSet const& pset,
  std::unique_ptr<geo::GeoObjectSorter const> sorter)
{
  auto const bForceReload = true;

  //
  // create the geometry object
  //
  auto geom = std::make_unique<geo::GeometryCore>(pset, std::move(sorter));

  //
  // extract of relevant configuration parameters
  //
  std::string relPath = pset.get<std::string>("RelativePath", "");
  const bool disableWiresInG4 = pset.get<bool>("DisableWiresInG4", false);
  const std::string GDMLFileName = pset.get<std::string>("GDML");
  //  const std::string ROOTFileName   = pset.get<std::string>("ROOT"                   );

  // add a final directory separator ("/") to relPath if not already there
  if (!relPath.empty() && (relPath.back() != '/')) relPath += '/';

  // We are going to find files now.
  // cet::search_path constructor decides if the constructor argument is a path
  // or an environment variable (in this case, the latter)
  cet::search_path sp("FW_SEARCH_PATH");

  //
  // "GDML" file (for GEANT4)
  //
  // this is our hint for the path; start with the relative path:
  std::string GDMLFilePathHint = relPath + GDMLFileName;

  // special if geometry with no wires is used for GEANT4 simulation
  if (disableWiresInG4) {
    GDMLFilePathHint.insert(std::min(GDMLFilePathHint.rfind(".gdml"), GDMLFilePathHint.length()),
                            "_nowires");
  } // if disable wires

  std::string GDMLFilePath;
  if (!sp.find_file(GDMLFilePathHint, GDMLFilePath)) {
    throw cet::exception("StaticLoadGeometry")
      << "Can't find geometry file '" << GDMLFilePathHint << "' (for GEANT4)!\n";
  }

  //
  // "ROOT" file (for geometry)
  //
  // this is our hint for the path; start with the relative path:
  std::string ROOTFilePathHint = relPath + GDMLFileName;

  std::string ROOTFilePath;
  if (!sp.find_file(ROOTFilePathHint, ROOTFilePath)) {
    throw cet::exception("StaticLoadGeometry")
      << "Can't find geometry file '" << ROOTFilePathHint << "' (for geometry)!\n";
  }

  //
  // initialize the geometry with the files we have found
  //
  geom->LoadGeometryFile(GDMLFilePath, ROOTFilePath, bForceReload);

  return geom;
} // lar::standalone::SetupGeometryWithChannelMapping()

//------------------------------------------------------------------------------
