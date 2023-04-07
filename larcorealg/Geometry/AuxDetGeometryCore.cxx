/**
 * @file   AuxDetGeometryCore.cxx
 * @brief  Access the description of auxiliary detector geometry - implementation file
 * @author brebel@fnal.gov
 * @see    AuxDetGeometryCore.h
 *
 */

// class header
#include "larcorealg/Geometry/AuxDetGeometryCore.h"

// lar includes
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetWireReadoutGeom.h"
#include "larcorealg/Geometry/GeoNodePath.h"
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/GeometryBuilderStandard.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include <TGeoManager.h>

// C/C++ includes
#include <algorithm> // std::for_each(), std::transform()
#include <cctype>    // ::tolower()
#include <cstddef>   // size_t
#include <memory>    // std::default_deleter<>
#include <string>
#include <utility> // std::swap()
#include <vector>

namespace geo {

  //......................................................................
  // Constructor.
  AuxDetGeometryCore::AuxDetGeometryCore(fhicl::ParameterSet const& pset)
    : fDetectorName(pset.get<std::string>("Name"))
    , fBuilderParameters(pset.get<fhicl::ParameterSet>("Builder", fhicl::ParameterSet()))
  {
    std::transform(fDetectorName.begin(), fDetectorName.end(), fDetectorName.begin(), ::tolower);
  }

  //......................................................................
  void AuxDetGeometryCore::ApplyChannelMap(std::unique_ptr<AuxDetWireReadoutGeom> pChannelMap)
  {
    pChannelMap->Initialize(fGeoData);
    fWireReadoutGeom = move(pChannelMap);
  }

  //......................................................................
  void AuxDetGeometryCore::LoadGeometryFile(std::string gdmlfile, std::string rootfile)
  {

    if (gdmlfile.empty()) {
      throw cet::exception("AuxDetGeometryCore") << "No GDML Geometry file specified!\n";
    }

    if (rootfile.empty()) {
      throw cet::exception("AuxDetGeometryCore") << "No ROOT Geometry file specified!\n";
    }

    ClearGeometry();

    // Open the GDML file, and convert it into ROOT TGeoManager format.
    // try to be efficient - if the GeometryCore object already imported
    // the file, then the gGeoManager will be non-null.  If not, import it.
    // Then lock the gGeoManager to prevent future imports.
    if (!gGeoManager) {
      // [20210701, petrillo@slac.stanford.edu]
      // same code, same comment as in `geo::GeometryCore::LoadGeometryFile()`.
      TGeoManager::LockDefaultUnits(false);
      TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);
      TGeoManager::LockDefaultUnits(true);
      TGeoManager::Import(rootfile.c_str());
      gGeoManager->LockGeometry();
    }

    GeometryBuilderStandard builder(
      fhicl::Table<GeometryBuilderStandard::Config>(fBuilderParameters, {"tool_type"}));
    GeoNodePath path{gGeoManager->GetTopNode()};

    AuxDets() = builder.extractAuxiliaryDetectors(path);

    fGDMLfile = gdmlfile;
    fROOTfile = rootfile;

    mf::LogInfo("AuxDetGeometryCore") << "New detector geometry loaded from "
                                      << "\n\t" << fROOTfile << "\n\t" << fGDMLfile << "\n";

  } // AuxDetGeometryCore::LoadGeometryFile()

  //......................................................................
  void AuxDetGeometryCore::ClearGeometry() { AuxDets().clear(); }

  //......................................................................
  unsigned int AuxDetGeometryCore::NAuxDetSensitive(size_t const& aid) const
  {
    if (aid > NAuxDets() - 1)
      throw cet::exception("Geometry")
        << "Requested AuxDet index " << aid << " is out of range: " << NAuxDets();

    return AuxDets()[aid].NSensitiveVolume();
  }

  //......................................................................
  //
  // Return the geometry description of the ith AuxDet.
  //
  // \param ad : input AuxDet number, starting from 0
  // \returns AuxDet geometry for ith AuxDet
  //
  // \throws geo::Exception if "ad" is outside allowed range
  //
  AuxDetGeo const& AuxDetGeometryCore::AuxDet(unsigned int const ad) const
  {
    if (ad >= NAuxDets())
      throw cet::exception("AuxDetGeometryCore") << "AuxDet " << ad << " does not exist\n";

    return AuxDets()[ad];
  }

  //......................................................................
  unsigned int AuxDetGeometryCore::FindAuxDetAtPosition(Point_t const& worldPos,
                                                        double tolerance) const
  {
    return fWireReadoutGeom->NearestAuxDet(worldPos, AuxDets(), tolerance);
  }

  //......................................................................
  AuxDetGeo const& AuxDetGeometryCore::PositionToAuxDet(Point_t const& worldLoc,
                                                        unsigned int& ad,
                                                        double tolerance) const
  {
    // locate the desired Auxiliary Detector
    ad = FindAuxDetAtPosition(worldLoc, tolerance);
    return AuxDet(ad);
  }

  //......................................................................
  void AuxDetGeometryCore::FindAuxDetSensitiveAtPosition(Point_t const& worldPos,
                                                         size_t& adg,
                                                         size_t& sv,
                                                         double tolerance) const
  {
    adg = FindAuxDetAtPosition(worldPos, tolerance);
    sv = fWireReadoutGeom->NearestSensitiveAuxDet(worldPos, AuxDets(), adg, tolerance);
  }

  //......................................................................
  AuxDetSensitiveGeo const& AuxDetGeometryCore::PositionToAuxDetSensitive(Point_t const& worldLoc,
                                                                          size_t& ad,
                                                                          size_t& sv,
                                                                          double tolerance) const
  {
    // locate the desired Auxiliary Detector
    FindAuxDetSensitiveAtPosition(worldLoc, ad, sv, tolerance);
    return AuxDet(ad).SensitiveVolume(sv);
  }

  //......................................................................
  uint32_t AuxDetGeometryCore::PositionToAuxDetChannel(Point_t const& worldLoc,
                                                       size_t& ad,
                                                       size_t& sv) const
  {
    return fWireReadoutGeom->PositionToAuxDetChannel(worldLoc, AuxDets(), ad, sv);
  }

  //......................................................................
  Point_t AuxDetGeometryCore::AuxDetChannelToPosition(std::string const& auxDetName,
                                                      uint32_t const channel) const
  {
    return fWireReadoutGeom->AuxDetChannelToPosition(channel, auxDetName, AuxDets());
  }

  //......................................................................
  AuxDetGeo const& AuxDetGeometryCore::ChannelToAuxDet(std::string const& auxDetName,
                                                       uint32_t const channel) const
  {
    size_t adIdx = fWireReadoutGeom->ChannelToAuxDet(AuxDets(), auxDetName, channel);
    return AuxDet(adIdx);
  }

  //......................................................................
  AuxDetSensitiveGeo const& AuxDetGeometryCore::ChannelToAuxDetSensitive(
    std::string const& auxDetName,
    uint32_t const channel) const
  {
    auto idx = fWireReadoutGeom->ChannelToSensitiveAuxDet(AuxDets(), auxDetName, channel);
    return AuxDet(idx.first).SensitiveVolume(idx.second);
  }

  //......................................................................

} // namespace geo
