/**
 * @file   AuxDetGeometryCore.cxx
 * @brief  Access the description of auxiliary detector geometry - implementation file
 * @see    AuxDetGeometryCore.h
 */

// class header
#include "larcorealg/Geometry/AuxDetGeometryCore.h"

// LArSoft includes
#include "larcorealg/CoreUtils/SearchPathPlusRelative.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetReadoutGeom.h"
#include "larcorealg/Geometry/GeoNodePath.h"
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/GeometryBuilderStandard.h"
#include "larcorealg/Geometry/details/maybe_default_detector_name.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TGeoManager.h"
#include "TROOT.h"

// C/C++ includes
#include <algorithm> // std::for_each(), std::transform()
#include <cctype>    // std::tolower()
#include <cstddef>   // std::size_t
#include <memory>
#include <string>

namespace geo {

  //......................................................................
  // Constructor.
  AuxDetGeometryCore::AuxDetGeometryCore(fhicl::ParameterSet const& pset,
                                         std::unique_ptr<AuxDetGeoObjectSorter> sorter,
                                         std::unique_ptr<AuxDetInitializer> initializer)
    : fSorter{std::move(sorter)}
    , fInitializer{std::move(initializer)}
    , fGDMLfile{lar::searchPathPlusRelative(pset.get<std::string>("RelativePath", ""),
                                            pset.get<std::string>("GDML"))}
    , fDetectorName{details::maybe_default_detector_name(pset, fGDMLfile)}
    , fBuilderParameters(pset.get<fhicl::ParameterSet>("Builder", {}))
    , fThrowIfAbsent(pset.get("ThrowIfAbsent", true))
  {
    std::transform(fDetectorName.begin(), fDetectorName.end(), fDetectorName.begin(), [](auto c) {
      return std::tolower(c);
    });

    LoadGeometryFile();
  }

  //......................................................................
  void AuxDetGeometryCore::ApplyChannelMap()
  {
    auto initializers = fInitializer ? fInitializer->init(fAuxDets) : AuxDetReadoutInitializers{};
    fReadoutGeom = std::make_unique<AuxDetReadoutGeom>(std::move(initializers));
  }

  //......................................................................
  void AuxDetGeometryCore::LoadGeometryFile()
  {

    if (fGDMLfile.empty()) {
      throw cet::exception("AuxDetGeometryCore") << "No GDML Geometry file specified!\n";
    }

    (void)gROOT; // <= Can be removed once ROOT 6.26/08 is adopted

    // Open the GDML file, and convert it into ROOT TGeoManager format.  Try to be
    // efficient - if the GeometryCore object already imported the file, then the
    // gGeoManager will be non-null.  If not, import it.  Then lock the gGeoManager to
    // prevent future imports.
    if (!gGeoManager) {
      // [20210701, petrillo@slac.stanford.edu]
      // same code, same comment as in `geo::GeometryCore::LoadGeometryFile()`.
      TGeoManager::LockDefaultUnits(false);
      TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);
      TGeoManager::LockDefaultUnits(true);
      TGeoManager::Import(fGDMLfile.c_str());
      gGeoManager->LockGeometry();
    }

    GeometryBuilderStandard builder(
      fhicl::Table<GeometryBuilderStandard::Config>(fBuilderParameters, {"tool_type"}));
    GeoNodePath path{gGeoManager->GetTopNode()};

    fAuxDets = builder.extractAuxiliaryDetectors(path);
    if (fSorter) { fSorter->sort(fAuxDets); }

    ApplyChannelMap();

    mf::LogInfo("AuxDetGeometryCore") << "New detector geometry loaded from\n\t" << fGDMLfile;
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
  AuxDetGeo const& AuxDetGeometryCore::AuxDet(std::size_t const ad) const
  {
    if (ad >= NAuxDets()) {
      throw cet::exception("AuxDetGeometryCore") << "AuxDet " << ad << " does not exist\n";
    }
    return fAuxDets[ad];
  }

  //......................................................................
  std::size_t AuxDetGeometryCore::NAuxDetSensitive(std::size_t const aid) const
  {
    if (aid >= NAuxDets()) {
      throw cet::exception("AuxDetGeometry")
        << "Requested AuxDet index " << aid << " is out of range: " << NAuxDets();
    }
    return fAuxDets[aid].NSensitiveVolume();
  }

  //......................................................................
  std::size_t AuxDetGeometryCore::FindAuxDetAtPosition(Point_t const& point, double tolerance) const
  {
    return fReadoutGeom->NearestAuxDet(point, fAuxDets, tolerance, fThrowIfAbsent);
  }

  //......................................................................
  AuxDetGeo const& AuxDetGeometryCore::PositionToAuxDet(Point_t const& point,
                                                        double tolerance) const
  {
    return AuxDet(FindAuxDetAtPosition(point, tolerance));
  }

  //......................................................................
  void AuxDetGeometryCore::FindAuxDetSensitiveAtPosition(Point_t const& point,
                                                         std::size_t& adg,
                                                         std::size_t& sv,
                                                         double tolerance) const
  {
    adg = FindAuxDetAtPosition(point, tolerance);
    sv = fReadoutGeom->NearestSensitiveAuxDet(point, fAuxDets, tolerance, fThrowIfAbsent);
  }

  //......................................................................
  Point_t AuxDetGeometryCore::AuxDetChannelToPosition(std::string const& auxDetName,
                                                      uint32_t const channel) const
  {
    return fReadoutGeom->AuxDetChannelToPosition(channel, auxDetName, fAuxDets);
  }

  //......................................................................
  AuxDetSensitiveGeo const& AuxDetGeometryCore::ChannelToAuxDetSensitive(
    std::string const& auxDetName,
    uint32_t const channel) const
  {
    auto idx = fReadoutGeom->ChannelToSensitiveAuxDet(auxDetName, channel);
    return AuxDet(idx.first).SensitiveVolume(idx.second);
  }

} // namespace geo
