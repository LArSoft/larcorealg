/**
 * @file   larcorealg/Geometry/GeometryCore.cxx
 * @brief  Access the description of detector geometry - implementation file
 * @author brebel@fnal.gov
 * @see    larcorealg/Geometry/GeometryCore.h
 *
 */

// class header
#include "larcorealg/Geometry/GeometryCore.h"

// lar includes
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi<>
#include "larcorealg/CoreUtils/DereferenceIterator.h" // lar::util::dereferenceIteratorLoop()
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/GeometryBuilderStandard.h"
#include "larcorealg/Geometry/Decomposer.h" // geo::vect::dot()
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h" // geo::vect

// Framework includes
#include "cetlib/pow.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TGeoVolume.h>
// #include <Rtypes.h>

// C/C++ includes
#include <cstddef> // size_t
#include <cctype> // ::tolower()
#include <cmath> // std::abs() ...
#include <sstream> // std::ostringstream
#include <vector>
#include <algorithm> // std::for_each(), std::transform()
#include <iterator> // std::back_inserter()
#include <utility> // std::swap()
#include <limits> // std::numeric_limits<>
#include <numeric> // std::accumulate

namespace geo {

  //......................................................................
  lar::util::RealComparisons<geo::Length_t> GeometryCore::coordIs{ 1e-8 };


  //......................................................................
  // Constructor.
  GeometryCore::GeometryCore(
    fhicl::ParameterSet const& pset
    )
    : fSurfaceY         (pset.get< double            >("SurfaceY"               ))
    , fDetectorName     (pset.get< std::string       >("Name"                   ))
    , fMinWireZDist     (pset.get< double            >("MinWireZDist",     3.0  ))
    , fPositionWiggle   (pset.get< double            >("PositionEpsilon",  1.e-4))
    , fBuilderParameters
      (pset.get<fhicl::ParameterSet>("Builder", fhicl::ParameterSet()))
  {
    std::transform(fDetectorName.begin(), fDetectorName.end(),
      fDetectorName.begin(), ::tolower);
  } // GeometryCore::GeometryCore()


  //......................................................................
  GeometryCore::~GeometryCore() {
    ClearGeometry();
  } // GeometryCore::~GeometryCore()


  //......................................................................
  void GeometryCore::ApplyChannelMap
    (std::unique_ptr<geo::ChannelMapAlg> pChannelMap)
  {
    SortGeometry(pChannelMap->Sorter());
    UpdateAfterSorting(); // after channel mapping has sorted objects, set their IDs
    pChannelMap->Initialize(fGeoData);
    fChannelMapAlg = move(pChannelMap);
  } // GeometryCore::ApplyChannelMap()

  //......................................................................
  void GeometryCore::LoadGeometryFile(
    std::string gdmlfile, std::string rootfile,
    geo::GeometryBuilder& builder,
    bool bForceReload /* = false*/
  ) {

    if (gdmlfile.empty()) {
      throw cet::exception("GeometryCore")
        << "No GDML Geometry file specified!\n";
    }

    if (rootfile.empty()) {
      throw cet::exception("GeometryCore")
        << "No ROOT Geometry file specified!\n";
    }

    ClearGeometry();

    // Open the GDML file, and convert it into ROOT TGeoManager format.
    // Then lock the gGeoManager to prevent future imports, for example
    // in AuxDetGeometry
    if( !gGeoManager || bForceReload ){
      if (gGeoManager) TGeoManager::UnlockGeometry();
      TGeoManager::Import(rootfile.c_str());
      gGeoManager->LockGeometry();
    }

    BuildGeometry(builder);

    fGDMLfile = gdmlfile;
    fROOTfile = rootfile;

    mf::LogInfo("GeometryCore") << "New detector geometry loaded from "
                                << "\n\t" << fROOTfile
                                << "\n\t" << fGDMLfile << "\n";

  } // GeometryCore::LoadGeometryFile()


  //......................................................................
  void GeometryCore::LoadGeometryFile(
    std::string gdmlfile, std::string rootfile,
    bool bForceReload /* = false*/
  ) {
    fhicl::Table<geo::GeometryBuilderStandard::Config> const builderConfig
      (fBuilderParameters, { "tool_type" });

    // this is a wink to the understanding that we might be using a art-based
    // service provider configuration sprinkled with tools.
    geo::GeometryBuilderStandard builder{builderConfig()};
    LoadGeometryFile(gdmlfile, rootfile, builder, bForceReload);
  } // GeometryCore::LoadGeometryFile()

  //......................................................................
  void GeometryCore::ClearGeometry()
  {
    Cryostats().clear();
    AuxDets().clear();
  }


  //......................................................................
  void GeometryCore::SortGeometry(geo::GeoObjectSorter const& sorter)
  {
    mf::LogInfo("GeometryCore") << "Sorting volumes...";

    sorter.SortAuxDets(AuxDets());
    sorter.SortCryostats(Cryostats());

    geo::CryostatID::CryostatID_t c = 0;
    for (geo::CryostatGeo& cryo: Cryostats())
    {
      cryo.SortSubVolumes(sorter);
      cryo.UpdateAfterSorting(geo::CryostatID(c));
      ++c;
    } // for

  } // GeometryCore::SortGeometry()


  //......................................................................
  void GeometryCore::UpdateAfterSorting() {

    for (size_t c = 0; c < Ncryostats(); ++c)
      Cryostats()[c].UpdateAfterSorting(geo::CryostatID(c));

    allViews.clear();
    for (geo::TPCGeo const& tpc: IterateTPCs()) {
      auto const& TPCviews = tpc.Views();
      allViews.insert(TPCviews.cbegin(), TPCviews.cend());
    }

  } // GeometryCore::UpdateAfterSorting()


  //......................................................................
  TGeoManager* GeometryCore::ROOTGeoManager() const
  {
    return gGeoManager;
  }

  //......................................................................
  unsigned int GeometryCore::Nchannels() const
  {
    return fChannelMapAlg->Nchannels();
  }

  //......................................................................
  unsigned int GeometryCore::Nchannels(readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->Nchannels(ropid);
  } // GeometryCore::Nchannels(ROPID)

  //......................................................................

  std::vector<raw::ChannelID_t> GeometryCore::ChannelsInTPCs() const
  {
    std::vector<raw::ChannelID_t> channels;
    channels.reserve(fChannelMapAlg->Nchannels());

    for (const readout::TPCsetID & ts: IterateTPCsetIDs())
    {
      for (auto const t: fChannelMapAlg->TPCsetToTPCs(ts))
      {
        for (auto const & wire: IterateWireIDs(t))
        {
          channels.push_back(fChannelMapAlg->PlaneWireToChannel(wire));
        }
      }
    }
    std::sort(channels.begin(), channels.end());
    auto last = std::unique(channels.begin(), channels.end());
    channels.erase(last, channels.end());
    return channels;
  }


  //......................................................................
  unsigned int GeometryCore::NOpDets() const
  {
    int N=0;
    for(size_t cstat=0; cstat!=Ncryostats(); cstat++)
      N += this->Cryostat(cstat).NOpDet();
    return N;
  }

  //......................................................................
  unsigned int GeometryCore::NOpChannels() const
  {
    return fChannelMapAlg->NOpChannels(this->NOpDets());
  }

  //......................................................................
  unsigned int GeometryCore::MaxOpChannel() const
  {
    return fChannelMapAlg->MaxOpChannel(this->NOpDets());
  }

  //......................................................................
  unsigned int GeometryCore::NOpHardwareChannels(int opDet) const
  {
    return fChannelMapAlg->NOpHardwareChannels(opDet);
  }

  //......................................................................
  unsigned int GeometryCore::OpChannel(int detNum, int hardwareChannel) const
  {
    return fChannelMapAlg->OpChannel(detNum, hardwareChannel);
  }

  //......................................................................
  unsigned int GeometryCore::OpDetFromOpChannel(int opChannel) const
  {
    return fChannelMapAlg->OpDetFromOpChannel(opChannel);
  }

  //......................................................................
  unsigned int GeometryCore::HardwareChannelFromOpChannel(int opChannel) const
  {
    return fChannelMapAlg->HardwareChannelFromOpChannel(opChannel);
  }

  //......................................................................
  // Is this a valid OpChannel number?
  bool GeometryCore::IsValidOpChannel(int opChannel) const
  {
    return fChannelMapAlg->IsValidOpChannel(opChannel, this->NOpDets());
  }

  //......................................................................
  unsigned int GeometryCore::NAuxDetSensitive(size_t const& aid) const
  {
    if( aid > NAuxDets() - 1)
      throw cet::exception("Geometry") << "Requested AuxDet index " << aid
                                       << " is out of range: " << NAuxDets();

    return AuxDets()[aid].NSensitiveVolume();
  }

  //......................................................................
  // Number of different views, or wire orientations
  unsigned int GeometryCore::Nviews() const
  {
    return MaxPlanes();
  }

  //......................................................................
  //
  // Return the geometry description of the ith plane in the detector.
  //
  // \param cstat : input cryostat number, starting from 0
  // \returns cryostat geometry for ith cryostat
  //
  // \throws geo::Exception if "cstat" is outside allowed range
  //
  CryostatGeo const& GeometryCore::Cryostat(CryostatID const& cryoid) const {
    CryostatGeo const* pCryo = CryostatPtr(cryoid);
    if(!pCryo) {
      throw cet::exception("GeometryCore") << "Cryostat #"
                                           << cryoid.Cryostat
                                           << " does not exist\n";
    }
    return *pCryo;
  } // GeometryCore::Cryostat(CryostatID)

  //......................................................................
  //
  // Return the geometry description of the ith AuxDet.
  //
  // \param ad : input AuxDet number, starting from 0
  // \returns AuxDet geometry for ith AuxDet
  //
  // \throws geo::Exception if "ad" is outside allowed range
  //
  const AuxDetGeo& GeometryCore::AuxDet(unsigned int const ad) const
  {
    if(ad >= NAuxDets())
    throw cet::exception("GeometryCore") << "AuxDet "
    << ad
    << " does not exist\n";

    return AuxDets()[ad];
  }


  //......................................................................
  geo::TPCID GeometryCore::FindTPCAtPosition(geo::Point_t const& point) const {

    // first find the cryostat
    geo::CryostatGeo const* cryo = PositionToCryostatPtr(point);
    if (!cryo) return {};

    // then ask it about the TPC
    geo::TPCID tpcid = cryo->PositionToTPCID(point, 1. + fPositionWiggle);
    if (tpcid) return tpcid;

    // return an invalid TPC ID with cryostat information set:
    tpcid.Cryostat = cryo->ID().Cryostat;
    tpcid.markInvalid();
    return tpcid;

  } // GeometryCore::FindTPCAtPosition()


  //......................................................................
  geo::CryostatGeo const* GeometryCore::PositionToCryostatPtr
    (geo::Point_t const& point) const
  {
    for (geo::CryostatGeo const& cryostat: IterateCryostats()) {
      if (cryostat.ContainsPosition(point, 1.0 + fPositionWiggle))
        return &cryostat;
    }
    return nullptr;
  } // GeometryCore::PositionToCryostatPtr()


  //......................................................................
  geo::CryostatID GeometryCore::PositionToCryostatID
    (geo::Point_t const& point) const
  {
    geo::CryostatGeo const* cryo = PositionToCryostatPtr(point);
    return cryo? cryo->ID(): geo::CryostatID{};
  } // GeometryCore::PositionToCryostatID()


  //......................................................................
  geo::CryostatID::CryostatID_t GeometryCore::FindCryostatAtPosition
    (geo::Point_t const& worldLoc) const
  {
    geo::CryostatGeo const* cryo = PositionToCryostatPtr(worldLoc);
    return cryo? cryo->ID().Cryostat: geo::CryostatID::InvalidID;
  } // GeometryCore::FindCryostatAtPosition(Point)


  //......................................................................
  geo::CryostatID::CryostatID_t GeometryCore::FindCryostatAtPosition
    (double const worldLoc[3]) const
  {
    return FindCryostatAtPosition(geo::vect::makePointFromCoords(worldLoc));
  } // GeometryCore::FindCryostatAtPosition(double[])


  //......................................................................
  geo::TPCGeo const* GeometryCore::PositionToTPCptr
    (geo::Point_t const& point) const
  {
    geo::CryostatGeo const* cryo = PositionToCryostatPtr(point);
    return cryo? cryo->PositionToTPCptr(point, 1. + fPositionWiggle): nullptr;
  } // GeometryCore::PositionToTPCptr()


  //......................................................................
  geo::TPCGeo const& GeometryCore::PositionToTPC
    (geo::Point_t const& point) const
  {
    geo::TPCGeo const* tpc = PositionToTPCptr(point);
    if (!tpc) {
      throw cet::exception("GeometryCore")
        << "Can't find any TPC at position " << point << "\n";
    }
    return *tpc;
  } // GeometryCore::PositionToTPC()


  //......................................................................
  TPCGeo const& GeometryCore::PositionToTPC
    (double const worldLoc[3], TPCID& tpcid) const
  {
    geo::TPCGeo const& TPC = PositionToTPC(worldLoc);
    tpcid = TPC.ID();
    return TPC;
  } // GeometryCore::PositionToTPC(double*, TPCID&)


  //......................................................................
  TPCGeo const& GeometryCore::PositionToTPC
    (double const worldLoc[3], unsigned int &tpc, unsigned int &cstat) const
  {
    geo::TPCGeo const& TPC = PositionToTPC(worldLoc);
    cstat = TPC.ID().Cryostat;
    tpc = TPC.ID().TPC;
    return TPC;
  } // GeometryCore::PositionToTPC(double*, TPCID&)


  //......................................................................
  geo::TPCID GeometryCore::PositionToTPCID(geo::Point_t const& point) const {
    geo::TPCGeo const* tpc = PositionToTPCptr(point);
    return tpc? tpc->ID(): geo::TPCID{};
  } // GeometryCore::PositionToTPCID()


  //......................................................................
  geo::CryostatGeo const& GeometryCore::PositionToCryostat
    (geo::Point_t const& point) const
  {
    geo::CryostatGeo const* cstat = PositionToCryostatPtr(point);
    if (!cstat) {
      throw cet::exception("GeometryCore")
        << "Can't find any cryostat at position " << point << "\n";
    }
    return *cstat;
  } // GeometryCore::PositionToCryostat()


  //......................................................................
  const CryostatGeo& GeometryCore::PositionToCryostat
    (double const  worldLoc[3], geo::CryostatID& cid) const
  {
    geo::CryostatID::CryostatID_t cstat = FindCryostatAtPosition(worldLoc);

    if(cstat == geo::CryostatID::InvalidID)
      throw cet::exception("GeometryCore") << "Can't find Cryostat for position ("
                                       << worldLoc[0] << ","
                                       << worldLoc[1] << ","
                                       << worldLoc[2] << ")\n";
    cid = geo::CryostatID(cstat);
    return Cryostat(cid);
  } // GeometryCore::PositionToCryostat(double[3], CryostatID)

  const CryostatGeo& GeometryCore::PositionToCryostat
    (double const worldLoc[3], unsigned int &cstat) const
  {
    geo::CryostatID cid;
    geo::CryostatGeo const& cryo = PositionToCryostat(worldLoc, cid);
    cstat = cid.Cryostat;
    return cryo;
  } // GeometryCore::PositionToCryostat(double[3], unsigned int)

  //......................................................................
  unsigned int GeometryCore::FindAuxDetAtPosition(geo::Point_t const& point) const
  {
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    std::array<double, 3U> worldPos = {{ point.X(), point.Y(), point.Z() }};
    return fChannelMapAlg->NearestAuxDet(worldPos.data(), AuxDets());
  } // GeometryCore::FindAuxDetAtPosition()

  //......................................................................
  unsigned int GeometryCore::FindAuxDetAtPosition(double const  worldPos[3]) const
    { return FindAuxDetAtPosition(geo::vect::makePointFromCoords(worldPos)); }


  //......................................................................
  const AuxDetGeo& GeometryCore::PositionToAuxDet
    (geo::Point_t const& point, unsigned int &ad) const
  {
    // locate the desired Auxiliary Detector
    ad = FindAuxDetAtPosition(point);
    return AuxDet(ad);
  }

  //......................................................................
  const AuxDetGeo& GeometryCore::PositionToAuxDet
    (double const  worldLoc[3], unsigned int &ad) const
    { return PositionToAuxDet(geo::vect::makePointFromCoords(worldLoc), ad); }

  //......................................................................
  void GeometryCore::FindAuxDetSensitiveAtPosition
    (geo::Point_t const& point, std::size_t& adg, std::size_t& sv) const
  {
    adg = FindAuxDetAtPosition(point);
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    std::array<double, 3U> const worldPos = {{ point.X(), point.Y(), point.Z() }};
    sv  = fChannelMapAlg->NearestSensitiveAuxDet(worldPos.data(), AuxDets());
  } // GeometryCore::FindAuxDetAtPosition()

  //......................................................................
  void GeometryCore::FindAuxDetSensitiveAtPosition
    (double const worldPos[3], size_t& adg, size_t& sv) const
    { return FindAuxDetSensitiveAtPosition(geo::vect::makePointFromCoords(worldPos), adg, sv); }


  //......................................................................
  const AuxDetSensitiveGeo& GeometryCore::PositionToAuxDetSensitive
    (geo::Point_t const& point, size_t& ad, size_t& sv) const
  {
    // locate the desired Auxiliary Detector
    FindAuxDetSensitiveAtPosition(point, ad, sv);
    return AuxDet(ad).SensitiveVolume(sv);
  }

  //......................................................................
  const AuxDetSensitiveGeo& GeometryCore::PositionToAuxDetSensitive
    (double const worldLoc[3], size_t& ad, size_t& sv) const
    { return PositionToAuxDetSensitive(geo::vect::makePointFromCoords(worldLoc), ad, sv); }

  //......................................................................
  const AuxDetGeo& GeometryCore::ChannelToAuxDet(std::string const& auxDetName,
                                             uint32_t    const& channel) const
  {
    size_t adIdx = fChannelMapAlg->ChannelToAuxDet(AuxDets(), auxDetName, channel);
    return this->AuxDet(adIdx);
  }

  //......................................................................
  const AuxDetSensitiveGeo& GeometryCore::ChannelToAuxDetSensitive(std::string const& auxDetName,
                                                      uint32_t    const& channel) const
  {
    auto idx = fChannelMapAlg->ChannelToSensitiveAuxDet(AuxDets(), auxDetName, channel);
    return this->AuxDet(idx.first).SensitiveVolume(idx.second);
  }


  //......................................................................
  SigType_t GeometryCore::SignalType(raw::ChannelID_t const channel) const
  {
    return fChannelMapAlg->SignalTypeForChannel(channel);
  }

  //......................................................................
  SigType_t GeometryCore::SignalType(geo::PlaneID const& pid) const
  {
    // map wire plane -> readout plane -> first channel,
    // then use SignalType(channel)

    auto const ropid = WirePlaneToROP(pid);
    if (!ropid.isValid) {
      throw cet::exception("GeometryCore")
        << "SignalType(): Mapping of wire plane " << std::string(pid)
        << " to readout plane failed!\n";
    }
    return SignalType(ropid);

  } // GeometryCore::SignalType(PlaneID)


  //......................................................................
  View_t GeometryCore::View(raw::ChannelID_t const channel) const {
    return (channel == raw::InvalidChannelID)
      ? geo::kUnknown: View(ChannelToROP(channel));
  } // GeometryCore::View()

  //......................................................................
  View_t GeometryCore::View(geo::PlaneID const& pid) const
  {
    return pid? Plane(pid).View(): geo::kUnknown;
  }

  //--------------------------------------------------------------------
  bool GeometryCore::HasChannel(raw::ChannelID_t channel) const {
    return fChannelMapAlg->HasChannel(channel);
  } // GeometryCore::HasChannel()


  //......................................................................
  std::set<PlaneID> const& GeometryCore::PlaneIDs() const
  {
    return fChannelMapAlg->PlaneIDs();
  }

  //......................................................................
  const std::string GeometryCore::GetWorldVolumeName() const
  {
    // For now, and possibly forever, this is a constant (given the
    // definition of "nodeNames" above).
    return std::string("volWorld");
  }


  //......................................................................
  geo::BoxBoundedGeo GeometryCore::DetectorEnclosureBox
    (std::string const& name /* = "volDetEnclosure" */) const
  {
    auto const& path = FindDetectorEnclosure(name);
    if (path.empty()) {
      throw cet::exception("GeometryCore")
        << "DetectorEnclosureBox(): can't find enclosure volume '" << name << "'\n";
    }

    TGeoVolume const* pEncl = path.back()->GetVolume();
    auto const* pBox = dynamic_cast<TGeoBBox const*>(pEncl->GetShape());

    // check that this is indeed a box
    if (!pBox) {
      // at initialisation time we don't know yet our real ID
      throw cet::exception("GeometryCore") << "Detector enclosure '"
        << name << "' is not a box! (it is a " << pEncl->GetShape()->IsA()->GetName()
        << ")\n";
    }

    geo::LocalTransformation<TGeoHMatrix> trans(path, path.size() - 1);
    // get the half width, height, etc of the cryostat
    const double halfwidth  = pBox->GetDX();
    const double halfheight = pBox->GetDY();
    const double halflength = pBox->GetDZ();

    return {
      trans.LocalToWorld(geo::Point_t{ -halfwidth, -halfheight, -halflength }),
      trans.LocalToWorld(geo::Point_t{ +halfwidth, +halfheight, +halflength })
      };
  } // geo::GeometryCore::DetectorEnclosureBox()

  //......................................................................
  struct NodeNameMatcherClass {
    std::set<std::string> const* vol_names;

    NodeNameMatcherClass(std::set<std::string> const& names)
      : vol_names(&names) {}

    /// Returns whether the specified node matches a set of names
    bool operator() (TGeoNode const& node) const
      {
        if (!vol_names) return true;
        return vol_names->find(node.GetVolume()->GetName()) != vol_names->end();
      }

  }; // NodeNameMatcherClass

  struct CollectNodesByName {
    std::vector<TGeoNode const*> nodes;

    CollectNodesByName(std::set<std::string> const& names): matcher(names) {}

    /// If the name of the node matches, records the end node
    void operator() (TGeoNode const& node)
      { if (matcher(node)) nodes.push_back(&node); }

    void operator() (ROOTGeoNodeForwardIterator const& iter)
      { operator() (**iter); }

      protected:
    NodeNameMatcherClass matcher;
  }; // CollectNodesByName

  struct CollectPathsByName {
    std::vector<std::vector<TGeoNode const*>> paths;

    CollectPathsByName(std::set<std::string> const& names): matcher(names) {}

    /// If the name of the node matches, records the node full path
    void operator() (ROOTGeoNodeForwardIterator const& iter)
      { if (matcher(**iter)) paths.push_back(iter.get_path()); }

      protected:
    NodeNameMatcherClass matcher;
  }; // CollectPathsByName

  std::vector<TGeoNode const*> GeometryCore::FindAllVolumes
    (std::set<std::string> const& vol_names) const
  {
    CollectNodesByName node_collector(vol_names);

    ROOTGeoNodeForwardIterator iNode(ROOTGeoManager()->GetTopNode());
    TGeoNode const* pCurrentNode;

    while ((pCurrentNode = *iNode)) {
      node_collector(*pCurrentNode);
      ++iNode;
    } // while

    return node_collector.nodes;
  } // GeometryCore::FindAllVolumes()


  std::vector<std::vector<TGeoNode const*>> GeometryCore::FindAllVolumePaths
    (std::set<std::string> const& vol_names) const
  {
    CollectPathsByName path_collector(vol_names);

    ROOTGeoNodeForwardIterator iNode(ROOTGeoManager()->GetTopNode());

    while (*iNode) {
      path_collector(iNode);
      ++iNode;
    } // while

    return path_collector.paths;
  } // GeometryCore::FindAllVolumePaths()



  //......................................................................
  std::string GeometryCore::GetLArTPCVolumeName(geo::TPCID const& tpcid) const
  {
    return std::string(TPC(tpcid).ActiveVolume()->GetName());
  }

  //......................................................................
  std::string GeometryCore::GetCryostatVolumeName(geo::CryostatID const& cid) const
  {
    return std::string(Cryostat(cid).Volume()->GetName());
  }

  //......................................................................
  geo::Length_t GeometryCore::DetHalfWidth(geo::TPCID const& tpcid)  const
  {
    return TPC(tpcid).ActiveHalfWidth();
  }

  //......................................................................
  geo::Length_t GeometryCore::DetHalfHeight(geo::TPCID const& tpcid) const
  {
    return TPC(tpcid).ActiveHalfHeight();
  }

  //......................................................................
  geo::Length_t GeometryCore::DetLength(geo::TPCID const& tpcid) const
  {
    return TPC(tpcid).ActiveLength();
  }

  //......................................................................
  geo::Length_t GeometryCore::CryostatHalfWidth
    (geo::CryostatID const& cid) const
  {
    return Cryostat(cid).HalfWidth();
  }

  //......................................................................
  geo::Length_t GeometryCore::CryostatHalfHeight
    (geo::CryostatID const& cid) const
  {
    return Cryostat(cid).HalfHeight();
  }

  //......................................................................
  geo::Length_t GeometryCore::CryostatLength(geo::CryostatID const& cid) const
  {
    return Cryostat(cid).Length();
  }

  //......................................................................
  void GeometryCore::CryostatBoundaries
    (double* boundaries, geo::CryostatID const& cid) const
  {
    geo::CryostatGeo const& cryo = Cryostat(cid);
    cryo.Boundaries(boundaries);
  } // GeometryCore::CryostatBoundaries()

  //......................................................................
  // This method returns the distance between the specified planes.
  // p1 < p2
  double GeometryCore::PlanePitch(
    geo::TPCID const& tpcid,
    geo::PlaneID::PlaneID_t p1, geo::PlaneID::PlaneID_t p2
    ) const
  {
    return TPC(tpcid).PlanePitch(p1, p2);
  }

  double GeometryCore::PlanePitch
    (geo::PlaneID const& pid1, geo::PlaneID const& pid2) const
  {
    return PlanePitch(pid1.asTPCID(), pid1.Plane, pid2.Plane);
  }

  double GeometryCore::PlanePitch(unsigned int p1,
                              unsigned int p2,
                              unsigned int tpc,
                              unsigned int cstat) const
  {
    return PlanePitch(geo::TPCID(cstat, tpc), p1, p2);
  }

  //......................................................................
  // This method returns the distance between wires in a plane.
  geo::Length_t GeometryCore::WirePitch(geo::PlaneID const& planeid) const
  {
    return Plane(planeid).WirePitch();
  }

  //......................................................................
  // This method returns the distance between wires in the specified view
  // it assumes all planes of a given view have the same pitch
  geo::Length_t GeometryCore::WirePitch(geo::View_t view) const
  {
    // look in cryostat 0, tpc 0 to find the plane with the
    // specified view
    return TPC({ 0, 0 }).Plane(view).WirePitch();
  }

  //......................................................................
  // This method returns the distance between wires in the specified view
  // it assumes all planes of a given view have the same pitch
  double GeometryCore::WireAngleToVertical
    (geo::View_t view, geo::TPCID const& tpcid) const
  {
    // loop over the planes in cryostat 0, tpc 0 to find the plane with the
    // specified view
    geo::TPCGeo const& TPC = this->TPC(tpcid);
    for (unsigned int p = 0; p < TPC.Nplanes(); ++p) {
      geo::PlaneGeo const& plane = TPC.Plane(p);
      if (plane.View() == view) return plane.ThetaZ();
    } // for
    throw cet::exception("GeometryCore") << "WireAngleToVertical(): no view \""
      << geo::PlaneGeo::ViewName(view) << "\" (#" << ((int) view)
      << ") in " << std::string(tpcid);
  } // GeometryCore::WireAngleToVertical()

  //......................................................................
  unsigned int GeometryCore::MaxTPCs() const {
    unsigned int maxTPCs = 0;
    for (geo::CryostatGeo const& cryo: Cryostats()) {
      unsigned int maxTPCsInCryo = cryo.NTPC();
      if (maxTPCsInCryo > maxTPCs) maxTPCs = maxTPCsInCryo;
    } // for
    return maxTPCs;
  } // GeometryCore::MaxTPCs()

  //......................................................................
  unsigned int GeometryCore::TotalNTPC() const {
    // it looks like C++11 lambdas have made STL algorithms easier to use,
    // but only so much:
    return std::accumulate(Cryostats().begin(), Cryostats().end(), 0U,
      [](unsigned int sum, geo::CryostatGeo const& cryo)
        { return sum + cryo.NTPC(); }
      );
  } // GeometryCore::TotalNTPC()

  //......................................................................
  unsigned int GeometryCore::MaxPlanes() const {
    unsigned int maxPlanes = 0;
    for (geo::CryostatGeo const& cryo: Cryostats()) {
      unsigned int maxPlanesInCryo = cryo.MaxPlanes();
      if (maxPlanesInCryo > maxPlanes) maxPlanes = maxPlanesInCryo;
    } // for
    return maxPlanes;
  } // GeometryCore::MaxPlanes()

  //......................................................................
  unsigned int GeometryCore::MaxWires() const {
    unsigned int maxWires = 0;
    for (geo::CryostatGeo const& cryo: Cryostats()) {
      unsigned int maxWiresInCryo = cryo.MaxWires();
      if (maxWiresInCryo > maxWires) maxWires = maxWiresInCryo;
    } // for
    return maxWires;
  } // GeometryCore::MaxWires()

  //......................................................................
  TGeoVolume const* GeometryCore::WorldVolume() const {
    return gGeoManager->FindVolumeFast(GetWorldVolumeName().c_str());
  }

  //......................................................................
  geo::BoxBoundedGeo GeometryCore::WorldBox() const {

    TGeoVolume const* world = WorldVolume();
    if(!world) {
      throw cet::exception("GeometryCore")
        << "no world volume '" << GetWorldVolumeName() << "'\n";
    }
    TGeoShape const* s = world->GetShape();
    if(!s) {
      throw cet::exception("GeometryCore")
        << "world volume '" << GetWorldVolumeName() << "' is shapeless!!!\n";
    }

    double x1, x2, y1, y2, z1, z2;
    s->GetAxisRange(1, x1, x2);
    s->GetAxisRange(2, y1, y2);
    s->GetAxisRange(3, z1, z2);

    // geo::BoxBoundedGeo constructor will sort the coordinates as needed
    return geo::BoxBoundedGeo{ x1, x2, y1, y2, z1, z2 };
  } // GeometryCore::WorldBox()

  //......................................................................
  void GeometryCore::WorldBox(double* xlo, double* xhi,
                          double* ylo, double* yhi,
                          double* zlo, double* zhi) const
  {
    geo::BoxBoundedGeo const box = WorldBox();
    if (xlo) *xlo = box.MinX();
    if (ylo) *ylo = box.MinY();
    if (zlo) *zlo = box.MinZ();
    if (xhi) *xhi = box.MaxX();
    if (yhi) *yhi = box.MaxY();
    if (zhi) *zhi = box.MaxZ();
  }

  //......................................................................
  std::string GeometryCore::VolumeName(geo::Point_t const& point) const
  {
    // check that the given point is in the World volume at least
    TGeoVolume const*volWorld = WorldVolume();
    double halflength = ((TGeoBBox*)volWorld->GetShape())->GetDZ();
    double halfheight = ((TGeoBBox*)volWorld->GetShape())->GetDY();
    double halfwidth  = ((TGeoBBox*)volWorld->GetShape())->GetDX();
    if(std::abs(point.x()) > halfwidth  ||
       std::abs(point.y()) > halfheight ||
       std::abs(point.z()) > halflength
       ){
      mf::LogWarning("GeometryCoreBadInputPoint") << "point (" << point.x() << ","
                                              << point.y() << "," << point.z() << ") "
                                              << "is not inside the world volume "
                                              << " half width = " << halfwidth
                                              << " half height = " << halfheight
                                              << " half length = " << halflength
                                              << " returning unknown volume name";
      const std::string unknown("unknownVolume");
      return unknown;
    }

    return gGeoManager->FindNode(point.X(), point.Y(), point.Z())->GetName();
  }

  //......................................................................
  TGeoMaterial const* GeometryCore::Material(geo::Point_t const& point) const {
    auto const pNode = gGeoManager->FindNode(point.X(), point.Y(), point.Z());
    if (!pNode) return nullptr;
    auto const pMedium = pNode->GetMedium();
    return pMedium? pMedium->GetMaterial(): nullptr;
  }

  //......................................................................
  std::string GeometryCore::MaterialName(geo::Point_t const& point) const
  {
    // check that the given point is in the World volume at least
    geo::BoxBoundedGeo worldBox = WorldBox();
    if (!worldBox.ContainsPosition(point)) {
      mf::LogWarning("GeometryCoreBadInputPoint")
        << "point " << point << " is not inside the world volume "
        << worldBox.Min() << " -- " << worldBox.Max()
        << "; returning unknown material name";
      return { "unknownMaterial" };
    }
    auto const pMaterial = Material(point);
    if (!pMaterial) {
      mf::LogWarning("GeometryCoreBadInputPoint")
        << "material for point " << point
        << " not found! returning unknown material name";
      return { "unknownMaterial" };
    }
    return pMaterial->GetName();
  }


  //......................................................................
  std::vector<TGeoNode const*> GeometryCore::FindDetectorEnclosure
    (std::string const& name /* = "volDetEnclosure" */) const
  {
     std::vector<TGeoNode const*> path { ROOTGeoManager()->GetTopNode() };
     if (!FindFirstVolume(name, path)) path.clear();
     return path;
  } // FindDetectorEnclosure()

  //......................................................................
  bool GeometryCore::FindFirstVolume
    (std::string const& name, std::vector<const TGeoNode*>& path) const
  {
    assert(!path.empty());

    auto const* pCurrent = path.back();

    // first check the current layer
    if (strncmp(name.c_str(), pCurrent->GetName(), name.length()) == 0)
      return true;

    //explore the next layer down
    auto const* pCurrentVolume = pCurrent->GetVolume();
    unsigned int nd = pCurrentVolume->GetNdaughters();
    for(unsigned int i = 0; i < nd; ++i) {
      path.push_back(pCurrentVolume->GetNode(i));
      if (FindFirstVolume(name, path)) return true;
      path.pop_back();
    } // for
    return false;
  } // FindFirstVolume()


  //......................................................................
  void GeometryCore::BuildGeometry(geo::GeometryBuilder& builder)
  {
    geo::GeoNodePath path{ gGeoManager->GetTopNode() };
    Cryostats() = builder.extractCryostats(path);
    AuxDets() = builder.extractAuxiliaryDetectors(path);
  }

  //......................................................................
  //
  // Return the total mass of the detector
  //
  //
  double GeometryCore::TotalMass(std::string vol) const
  {
    //the TGeoNode::GetVolume() returns the TGeoVolume of the detector outline
    //and ROOT calculates the mass in kg for you
    TGeoVolume *gvol = gGeoManager->FindVolumeFast(vol.c_str());
    if(gvol) return gvol->Weight();

    throw cet::exception("GeometryCore") << "could not find specified volume '"
                                         << vol
                                         << " 'to determine total mass\n";
  }

  //......................................................................
  double GeometryCore::MassBetweenPoints
    (geo::Point_t const& p1, geo::Point_t const& p2) const
  {

    //The purpose of this method is to determine the column density
    //between the two points given.  Do that by starting at p1 and
    //stepping until you get to the node of p2.  calculate the distance
    //between the point just inside that node and p2 to get the last
    //bit of column density
    double columnD = 0.;

    //first initialize a track - get the direction cosines
    geo::Vector_t const dir = (p2 - p1).Unit();

    double const dxyz[3] = { dir.X(), dir.Y(), dir.Z() };
    double const cp1[3] = { p1.X(), p1.Y(), p1.Z() };
    gGeoManager->InitTrack(cp1, dxyz);

    //might be helpful to have a point to a TGeoNode
    TGeoNode *node = gGeoManager->GetCurrentNode();

    //check that the points are not in the same volume already.
    //if they are in different volumes, keep stepping until you
    //are in the same volume as the second point
    while(!gGeoManager->IsSameLocation(p2.X(), p2.Y(), p2.Z())){
      gGeoManager->FindNextBoundary();
      columnD += gGeoManager->GetStep()*node->GetMedium()->GetMaterial()->GetDensity();

      //the act of stepping puts you in the next node and returns that node
      node = gGeoManager->Step();
    }//end loop to get to volume of second point

    //now you are in the same volume as the last point, but not at that point.
    //get the distance between the current point and the last one
    geo::Point_t const last
      = geo::vect::makePointFromCoords(gGeoManager->GetCurrentPoint());
    double const lastStep = (p2 - last).R();
    columnD += lastStep*node->GetMedium()->GetMaterial()->GetDensity();

    return columnD;
  }

  //......................................................................
  std::string GeometryCore::Info(std::string indent /* = "" */) const {
    std::ostringstream sstr;
    Print(sstr, indent);
    return sstr.str();
  } // GeometryCore::Info()



  //......................................................................
  std::vector< geo::WireID > GeometryCore::ChannelToWire( raw::ChannelID_t channel ) const
  {
    return fChannelMapAlg->ChannelToWire(channel);
  }

  //--------------------------------------------------------------------
  readout::ROPID GeometryCore::ChannelToROP(raw::ChannelID_t channel) const
  {
    return fChannelMapAlg->ChannelToROP(channel);
  } // GeometryCore::ChannelToROP()


  //----------------------------------------------------------------------------
  geo::Length_t GeometryCore::WireCoordinate
    (geo::Point_t const& pos, geo::PlaneID const& planeid) const
  {
    return Plane(planeid).WireCoordinate(pos);
  }

  //----------------------------------------------------------------------------
  geo::Length_t GeometryCore::WireCoordinate
    (double YPos, double ZPos, geo::PlaneID const& planeid) const
  {
    return fChannelMapAlg->WireCoordinate(YPos, ZPos, planeid);
  }

  //----------------------------------------------------------------------------
  // The NearestWire and PlaneWireToChannel are attempts to speed
  // up the simulation by memorizing the computationally intensive
  // setup steps for some geometry calculations.  The results are
  // valid assuming the wire planes are comprised of straight,
  // parallel wires with constant pitch across the entire plane, with
  // a hierarchical numbering scheme - Ben J Oct 2011
  unsigned int GeometryCore::NearestWire
    (geo::Point_t const& point, geo::PlaneID const& planeid) const
  {
    return NearestWireID(point, planeid).Wire;
    // return fChannelMapAlg->NearestWire(worldPos, planeid);
  }

  //----------------------------------------------------------------------------
  unsigned int GeometryCore::NearestWire
    (const double worldPos[3], geo::PlaneID const& planeid) const
  {
    return NearestWire(geo::vect::makePointFromCoords(worldPos), planeid);
  }

  //----------------------------------------------------------------------------
  unsigned int GeometryCore::NearestWire
    (std::vector<double> const& worldPos, geo::PlaneID const& planeid) const
  {
    if(worldPos.size() > 3) throw cet::exception("GeometryCore") << "bad size vector for "
                                                             << "worldPos: "
                                                             << worldPos.size() << "\n";
    return NearestWire(worldPos.data(), planeid);
  }

  //----------------------------------------------------------------------------
  geo::WireID GeometryCore::NearestWireID
    (geo::Point_t const& worldPos, geo::PlaneID const& planeid) const
  {
    return Plane(planeid).NearestWireID(worldPos);
  }

  //----------------------------------------------------------------------------
  geo::WireID GeometryCore::NearestWireID
    (std::vector<double> const& worldPos, geo::PlaneID const& planeid) const
  {
    if(worldPos.size() > 3) throw cet::exception("GeometryCore") << "bad size vector for "
                                                             << "worldPos: "
                                                             << worldPos.size() << "\n";
    return NearestWireID(worldPos.data(), planeid);
  }

  //----------------------------------------------------------------------------
  geo::WireID GeometryCore::NearestWireID
    (const double worldPos[3], geo::PlaneID const& planeid) const
  {
    return NearestWireID(geo::vect::makePointFromCoords(worldPos), planeid);
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t GeometryCore::NearestChannel
    (const double worldPos[3], geo::PlaneID const& planeid) const
  {
    return NearestChannel(geo::vect::makePointFromCoords(worldPos), planeid);
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t GeometryCore::NearestChannel
    (std::vector<double> const& worldPos, geo::PlaneID const& planeid) const
  {
    if(worldPos.size() > 3) throw cet::exception("GeometryCore") << "bad size vector for "
                                                             << "worldPos: "
                                                             << worldPos.size() << "\n";
    return NearestChannel(worldPos.data(), planeid);
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t GeometryCore::NearestChannel
    (geo::Point_t const& worldPos, geo::PlaneID const& planeid) const
  {

    // This method is supposed to return a channel number rather than
    //  a wire number.  Perform the conversion here (although, maybe
    //  faster if we deal in wire numbers rather than channel numbers?)

    // NOTE on failure both NearestChannel() and upstream:
    // * according to documentation, should return invalid channel
    // * in the actual code throw an exception because of a BUG
    //
    // The following implementation automatically becomes in fact compliant to
    // the documentation if upstreams becomes compliant to.
    // When that happens, just delete this comment.
    geo::WireID const wireID = NearestWireID(worldPos, planeid);
    return wireID? PlaneWireToChannel(wireID): raw::InvalidChannelID;
  } // GeometryCore::NearestChannel()

  //--------------------------------------
  raw::ChannelID_t GeometryCore::PlaneWireToChannel(WireID const& wireid) const
  {
    return fChannelMapAlg->PlaneWireToChannel(wireid);
  }

  // Functions to allow determination if two wires intersect, and if so where.
  // This is useful information during 3D reconstruction.
  //......................................................................
  bool GeometryCore::ValueInRange(double value, double min, double max) const
  {
    if(min>max) std::swap(min,max);//protect against funny business due to wire angles
    if (std::abs(value-min)<1e-6||std::abs(value-max)<1e-6) return true;
    return (value>=min) && (value<=max);
  }

  //......................................................................
  void GeometryCore::WireEndPoints
    (geo::WireID const& wireid, double *xyzStart, double *xyzEnd) const
  {
    Segment_t result = WireEndPoints(wireid);

    xyzStart[0] = result.start().X();
    xyzStart[1] = result.start().Y();
    xyzStart[2] = result.start().Z();
    xyzEnd[0]   = result.end().X();
    xyzEnd[1]   = result.end().Y();
    xyzEnd[2]   = result.end().Z();

    if(xyzEnd[2]<xyzStart[2]){
      //ensure that "End" has higher z-value than "Start"
      std::swap(xyzStart[0],xyzEnd[0]);
      std::swap(xyzStart[1],xyzEnd[1]);
      std::swap(xyzStart[2],xyzEnd[2]);
    }
    if(xyzEnd[1]<xyzStart[1] && std::abs(xyzEnd[2]-xyzStart[2])<0.01){
      // if wire is vertical ensure that "End" has higher y-value than "Start"
      std::swap(xyzStart[0],xyzEnd[0]);
      std::swap(xyzStart[1],xyzEnd[1]);
      std::swap(xyzStart[2],xyzEnd[2]);
    }

  } // GeometryCore::WireEndPoints(WireID, 2x double*)

  //Changed to use WireIDsIntersect(). Apr, 2015 T.Yang
  //......................................................................
  bool GeometryCore::ChannelsIntersect(raw::ChannelID_t c1,
                                   raw::ChannelID_t c2,
                                   double &y,
                                   double &z) const
  {

    // [GP] these errors should be exceptions, and this function is deprecated
    // because it violates interoperability
    std::vector<geo::WireID> chan1wires = ChannelToWire(c1);
    if (chan1wires.empty()) {
      mf::LogError("ChannelsIntersect")
        << "1st channel " << c1 << " maps to no wire (is it a real one?)";
      return false;
    }
    std::vector<geo::WireID> chan2wires = ChannelToWire(c2);
    if (chan2wires.empty()) {
      mf::LogError("ChannelsIntersect")
        << "2nd channel " << c2 << " maps to no wire (is it a real one?)";
      return false;
    }

    if (chan1wires.size() > 1) {
      mf::LogWarning("ChannelsIntersect")
        << "1st channel " << c1 << " maps to " << chan2wires.size()
        << " wires; using the first!";
      return false;
    }
    if (chan2wires.size() > 1) {
      mf::LogError("ChannelsIntersect")
        << "2nd channel " << c2 << " maps to " << chan2wires.size()
        << " wires; using the first!";
      return false;
    }

    geo::WireIDIntersection widIntersect;
    if (this->WireIDsIntersect(chan1wires[0],chan2wires[0],widIntersect)){
      y = widIntersect.y;
      z = widIntersect.z;
      return true;
    }
    else{
      y = widIntersect.y;
      z = widIntersect.z;
      return false;
    }
  }


  //......................................................................
  bool GeometryCore::IntersectLines(
    double A_start_x, double A_start_y, double A_end_x, double A_end_y,
    double B_start_x, double B_start_y, double B_end_x, double B_end_y,
    double& x, double& y
  ) const {

    // Equation from http://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    // T.Yang Nov, 2014
    // Notation: x => coordinate orthogonal to the drift direction and to the beam direction
    //           y => direction orthogonal to the previous and to beam direction

    double const denom = (A_start_x - A_end_x)*(B_start_y - B_end_y)
      - (A_start_y - A_end_y)*(B_start_x - B_end_x);

    if (coordIs.zero(denom)) return false;

    double const A = (A_start_x * A_end_y - A_start_y * A_end_x) / denom;
    double const B = (B_start_x * B_end_y - B_start_y * B_end_x) / denom;

    x = (B_start_x - B_end_x) * A - (A_start_x - A_end_x) * B;
    y = (B_start_y - B_end_y) * A - (A_start_y - A_end_y) * B;

    return true;

  } // GeometryCore::IntersectLines()

  //......................................................................
  bool GeometryCore::IntersectSegments(
    double A_start_x, double A_start_y, double A_end_x, double A_end_y,
    double B_start_x, double B_start_y, double B_end_x, double B_end_y,
    double& x, double& y
  ) const {

    bool bCross = IntersectLines(
      A_start_x, A_start_y, A_end_x, A_end_y,
      B_start_x, B_start_y, B_end_x, B_end_y,
      x, y
      );

    if (bCross) {
      mf::LogWarning("IntersectSegments") << "The segments are parallel!";
      return false;
    }

    return PointWithinSegments(
      A_start_x, A_start_y, A_end_x, A_end_y,
      B_start_x, B_start_y, B_end_x, B_end_y,
      x, y
      );

  } // GeometryCore::IntersectSegments()

  //......................................................................
  bool GeometryCore::WireIDsIntersect(
    const geo::WireID& wid1, const geo::WireID& wid2,
    geo::WireIDIntersection & widIntersect
  ) const {

    static_assert(
      std::numeric_limits<decltype(widIntersect.y)>::has_infinity,
      "the vector coordinate type can't represent infinity!"
      );
    constexpr auto infinity
      = std::numeric_limits<decltype(widIntersect.y)>::infinity();

    if (!WireIDIntersectionCheck(wid1, wid2)) {
      widIntersect.y = widIntersect.z = infinity;
      widIntersect.TPC = geo::TPCID::InvalidID;
      return false;
    }

    // get the endpoints to see if wires intersect
    Segment_t const w1 = WireEndPoints(wid1);
    Segment_t const w2 = WireEndPoints(wid2);

    // TODO extract the coordinates in the right way;
    // is it any worth, since then the result is in (y, z), whatever it means?
    bool const cross = IntersectLines(
      w1.start()[1], w1.start()[2], w1.end()[1], w1.end()[2],
      w2.start()[1], w2.start()[2], w2.end()[1], w2.end()[2],
      widIntersect.y, widIntersect.z
      );
    if (!cross) {
      widIntersect.y = widIntersect.z = infinity;
      widIntersect.TPC = geo::TPCID::InvalidID;
      return false;
    }
    bool const within = PointWithinSegments(
      w1.start()[1], w1.start()[2], w1.end()[1], w1.end()[2],
      w2.start()[1], w2.start()[2], w2.end()[1], w2.end()[2],
      widIntersect.y, widIntersect.z
      );

    widIntersect.TPC = (within? wid1.TPC: geo::TPCID::InvalidID);

    // return whether the intersection is within the length of both wires
    return within;

  } // GeometryCore::WireIDsIntersect(WireIDIntersection)


  //......................................................................
  bool GeometryCore::WireIDsIntersect(
    const geo::WireID& wid1, const geo::WireID& wid2,
    geo::Point_t& intersection
  )
    const
  {
    //
    // This is not a real 3D intersection: the wires do not cross, since they
    // are required to belong to two different planes.
    //
    // After Christopher Backhouse suggestion, we take the point on the first
    // wire which is closest to the other one.
    //
    //
    static_assert(
      std::numeric_limits<decltype(intersection.X())>::has_infinity,
      "the vector coordinate type can't represent infinity!"
      );
    constexpr auto infinity
      = std::numeric_limits<decltype(intersection.X())>::infinity();

    if (!WireIDIntersectionCheck(wid1, wid2)) {
      intersection = { infinity, infinity, infinity };
      return false;
    }

    /*
     * the point on the first wire:
     *
     *     p1(t) = c1 + t w1 (c1 the center of the wire, w1 its direction[1])
     *
     * has the minimal distance from the other wire (c2 + u w2, with the same
     * notation) at
     *
     *     t = [(dc,w1) - (dc,w2)(w1,w2)] / [ 1 - (w1,w2)^2 ]
     *     u = [-(dc,w2) + (dc,w1)(w1,w2)] / [ 1 - (w1,w2)^2 ]
     *
     * (where (a,b) is a scalar product and dc = (c2 - c1) ).
     * If w1 and w2 are unit vectors, t and u are in fact the distance of the
     * point from the center of the respective wires in "standard" geometry
     * units.
     *
     */

    geo::WireGeo const& wire1 = Wire(wid1);
    decltype(auto) c1 = wire1.GetCenter();
    decltype(auto) w1 = wire1.Direction();

    geo::WireGeo const& wire2 = Wire(wid2);
    decltype(auto) c2 = wire2.GetCenter();
    decltype(auto) w2 = wire2.Direction();

    auto const dc = c2 - c1;

    // note: we are not checking that w1 and w2 are not parallel.
    using geo::vect::dot;
    double const w1w2 = dot(w1, w2); // this is cos(angle), angle between wires
    double const cscAngle2 = 1.0 / (1.0 - cet::square(w1w2)); // this is 1/sin^2(angle)
    double const dcw1 = dot(dc, w1);
    double const dcw2 = dot(dc, w2);
    double const t = (dcw1 - (dcw2 * w1w2)) * cscAngle2;
    double const u = (-dcw2 + (dcw1 * w1w2)) * cscAngle2;

    intersection = c1 + t * w1;

    bool const within
      = (std::abs(t) <= wire1.HalfL()) && (std::abs(u) <= wire2.HalfL());

    // return whether the intersection is within the length of both wires
    return within;

  } // GeometryCore::WireIDsIntersect(Point3D_t)


  //......................................................................
  bool GeometryCore::WireIDsIntersect
    (const geo::WireID& wid1, const geo::WireID& wid2, TVector3& intersection)
    const
  {
    geo::Point_t p;
    bool res = WireIDsIntersect(wid1, wid2, p);
    intersection = geo::vect::toTVector3(p);
    return res;
  } // GeometryCore::WireIDsIntersect(TVector3)


  //----------------------------------------------------------------------------
  geo::PlaneID GeometryCore::ThirdPlane
    (geo::PlaneID const& pid1, geo::PlaneID const& pid2) const
  {
    // how many planes in the TPC pid1 belongs to:
    const unsigned int nPlanes = Nplanes(pid1);
    if(nPlanes != 3) {
      throw cet::exception("GeometryCore")
        << "ThirdPlane() supports only TPCs with 3 planes, and I see "
        << nPlanes << " instead\n";
    }

    geo::PlaneID::PlaneID_t target_plane = nPlanes;
    for (geo::PlaneID::PlaneID_t iPlane = 0; iPlane < nPlanes; ++iPlane){
      if ((iPlane == pid1.Plane) || (iPlane == pid2.Plane)) continue;
      if (target_plane != nPlanes) {
        throw cet::exception("GeometryCore")
          << "ThirdPlane() found too many planes that are not "
          << std::string(pid1) << " nor " << std::string(pid2)
          << "! (first " << target_plane << ", then " << iPlane << ")\n";
      } // if we had a target already
      target_plane = iPlane;
    } // for
    if (target_plane == nPlanes) {
      throw cet::exception("GeometryCore")
        << "ThirdPlane() can't find a plane that is not " << std::string(pid1)
        << " nor " << std::string(pid2) << "!\n";
    }

    return geo::PlaneID(pid1, target_plane);
  } // GeometryCore::ThirdPlane()


  void GeometryCore::CheckIndependentPlanesOnSameTPC
    (geo::PlaneID const& pid1, geo::PlaneID const& pid2, const char* caller)
  {
    if(pid1.asTPCID() != pid2.asTPCID()) {
      throw cet::exception("GeometryCore")
        << caller << " needs two planes on the same TPC (got "
        << std::string(pid1) << " and " << std::string(pid2) << ")\n";
    }
    if(pid1 == pid2) { // was: return 999;
      throw cet::exception("GeometryCore")
        << caller << " needs two different planes, got "
        << std::string(pid1) << " twice\n";
    }
  } // GeometryCore::CheckIndependentPlanesOnSameTPC()


  double GeometryCore::ThirdPlaneSlope(
    geo::PlaneID const& pid1, double slope1,
    geo::PlaneID const& pid2, double slope2,
    geo::PlaneID const& output_plane
  ) const {

    CheckIndependentPlanesOnSameTPC(pid1, pid2, "ThirdPlaneSlope()");

    geo::TPCGeo const& TPC = this->TPC(pid1);

    // We need the "wire coordinate direction" for each plane.
    // This is perpendicular to the wire orientation.
    // PlaneGeo::PhiZ() defines the right orientation too.
    return ComputeThirdPlaneSlope(
      TPC.Plane(pid1).PhiZ(), slope1,
      TPC.Plane(pid2).PhiZ(), slope2,
      TPC.Plane(output_plane).PhiZ()
      );
  } // ThirdPlaneSlope()


  double GeometryCore::ThirdPlaneSlope(
    geo::PlaneID const& pid1, double slope1,
    geo::PlaneID const& pid2, double slope2
  ) const {
    geo::PlaneID target_plane = ThirdPlane(pid1, pid2);
    return ThirdPlaneSlope(pid1, slope1, pid2, slope2, target_plane);
  } // ThirdPlaneSlope()


  double GeometryCore::ThirdPlane_dTdW(
    geo::PlaneID const& pid1, double slope1,
    geo::PlaneID const& pid2, double slope2,
    geo::PlaneID const& output_plane
  ) const {

    CheckIndependentPlanesOnSameTPC(pid1, pid2, "ThirdPlane_dTdW()");

    geo::TPCGeo const& TPC = this->TPC(pid1);

    double angle[3], pitch[3];
    geo::PlaneGeo const* const planes[3]
      = { &TPC.Plane(pid1), &TPC.Plane(pid2), &TPC.Plane(output_plane) };

    // We need wire pitch and "wire coordinate direction" for each plane.
    // The latter is perpendicular to the wire orientation.
    // PlaneGeo::PhiZ() defines the right orientation too.
    for (size_t i = 0; i < 3; ++i) {
      angle[i] = planes[i]->PhiZ();
      pitch[i] = planes[i]->WirePitch();
    } // for

    return ComputeThirdPlane_dTdW(
      angle[0], pitch[0], slope1,
      angle[1], pitch[1], slope2,
      angle[2], pitch[2]
      );

  } // GeometryCore::ThirdPlane_dTdW()


  double GeometryCore::ThirdPlane_dTdW(
    geo::PlaneID const& pid1, double slope1,
    geo::PlaneID const& pid2, double slope2
  ) const {
    geo::PlaneID target_plane = ThirdPlane(pid1, pid2);
    return ThirdPlane_dTdW(pid1, slope1, pid2, slope2, target_plane);
  } // ThirdPlane_dTdW()


  // Given slopes dTime/dWire in two planes, return with the slope in the 3rd plane.
  // Requires slopes to be in the same metrics,
  // e.g. converted in a distances ratio.
  // B. Baller August 2014
  // Rewritten by T. Yang Apr 2015 using the equation in H. Greenlee's talk:
  // https://cdcvs.fnal.gov/redmine/attachments/download/1821/larsoft_apr20_2011.pdf
  // slide 2
  double GeometryCore::ComputeThirdPlaneSlope
    (double angle1, double slope1, double angle2, double slope2, double angle3)
  {
    // note that, if needed, the trigonometric functions can be pre-calculated.

    // Can't resolve very small slopes
    if ((std::abs(slope1) < 0.001) && (std::abs(slope2)) < 0.001) return 0.001;

    // We need the "wire coordinate direction" for each plane.
    // This is perpendicular to the wire orientation.
    double slope3 = 0.001;
    if (std::abs(slope1) > 0.001 && std::abs(slope2) > 0.001) {
      slope3
        = (
          + (1./slope1)*std::sin(angle3-angle2)
          - (1./slope2)*std::sin(angle3-angle1)
        ) / std::sin(angle1-angle2)
        ;
    }
    if (slope3 != 0.) slope3 = 1./slope3;
    else slope3 = 999.;

    return slope3;
  } // GeometryCore::ComputeThirdPlaneSlope()


  double GeometryCore::ComputeThirdPlane_dTdW(
    double angle1, double pitch1, double dTdW1,
    double angle2, double pitch2, double dTdW2,
    double angle_target, double pitch_target
      )
  {
    // we need to convert dt/dw into homogeneous coordinates, and then back;
    // slope = [dT * (TDCperiod / driftVelocity)] / [dW * wirePitch]
    // The coefficient of dT is assumed to be the same for all the planes,
    // and it finally cancels out. Pitches cancel out only if they are all
    // the same.
    return pitch_target * ComputeThirdPlaneSlope
      (angle1, dTdW1 / pitch1, angle2, dTdW2 / pitch2, angle_target);
  } // GeometryCore::ComputeThirdPlane_dTdW()


  //......................................................................
  // This function is called if it is determined that two wires in a single TPC must overlap.
  // To determine the yz coordinate of the wire intersection, we need to know the
  // endpoints of both wires in xyz-space, and also their orientation (angle), and the
  // inner dimensions of the TPC frame.
  // Note: This calculation is entirely dependent  on an accurate GDML description of the TPC!
  // Mitch - Feb., 2011
  // Changed to use WireIDsIntersect(). It does not check whether the intersection is on both wires (the same as the old behavior). T. Yang - Apr, 2015
  //--------------------------------------------------------------------
  bool GeometryCore::IntersectionPoint(geo::WireID const& wid1,
                                   geo::WireID const& wid2,
                                   double &y, double &z) const
  {
    geo::WireIDIntersection widIntersect;
    bool const found = WireIDsIntersect(wid1, wid2, widIntersect);
    y = widIntersect.y;
    z = widIntersect.z;
    return found;
  }

  //============================================================================
  //===  TPC set information
  //===
  //--------------------------------------------------------------------
  unsigned int GeometryCore::NTPCsets(readout::CryostatID const& cryoid) const {
    return fChannelMapAlg->NTPCsets(cryoid);
  } // GeometryCore::NTPCsets()


  //--------------------------------------------------------------------
  unsigned int GeometryCore::MaxTPCsets() const {
    return fChannelMapAlg->MaxTPCsets();
  } // GeometryCore::MaxTPCsets()


  //--------------------------------------------------------------------
  bool GeometryCore::HasTPCset(readout::TPCsetID const& tpcsetid) const {
    return fChannelMapAlg->HasTPCset(tpcsetid);
  } // GeometryCore::HasTPCset()


  //--------------------------------------------------------------------
  readout::TPCsetID GeometryCore::FindTPCsetAtPosition
    (double const worldLoc[3]) const
  {
    return TPCtoTPCset(FindTPCAtPosition(worldLoc));
  } // GeometryCore::FindTPCsetAtPosition()


  //--------------------------------------------------------------------
  readout::TPCsetID GeometryCore::TPCtoTPCset(geo::TPCID const& tpcid) const
  {
    return fChannelMapAlg->TPCtoTPCset(tpcid);
  } // GeometryCore::TPCtoTPCset()


  //--------------------------------------------------------------------
  std::vector<geo::TPCID> GeometryCore::TPCsetToTPCs
    (readout::TPCsetID const& tpcsetid) const
  {
    return fChannelMapAlg->TPCsetToTPCs(tpcsetid);
  } // GeometryCore::TPCsetToTPCs()


  //============================================================================
  //===  Readout plane information
  //===
  //--------------------------------------------------------------------
  unsigned int GeometryCore::NROPs(readout::TPCsetID const& tpcsetid) const {
    return fChannelMapAlg->NROPs(tpcsetid);
  } // GeometryCore::NROPs()


  //--------------------------------------------------------------------
  unsigned int GeometryCore::MaxROPs() const {
    return fChannelMapAlg->MaxROPs();
  } // GeometryCore::MaxROPs()


  //--------------------------------------------------------------------
  bool GeometryCore::HasROP(readout::ROPID const& ropid) const {
    return fChannelMapAlg->HasROP(ropid);
  } // GeometryCore::HasROP()


  //--------------------------------------------------------------------
  readout::ROPID GeometryCore::WirePlaneToROP(geo::PlaneID const& planeid) const
  {
    return fChannelMapAlg->WirePlaneToROP(planeid);
  } // GeometryCore::WirePlaneToROP()


  //--------------------------------------------------------------------
  std::vector<geo::PlaneID> GeometryCore::ROPtoWirePlanes
    (readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->ROPtoWirePlanes(ropid);
  } // GeometryCore::ROPtoWirePlanes()


  //--------------------------------------------------------------------
  std::vector<geo::TPCID> GeometryCore::ROPtoTPCs
    (readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->ROPtoTPCs(ropid);
  } // GeometryCore::ROPtoTPCs()


  //--------------------------------------------------------------------
  raw::ChannelID_t GeometryCore::FirstChannelInROP
    (readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->FirstChannelInROP(ropid);
  } // GeometryCore::FirstChannelInROP()


  //--------------------------------------------------------------------
  geo::View_t GeometryCore::View(readout::ROPID const& ropid) const {
    return View(fChannelMapAlg->FirstWirePlaneInROP(ropid));
  } // GeometryCore::View()


  //--------------------------------------------------------------------
  geo::SigType_t GeometryCore::SignalType(readout::ROPID const& ropid) const {
    return fChannelMapAlg->SignalTypeForROPID(ropid);
  } // GeometryCore::SignalType(ROPID)




  //============================================================================
  //--------------------------------------------------------------------
  // Convert OpDet, Cryo into unique OpDet number
  unsigned int GeometryCore::OpDetFromCryo(unsigned int o, unsigned int c ) const
  {
    static bool Loaded=false;
    static std::vector<unsigned int> LowestID;
    static unsigned int NCryo;
    // If not yet loaded static parameters, do it
    if(Loaded == false){

      Loaded = true;

      // Store the lowest ID for each cryostat
      NCryo=Ncryostats();
      LowestID.resize(NCryo + 1);
      LowestID.at(0)=0;
      for(size_t cryo=0; cryo!=NCryo; ++cryo){
        LowestID.at(cryo+1)=LowestID.at(cryo)+Cryostat(c).NOpDet();
      }

    }

    if( (c < NCryo) && (o < Cryostat(c).NOpDet())){
      return LowestID.at(c)+o;
    }
    else{
      throw cet::exception("OpDetCryoToOpID Error") << "Coordinates c=" << c
                                                    << ", o=" << o
                                                    << " out of range. Abort\n";
    }

    // if all is well, we never get to this point in the method
    // but still a good idea to be sure to always return something.

    return INT_MAX;
  }

  
  //--------------------------------------------------------------------
  // Return gdml string which gives sensitive opdet name
  std::string GeometryCore::OpDetGeoName(unsigned int c) const
  {
    return Cryostat(c).OpDetGeoName();
  }

  //--------------------------------------------------------------------
  geo::OpDetGeo const& GeometryCore::OpDet(geo::OpDetID const& oid) const {
    
    if (!oid.isValid) {
      throw cet::exception("GeometryCore")
        << "geo::GeometryCore::OpDet(): invalid optical detector specified.\n";
    }
    
    return Cryostat(oid).OpDet(oid);
    
  } // GeometryCore::OpDet()
  
  
  //--------------------------------------------------------------------
  const OpDetGeo& GeometryCore::OpDetGeoFromOpChannel(unsigned int OpChannel) const
  {
    return this->OpDetGeoFromOpDet(this->OpDetFromOpChannel(OpChannel));
  }

  //--------------------------------------------------------------------
  const OpDetGeo& GeometryCore::OpDetGeoFromOpDet(unsigned int OpDet) const
  {
    static bool Loaded=false;
    static std::vector<unsigned int> LowestID;
    static size_t NCryo;
    // If not yet loaded static parameters, do it
    if(Loaded == false){

      Loaded = true;

      // Store the lowest ID for each cryostat
      NCryo=Ncryostats();
      LowestID.resize(NCryo + 1);
      LowestID[0] = 0;
      for(size_t cryo = 0; cryo != NCryo; ++cryo){
        LowestID[cryo+1] = LowestID[cryo] + Cryostat(cryo).NOpDet();
      }

    }

    for(size_t i=0; i!=NCryo; ++i){
      if( (OpDet >= LowestID[i]) && (OpDet < LowestID[i+1]) ){
        int c = i;
        int o = OpDet-LowestID[i];
        return this->Cryostat(c).OpDet(o);
      }
    }
    // If we made it here, we didn't find the right combination. abort
    throw cet::exception("OpID To OpDetCryo error")<<"OpID out of range, "<< OpDet << "\n";

    // Will not reach due to exception
    return this->Cryostat(0).OpDet(0);
  }


  //--------------------------------------------------------------------
  // Find the closest OpChannel to this point, in the appropriate cryostat
  unsigned int GeometryCore::GetClosestOpDet(geo::Point_t const& point) const
  {
    geo::CryostatGeo const* cryo = PositionToCryostatPtr(point);
    if (!cryo) return std::numeric_limits<unsigned int>::max();
    int o = cryo->GetClosestOpDet(point);
    return OpDetFromCryo(o, cryo->ID().Cryostat);
  }


  //--------------------------------------------------------------------
  // Find the closest OpChannel to this point, in the appropriate cryostat
  unsigned int GeometryCore::GetClosestOpDet(double const* point) const
    { return GetClosestOpDet(geo::vect::makePointFromCoords(point)); }


  //--------------------------------------------------------------------
  geo::TPCGeo const& GeometryCore::GetTPCclosestToOpDet
    (geo::OpDetID const& oid) const
  {
    if (!oid.isValid) {
      throw cet::exception("GeometryCore")
        << "geo::GeometryCore::GetTPCclosestToOpDet(): "
        "invalid optical detector specified.\n";
    }
    
    return Cryostat(oid).GetTPCclosestToOpDet(oid);
    
  } // GeometryCore::GetTPCclosestToOpDet()
  
  //--------------------------------------------------------------------
  bool GeometryCore::WireIDIntersectionCheck
    (const geo::WireID& wid1, const geo::WireID& wid2) const
  {
    if (wid1.asTPCID() != wid2) {
      mf::LogError("WireIDIntersectionCheck")
        << "Comparing two wires on different TPCs: return failure.";
      return false;
    }
    if (wid1.Plane == wid2.Plane) {
      mf::LogError("WireIDIntersectionCheck")
        << "Comparing two wires in the same plane: return failure";
      return false;
    }
    if (!HasWire(wid1)) {
      mf::LogError("WireIDIntersectionCheck")
        << "1st wire " << wid1 << " does not exist (max wire number: "
        << Nwires(wid1.planeID()) << ")";
      return false;
    }
    if (!HasWire(wid2)) {
      mf::LogError("WireIDIntersectionCheck")
        << "2nd wire " << wid2 << " does not exist (max wire number: "
        << Nwires(wid2.planeID()) << ")";
      return false;
    }
    return true;
  } // GeometryCore::WireIDIntersectionCheck()


  //--------------------------------------------------------------------
  bool GeometryCore::PointWithinSegments(
    double A_start_x, double A_start_y, double A_end_x, double A_end_y,
    double B_start_x, double B_start_y, double B_end_x, double B_end_y,
    double x, double y
  ) {
    return coordIs.withinSorted(x, A_start_x, A_end_x)
      && coordIs.withinSorted(y, A_start_y, A_end_y)
      && coordIs.withinSorted(x, B_start_x, B_end_x)
      && coordIs.withinSorted(y, B_start_y, B_end_y)
      ;

  } // GeometryCore::PointWithinSegments()

  //--------------------------------------------------------------------
  constexpr details::geometry_iterator_types::BeginPos_t
    details::geometry_iterator_types::begin_pos;
  constexpr details::geometry_iterator_types::EndPos_t
    details::geometry_iterator_types::end_pos;
  constexpr details::geometry_iterator_types::UndefinedPos_t
    details::geometry_iterator_types::undefined_pos;

  //--------------------------------------------------------------------
  //--- ROOTGeoNodeForwardIterator
  //---

  ROOTGeoNodeForwardIterator& ROOTGeoNodeForwardIterator::operator++ () {
    if (current_path.empty()) return *this;
    if (current_path.size() == 1) { current_path.pop_back(); return *this; }

    // I am done; all my descendants were also done already;
    // first look at my younger siblings
    NodeInfo_t& current = current_path.back();
    NodeInfo_t const& parent = current_path[current_path.size() - 2];
    if (++(current.sibling) < parent.self->GetNdaughters()) {
      // my next sibling exists, let's parse his descendents
      current.self = parent.self->GetDaughter(current.sibling);
      reach_deepest_descendant();
    }
    else current_path.pop_back(); // no sibling, it's time for mum
    return *this;
  } // ROOTGeoNodeForwardIterator::operator++


  //--------------------------------------------------------------------
  std::vector<TGeoNode const*> ROOTGeoNodeForwardIterator::get_path() const {

    std::vector<TGeoNode const*> node_path(current_path.size());

    std::transform(current_path.begin(), current_path.end(), node_path.begin(),
      [](NodeInfo_t const& node_info){ return node_info.self; });
    return node_path;

  } // ROOTGeoNodeForwardIterator::path()


  //--------------------------------------------------------------------
  void ROOTGeoNodeForwardIterator::reach_deepest_descendant() {
    Node_t descendent = current_path.back().self;
    while (descendent->GetNdaughters() > 0) {
      descendent = descendent->GetDaughter(0);
      current_path.emplace_back(descendent, 0U);
    } // while
  } // ROOTGeoNodeForwardIterator::reach_deepest_descendant()

  //--------------------------------------------------------------------
  void ROOTGeoNodeForwardIterator::init(TGeoNode const* start_node) {
    current_path.clear();
    if (!start_node) return;
    current_path.emplace_back(start_node, 0U);
    reach_deepest_descendant();
  } // ROOTGeoNodeForwardIterator::init()

  //--------------------------------------------------------------------

} // namespace geo
