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
#include "larcorealg/CoreUtils/NumericUtils.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/Decomposer.h" // geo::vect::dot()
#include "larcorealg/Geometry/GeometryBuilderStandard.h"
#include "larcorealg/Geometry/Intersections.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"                // geo::vect
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h"        // geo::vect
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi<>

// Framework includes
#include "cetlib/pow.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>

// C/C++ includes
#include <algorithm> // std::for_each(), std::transform()
#include <cctype>    // ::tolower()
#include <cmath>     // std::abs() ...
#include <cstddef>   // size_t
#include <iterator>  // std::back_inserter()
#include <limits>    // std::numeric_limits<>
#include <numeric>   // std::accumulate
#include <sstream>   // std::ostringstream
#include <tuple>
#include <utility> // std::swap()
#include <vector>

namespace {
  /// Throws an exception ("GeometryCore" category) unless pid1 and pid2
  /// are on different planes of the same TPC (ID validity is not checked)
  void CheckIndependentPlanesOnSameTPC(geo::PlaneID const& pid1,
                                       geo::PlaneID const& pid2,
                                       const char* caller)
  {
    if (pid1.asTPCID() != pid2.asTPCID()) {
      throw cet::exception("GeometryCore")
        << caller << " needs two planes on the same TPC (got " << std::string(pid1) << " and "
        << std::string(pid2) << ")\n";
    }
    if (pid1 == pid2) { // was: return 999;
      throw cet::exception("GeometryCore")
        << caller << " needs two different planes, got " << std::string(pid1) << " twice\n";
    }
  }
}

namespace geo {

  //......................................................................
  // Constructor.
  GeometryCore::GeometryCore(fhicl::ParameterSet const& pset)
    : fSurfaceY(pset.get<double>("SurfaceY"))
    , fDetectorName(pset.get<std::string>("Name"))
    , fMinWireZDist(pset.get<double>("MinWireZDist", 3.0))
    , fPositionWiggle(pset.get<double>("PositionEpsilon", 1.e-4))
    , fBuilderParameters(pset.get<fhicl::ParameterSet>("Builder", fhicl::ParameterSet()))
  {
    std::transform(fDetectorName.begin(), fDetectorName.end(), fDetectorName.begin(), ::tolower);
  }

  //......................................................................
  void GeometryCore::ApplyChannelMap(std::unique_ptr<ChannelMapAlg> pChannelMap)
  {
    SortGeometry(pChannelMap->Sorter());
    UpdateAfterSorting(); // after channel mapping has sorted objects, set their IDs
    pChannelMap->Initialize(fGeoData);
    fChannelMapAlg = move(pChannelMap);
  }

  //......................................................................
  void GeometryCore::LoadGeometryFile(std::string gdmlfile,
                                      std::string rootfile,
                                      GeometryBuilder& builder,
                                      bool bForceReload /* = false*/
  )
  {
    if (gdmlfile.empty()) {
      throw cet::exception("GeometryCore") << "No GDML Geometry file specified!\n";
    }

    if (rootfile.empty()) {
      throw cet::exception("GeometryCore") << "No ROOT Geometry file specified!\n";
    }

    ClearGeometry();

    // Open the GDML file, and convert it into ROOT TGeoManager format.
    // Then lock the gGeoManager to prevent future imports, for example
    // in AuxDetGeometry
    if (!gGeoManager || bForceReload) {
      if (gGeoManager)
        TGeoManager::UnlockGeometry();
      else { // very first time (or so it should)
        // [20210630, petrillo@slac.stanford.edu]
        // ROOT 6.22.08 allows us to choose the representation of lengths
        // in the geometry objects parsed from GDML.
        // In LArSoft we want them to be centimeters (ROOT standard).
        // This was tracked as Redmine issue #25990, and I leave this mark
        // because I feel that we'll be back to it not too far in the future.
        // Despite the documentation (ROOT 6.22/08),
        // it seems the units are locked from the beginning,
        // so we unlock without prejudice.
        TGeoManager::LockDefaultUnits(false);
        TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);
        TGeoManager::LockDefaultUnits(true);
      }
      TGeoManager::Import(rootfile.c_str());
      gGeoManager->LockGeometry();
    }

    BuildGeometry(builder);

    fGDMLfile = move(gdmlfile);
    fROOTfile = move(rootfile);

    mf::LogInfo("GeometryCore") << "New detector geometry loaded from "
                                << "\n\t" << fROOTfile << "\n\t" << fGDMLfile << "\n";
  } // GeometryCore::LoadGeometryFile()

  //......................................................................
  void GeometryCore::LoadGeometryFile(std::string gdmlfile,
                                      std::string rootfile,
                                      bool bForceReload /* = false*/
  )
  {
    fhicl::Table<GeometryBuilderStandard::Config> const builderConfig(fBuilderParameters,
                                                                      {"tool_type"});
    // this is a wink to the understanding that we might be using an art-based
    // service provider configuration sprinkled with tools.
    GeometryBuilderStandard builder{builderConfig()};
    LoadGeometryFile(gdmlfile, rootfile, builder, bForceReload);
  }

  //......................................................................
  void GeometryCore::ClearGeometry()
  {
    fGeoData = {};
  }

  //......................................................................
  void GeometryCore::SortGeometry(GeoObjectSorter const& sorter)
  {
    mf::LogInfo("GeometryCore") << "Sorting volumes...";

    sorter.SortAuxDets(AuxDets());
    sorter.SortCryostats(Cryostats());

    CryostatID::CryostatID_t c = 0;
    for (CryostatGeo& cryo : Cryostats()) {
      cryo.SortSubVolumes(sorter);
      cryo.UpdateAfterSorting(CryostatID(c));
      ++c;
    }
  }

  //......................................................................
  void GeometryCore::UpdateAfterSorting()
  {
    for (size_t c = 0; c < Ncryostats(); ++c)
      Cryostats()[c].UpdateAfterSorting(CryostatID(c));

    allViews.clear();
    for (auto const& tpc : Iterate<TPCGeo>()) {
      auto const& TPCviews = tpc.Views();
      allViews.insert(TPCviews.cbegin(), TPCviews.cend());
    }
  }

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
  }

  //......................................................................
  std::vector<raw::ChannelID_t> GeometryCore::ChannelsInTPCs() const
  {
    std::vector<raw::ChannelID_t> channels;
    channels.reserve(fChannelMapAlg->Nchannels());

    for (auto const& ts : Iterate<readout::TPCsetID>()) {
      for (auto const t : fChannelMapAlg->TPCsetToTPCs(ts)) {
        for (auto const& wire : Iterate<WireID>(t)) {
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
    int N = 0;
    for (size_t cstat = 0; cstat != Ncryostats(); ++cstat)
      N += this->Cryostat(CryostatID(cstat)).NOpDet();
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
  bool GeometryCore::IsValidOpChannel(int opChannel) const
  {
    return fChannelMapAlg->IsValidOpChannel(opChannel, this->NOpDets());
  }

  //......................................................................
  unsigned int GeometryCore::NAuxDetSensitive(size_t const& aid) const
  {
    if (aid < NAuxDets()) { return AuxDets()[aid].NSensitiveVolume(); }
    throw cet::exception("Geometry")
      << "Requested AuxDet index " << aid << " is out of range: " << NAuxDets();
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
  CryostatGeo const& GeometryCore::Cryostat(CryostatID const& cryoid) const
  {
    if (auto pCryo = CryostatPtr(cryoid)) { return *pCryo; }
    throw cet::exception("GeometryCore") << "Cryostat #" << cryoid.Cryostat << " does not exist\n";
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
  const AuxDetGeo& GeometryCore::AuxDet(unsigned int const ad) const
  {
    if (ad >= NAuxDets())
      throw cet::exception("GeometryCore") << "AuxDet " << ad << " does not exist\n";
    return AuxDets()[ad];
  }

  //......................................................................
  TPCID GeometryCore::FindTPCAtPosition(Point_t const& point) const
  {
    // first find the cryostat
    CryostatGeo const* cryo = PositionToCryostatPtr(point);
    if (!cryo) return {};

    // then ask it about the TPC
    TPCID tpcid = cryo->PositionToTPCID(point, 1. + fPositionWiggle);
    if (tpcid) return tpcid;

    // return an invalid TPC ID with cryostat information set:
    tpcid.Cryostat = cryo->ID().Cryostat;
    tpcid.markInvalid();
    return tpcid;
  }

  //......................................................................
  CryostatGeo const* GeometryCore::PositionToCryostatPtr(Point_t const& point) const
  {
    for (auto const& cryostat : Iterate<CryostatGeo>()) {
      if (cryostat.ContainsPosition(point, 1.0 + fPositionWiggle)) return &cryostat;
    }
    return nullptr;
  }

  //......................................................................
  CryostatID GeometryCore::PositionToCryostatID(Point_t const& point) const
  {
    CryostatGeo const* cryo = PositionToCryostatPtr(point);
    return cryo ? cryo->ID() : CryostatID{};
  }

  //......................................................................
  TPCGeo const* GeometryCore::PositionToTPCptr(Point_t const& point) const
  {
    CryostatGeo const* cryo = PositionToCryostatPtr(point);
    return cryo ? cryo->PositionToTPCptr(point, 1. + fPositionWiggle) : nullptr;
  }

  //......................................................................
  TPCGeo const& GeometryCore::PositionToTPC(Point_t const& point) const
  {
    if (auto tpc = PositionToTPCptr(point)) { return *tpc; }
    throw cet::exception("GeometryCore") << "Can't find any TPC at position " << point << "\n";
  }

  //......................................................................
  TPCID GeometryCore::PositionToTPCID(Point_t const& point) const
  {
    TPCGeo const* tpc = PositionToTPCptr(point);
    return tpc ? tpc->ID() : TPCID{};
  }

  //......................................................................
  void GeometryCore::GetEndID(TPCID& id) const
  {
    if (MaxTPCs() == 0) {
      GetBeginID(id);
      id.markInvalid();
    }
    else {
      GetEndID(id.asCryostatID());
      id.deepestIndex() = 0;
    }
  }

  //......................................................................
  TPCID GeometryCore::GetEndTPCID(CryostatID const& id) const
  {
    if (CryostatGeo const* cryo = CryostatPtr(id); cryo && cryo->NTPC() > 0)
      return {id.Cryostat + 1, 0};
    TPCID tpcid = GetBeginTPCID(id);
    tpcid.markInvalid();
    return tpcid;
  }

  //......................................................................
  CryostatGeo const& GeometryCore::PositionToCryostat(Point_t const& point) const
  {
    if (auto cstat = PositionToCryostatPtr(point)) { return *cstat; }
    throw cet::exception("GeometryCore") << "Can't find any cryostat at position " << point << "\n";
  }

  //......................................................................
  unsigned int GeometryCore::FindAuxDetAtPosition(Point_t const& point, double tolerance) const
  {
    return fChannelMapAlg->NearestAuxDet(point, AuxDets(), tolerance);
  }

  //......................................................................
  const AuxDetGeo& GeometryCore::PositionToAuxDet(Point_t const& point,
                                                  unsigned int& ad,
                                                  double tolerance) const
  {
    // locate the desired Auxiliary Detector
    ad = FindAuxDetAtPosition(point, tolerance);
    return AuxDet(ad);
  }

  //......................................................................
  void GeometryCore::FindAuxDetSensitiveAtPosition(Point_t const& point,
                                                   std::size_t& adg,
                                                   std::size_t& sv,
                                                   double tolerance) const
  {
    adg = FindAuxDetAtPosition(point, tolerance);
    sv = fChannelMapAlg->NearestSensitiveAuxDet(point, AuxDets(), tolerance);
  }

  //......................................................................
  const AuxDetSensitiveGeo& GeometryCore::PositionToAuxDetSensitive(Point_t const& point,
                                                                    size_t& ad,
                                                                    size_t& sv,
                                                                    double tolerance) const
  {
    // locate the desired Auxiliary Detector
    FindAuxDetSensitiveAtPosition(point, ad, sv, tolerance);
    return AuxDet(ad).SensitiveVolume(sv);
  }

  //......................................................................
  const AuxDetGeo& GeometryCore::ChannelToAuxDet(std::string const& auxDetName,
                                                 uint32_t const& channel) const
  {
    size_t adIdx = fChannelMapAlg->ChannelToAuxDet(AuxDets(), auxDetName, channel);
    return this->AuxDet(adIdx);
  }

  //......................................................................
  const AuxDetSensitiveGeo& GeometryCore::ChannelToAuxDetSensitive(std::string const& auxDetName,
                                                                   uint32_t const& channel) const
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
  SigType_t GeometryCore::SignalType(PlaneID const& pid) const
  {
    // map wire plane -> readout plane -> first channel,
    // then use SignalType(channel)

    auto const ropid = WirePlaneToROP(pid);
    if (!ropid.isValid) {
      throw cet::exception("GeometryCore") << "SignalType(): Mapping of wire plane "
                                           << std::string(pid) << " to readout plane failed!\n";
    }
    return SignalType(ropid);
  }

  //......................................................................
  View_t GeometryCore::View(raw::ChannelID_t const channel) const
  {
    return (channel == raw::InvalidChannelID) ? kUnknown : View(ChannelToROP(channel));
  }

  //......................................................................
  View_t GeometryCore::View(PlaneID const& pid) const
  {
    return pid ? Plane(pid).View() : kUnknown;
  }

  //--------------------------------------------------------------------
  bool GeometryCore::HasChannel(raw::ChannelID_t channel) const
  {
    return fChannelMapAlg->HasChannel(channel);
  }

  //......................................................................
  const std::string GeometryCore::GetWorldVolumeName() const
  {
    // For now, and possibly forever, this is a constant (given the
    // definition of "nodeNames" above).
    return std::string("volWorld");
  }

  //......................................................................
  BoxBoundedGeo GeometryCore::DetectorEnclosureBox(
    std::string const& name /* = "volDetEnclosure" */) const
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
      throw cet::exception("GeometryCore")
        << "Detector enclosure '" << name << "' is not a box! (it is a "
        << pEncl->GetShape()->IsA()->GetName() << ")\n";
    }

    LocalTransformation<TGeoHMatrix> trans(path, path.size() - 1);
    // get the half width, height, etc of the cryostat
    const double halfwidth = pBox->GetDX();
    const double halfheight = pBox->GetDY();
    const double halflength = pBox->GetDZ();

    return {trans.LocalToWorld(Point_t{-halfwidth, -halfheight, -halflength}),
            trans.LocalToWorld(Point_t{+halfwidth, +halfheight, +halflength})};
  }

  //......................................................................
  /** **************************************************************************
   * @brief Iterator to navigate through all the nodes
   *
   * Note that this is not a fully standard forward iterator in that it lacks
   * of the postfix operator. The reason is that it's too expensive and it
   * should be avoided.
   * Also I did not bother declaring the standard type definitions
   * (that's just laziness).
   *
   * An example of iteration:
   *
   *     TGeoNode const* pCurrentNode;
   *
   *     ROOTGeoNodeForwardIterator iNode(geom->ROOTGeoManager()->GetTopNode());
   *     while ((pCurrentNode = *iNode)) {
   *       // do something with pCurrentNode
   *       ++iNode;
   *     } // while
   *
   * These iterators are one use only, and they can't be reset after a loop
   * is completed.
   */
  class ROOTGeoNodeForwardIterator {
  public:
    explicit ROOTGeoNodeForwardIterator(TGeoNode const* start_node);

    /// Returns the pointer to the current node, or nullptr if none
    TGeoNode const* operator*() const
    {
      return current_path.empty() ? nullptr : current_path.back().self;
    }

    /// Points to the next node, or to nullptr if there are no more
    ROOTGeoNodeForwardIterator& operator++();

    /// Returns the full path of the current node
    std::vector<TGeoNode const*> get_path() const;

  private:
    struct NodeInfo_t {
      TGeoNode const* self;
      int sibling;
    };

    /// which node, which sibling?
    std::vector<NodeInfo_t> current_path;

    void reach_deepest_descendant();
  }; // class ROOTGeoNodeForwardIterator

  struct NodeNameMatcherClass {
    std::set<std::string> const* vol_names;

    NodeNameMatcherClass(std::set<std::string> const& names) : vol_names(&names) {}

    /// Returns whether the specified node matches a set of names
    bool operator()(TGeoNode const& node) const
    {
      if (!vol_names) return true;
      return vol_names->find(node.GetVolume()->GetName()) != vol_names->end();
    }

  }; // NodeNameMatcherClass

  struct CollectNodesByName {
    std::vector<TGeoNode const*> nodes;

    CollectNodesByName(std::set<std::string> const& names) : matcher(names) {}

    /// If the name of the node matches, records the end node
    void operator()(TGeoNode const& node)
    {
      if (matcher(node)) nodes.push_back(&node);
    }

    void operator()(ROOTGeoNodeForwardIterator const& iter) { operator()(**iter); }

  protected:
    NodeNameMatcherClass matcher;
  }; // CollectNodesByName

  struct CollectPathsByName {
    std::vector<std::vector<TGeoNode const*>> paths;

    CollectPathsByName(std::set<std::string> const& names) : matcher(names) {}

    /// If the name of the node matches, records the node full path
    void operator()(ROOTGeoNodeForwardIterator const& iter)
    {
      if (matcher(**iter)) paths.push_back(iter.get_path());
    }

  protected:
    NodeNameMatcherClass matcher;
  }; // CollectPathsByName

  //......................................................................
  std::vector<TGeoNode const*> GeometryCore::FindAllVolumes(
    std::set<std::string> const& vol_names) const
  {
    CollectNodesByName node_collector(vol_names);

    ROOTGeoNodeForwardIterator iNode{ROOTGeoManager()->GetTopNode()};
    TGeoNode const* pCurrentNode;

    while ((pCurrentNode = *iNode)) {
      node_collector(*pCurrentNode);
      ++iNode;
    }

    return node_collector.nodes;
  }

  //......................................................................
  std::vector<std::vector<TGeoNode const*>> GeometryCore::FindAllVolumePaths(
    std::set<std::string> const& vol_names) const
  {
    CollectPathsByName path_collector(vol_names);

    ROOTGeoNodeForwardIterator iNode(ROOTGeoManager()->GetTopNode());

    while (*iNode) {
      path_collector(iNode);
      ++iNode;
    }

    return path_collector.paths;
  }

  //......................................................................
  std::string GeometryCore::GetLArTPCVolumeName(TPCID const& tpcid) const
  {
    return TPC(tpcid).ActiveVolume()->GetName();
  }

  //......................................................................
  std::string GeometryCore::GetCryostatVolumeName(CryostatID const& cid) const
  {
    return Cryostat(cid).Volume()->GetName();
  }

  //......................................................................
  Length_t GeometryCore::DetHalfWidth(TPCID const& tpcid) const
  {
    return TPC(tpcid).ActiveHalfWidth();
  }

  //......................................................................
  Length_t GeometryCore::DetHalfHeight(TPCID const& tpcid) const
  {
    return TPC(tpcid).ActiveHalfHeight();
  }

  //......................................................................
  Length_t GeometryCore::DetLength(TPCID const& tpcid) const
  {
    return TPC(tpcid).ActiveLength();
  }

  //......................................................................
  Length_t GeometryCore::CryostatHalfWidth(CryostatID const& cid) const
  {
    return Cryostat(cid).HalfWidth();
  }

  //......................................................................
  Length_t GeometryCore::CryostatHalfHeight(CryostatID const& cid) const
  {
    return Cryostat(cid).HalfHeight();
  }

  //......................................................................
  Length_t GeometryCore::CryostatLength(CryostatID const& cid) const
  {
    return Cryostat(cid).Length();
  }

  //......................................................................
  void GeometryCore::GetEndID(PlaneID& id) const
  {
    if (MaxPlanes() == 0) {
      GetBeginID(id);
      id.markInvalid();
    }
    else {
      GetEndID(id.asTPCID());
      id.deepestIndex() = 0;
    }
  }

  //......................................................................
  PlaneID GeometryCore::GetEndPlaneID(CryostatID const& id) const
  {
    CryostatGeo const* cryo = CryostatPtr(id);
    return (cryo && cryo->MaxPlanes() > 0) ? PlaneID{GetEndTPCID(id), 0} : GetBeginPlaneID(id);
  }

  //......................................................................
  PlaneID GeometryCore::GetEndPlaneID(TPCID const& id) const
  {
    if (TPCGeo const* TPC = TPCPtr(id); TPC && TPC->Nplanes() > 0) return {GetNextID(id), 0};
    PlaneID pid = GetBeginPlaneID(id);
    pid.markInvalid();
    return pid;
  }

  //......................................................................
  // This method returns the distance between the specified planes.
  // p1 < p2
  double GeometryCore::PlanePitch(TPCID const& tpcid,
                                  PlaneID::PlaneID_t p1,
                                  PlaneID::PlaneID_t p2) const
  {
    return TPC(tpcid).PlanePitch(p1, p2);
  }

  double GeometryCore::PlanePitch(PlaneID const& pid1, PlaneID const& pid2) const
  {
    return PlanePitch(pid1.asTPCID(), pid1.Plane, pid2.Plane);
  }

  //......................................................................
  // This method returns the distance between wires in a plane.
  Length_t GeometryCore::WirePitch(PlaneID const& planeid) const
  {
    return Plane(planeid).WirePitch();
  }

  //......................................................................
  // This method returns the distance between wires in the specified view
  // it assumes all planes of a given view have the same pitch
  Length_t GeometryCore::WirePitch(View_t view) const
  {
    // look in cryostat 0, tpc 0 to find the plane with the
    // specified view
    return TPC({0, 0}).Plane(view).WirePitch();
  }

  //......................................................................
  // This method returns the distance between wires in the specified view
  // it assumes all planes of a given view have the same pitch
  double GeometryCore::WireAngleToVertical(View_t view, TPCID const& tpcid) const
  {
    // loop over the planes in cryostat 0, tpc 0 to find the plane with the
    // specified view
    TPCGeo const& TPC = this->TPC(tpcid);
    for (unsigned int p = 0; p < TPC.Nplanes(); ++p) {
      PlaneGeo const& plane = TPC.Plane(p);
      if (plane.View() == view) return plane.ThetaZ();
    } // for
    throw cet::exception("GeometryCore")
      << "WireAngleToVertical(): no view \"" << PlaneGeo::ViewName(view) << "\" (#" << ((int)view)
      << ") in " << std::string(tpcid);
  }

  //......................................................................
  unsigned int GeometryCore::MaxTPCs() const
  {
    unsigned int maxTPCs = 0;
    for (CryostatGeo const& cryo : Cryostats()) {
      unsigned int maxTPCsInCryo = cryo.NTPC();
      if (maxTPCsInCryo > maxTPCs) maxTPCs = maxTPCsInCryo;
    }
    return maxTPCs;
  }

  //......................................................................
  unsigned int GeometryCore::TotalNTPC() const
  {
    // it looks like C++11 lambdas have made STL algorithms easier to use,
    // but only so much:
    return std::accumulate(
      Cryostats().begin(), Cryostats().end(), 0U, [](unsigned int sum, CryostatGeo const& cryo) {
        return sum + cryo.NTPC();
      });
  }

  //......................................................................
  unsigned int GeometryCore::MaxPlanes() const
  {
    unsigned int maxPlanes = 0;
    for (CryostatGeo const& cryo : Cryostats()) {
      unsigned int maxPlanesInCryo = cryo.MaxPlanes();
      if (maxPlanesInCryo > maxPlanes) maxPlanes = maxPlanesInCryo;
    }
    return maxPlanes;
  }

  //......................................................................
  unsigned int GeometryCore::MaxWires() const
  {
    unsigned int maxWires = 0;
    for (CryostatGeo const& cryo : Cryostats()) {
      unsigned int maxWiresInCryo = cryo.MaxWires();
      if (maxWiresInCryo > maxWires) maxWires = maxWiresInCryo;
    }
    return maxWires;
  }

  //......................................................................
  void GeometryCore::GetEndID(WireID& id) const
  {
    if (MaxWires() == 0) {
      GetBeginID(id);
      id.markInvalid();
    }
    else {
      GetEndID(id.asPlaneID());
      id.deepestIndex() = 0;
    }
  }

  //......................................................................
  WireID GeometryCore::GetEndWireID(CryostatID const& id) const
  {
    CryostatGeo const* cryo = CryostatPtr(id);
    if (cryo && cryo->MaxWires() > 0) return {GetEndPlaneID(id), 0};
    WireID wid = GetBeginWireID(id);
    wid.markInvalid();
    return wid;
  }

  //......................................................................
  WireID GeometryCore::GetEndWireID(TPCID const& id) const
  {
    TPCGeo const* TPC = TPCPtr(id);
    if (TPC && TPC->MaxWires() > 0) return {GetEndPlaneID(id), 0};
    WireID wid = GetBeginWireID(id);
    wid.markInvalid();
    return wid;
  }

  //......................................................................
  WireID GeometryCore::GetEndWireID(PlaneID const& id) const
  {
    if (PlaneGeo const* plane = PlanePtr(id); plane && plane->Nwires() > 0)
      return {GetNextID(id), 0};
    WireID wid = GetBeginWireID(id);
    wid.markInvalid();
    return wid;
  }

  //......................................................................
  TGeoVolume const* GeometryCore::WorldVolume() const
  {
    return gGeoManager->FindVolumeFast(GetWorldVolumeName().c_str());
  }

  //......................................................................
  BoxBoundedGeo GeometryCore::WorldBox() const
  {
    TGeoVolume const* world = WorldVolume();
    if (!world) {
      throw cet::exception("GeometryCore") << "no world volume '" << GetWorldVolumeName() << "'\n";
    }
    TGeoShape const* s = world->GetShape();
    if (!s) {
      throw cet::exception("GeometryCore")
        << "world volume '" << GetWorldVolumeName() << "' is shapeless!!!\n";
    }

    double x1, x2, y1, y2, z1, z2;
    s->GetAxisRange(1, x1, x2);
    s->GetAxisRange(2, y1, y2);
    s->GetAxisRange(3, z1, z2);

    // BoxBoundedGeo constructor will sort the coordinates as needed
    return BoxBoundedGeo{x1, x2, y1, y2, z1, z2};
  }

  //......................................................................
  void GeometryCore::WorldBox(double* xlo,
                              double* xhi,
                              double* ylo,
                              double* yhi,
                              double* zlo,
                              double* zhi) const
  {
    BoxBoundedGeo const box = WorldBox();
    if (xlo) *xlo = box.MinX();
    if (ylo) *ylo = box.MinY();
    if (zlo) *zlo = box.MinZ();
    if (xhi) *xhi = box.MaxX();
    if (yhi) *yhi = box.MaxY();
    if (zhi) *zhi = box.MaxZ();
  }

  //......................................................................
  std::string GeometryCore::VolumeName(Point_t const& point) const
  {
    // check that the given point is in the World volume at least
    TGeoVolume const* volWorld = WorldVolume();
    double halflength = ((TGeoBBox*)volWorld->GetShape())->GetDZ();
    double halfheight = ((TGeoBBox*)volWorld->GetShape())->GetDY();
    double halfwidth = ((TGeoBBox*)volWorld->GetShape())->GetDX();
    if (std::abs(point.x()) > halfwidth || std::abs(point.y()) > halfheight ||
        std::abs(point.z()) > halflength) {
      mf::LogWarning("GeometryCoreBadInputPoint")
        << "point (" << point.x() << "," << point.y() << "," << point.z() << ") "
        << "is not inside the world volume "
        << " half width = " << halfwidth << " half height = " << halfheight
        << " half length = " << halflength << " returning unknown volume name";
      return "unknownVolume";
    }

    return gGeoManager->FindNode(point.X(), point.Y(), point.Z())->GetName();
  }

  //......................................................................
  TGeoMaterial const* GeometryCore::Material(Point_t const& point) const
  {
    auto const pNode = gGeoManager->FindNode(point.X(), point.Y(), point.Z());
    if (!pNode) return nullptr;
    auto const pMedium = pNode->GetMedium();
    return pMedium ? pMedium->GetMaterial() : nullptr;
  }

  //......................................................................
  std::string GeometryCore::MaterialName(Point_t const& point) const
  {
    // check that the given point is in the World volume at least
    BoxBoundedGeo worldBox = WorldBox();
    if (!worldBox.ContainsPosition(point)) {
      mf::LogWarning("GeometryCoreBadInputPoint")
        << "point " << point << " is not inside the world volume " << worldBox.Min() << " -- "
        << worldBox.Max() << "; returning unknown material name";
      return {"unknownMaterial"};
    }
    auto const pMaterial = Material(point);
    if (!pMaterial) {
      mf::LogWarning("GeometryCoreBadInputPoint")
        << "material for point " << point << " not found! returning unknown material name";
      return {"unknownMaterial"};
    }
    return pMaterial->GetName();
  }

  //......................................................................
  std::vector<TGeoNode const*> GeometryCore::FindDetectorEnclosure(
    std::string const& name /* = "volDetEnclosure" */) const
  {
    std::vector<TGeoNode const*> path{ROOTGeoManager()->GetTopNode()};
    if (!FindFirstVolume(name, path)) path.clear();
    return path;
  }

  //......................................................................
  bool GeometryCore::FindFirstVolume(std::string const& name,
                                     std::vector<const TGeoNode*>& path) const
  {
    assert(!path.empty());

    auto const* pCurrent = path.back();

    // first check the current layer
    if (strncmp(name.c_str(), pCurrent->GetName(), name.length()) == 0) return true;

    //explore the next layer down
    auto const* pCurrentVolume = pCurrent->GetVolume();
    unsigned int nd = pCurrentVolume->GetNdaughters();
    for (unsigned int i = 0; i < nd; ++i) {
      path.push_back(pCurrentVolume->GetNode(i));
      if (FindFirstVolume(name, path)) return true;
      path.pop_back();
    }
    return false;
  }

  //......................................................................
  void GeometryCore::BuildGeometry(GeometryBuilder& builder)
  {
    GeoNodePath path{gGeoManager->GetTopNode()};
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
    TGeoVolume* gvol = gGeoManager->FindVolumeFast(vol.c_str());
    if (gvol) return gvol->Weight();

    throw cet::exception("GeometryCore")
      << "could not find specified volume '" << vol << " 'to determine total mass\n";
  }

  //......................................................................
  double GeometryCore::MassBetweenPoints(Point_t const& p1, Point_t const& p2) const
  {
    //The purpose of this method is to determine the column density
    //between the two points given.  Do that by starting at p1 and
    //stepping until you get to the node of p2.  calculate the distance
    //between the point just inside that node and p2 to get the last
    //bit of column density
    double columnD = 0.;

    //first initialize a track - get the direction cosines
    Vector_t const dir = (p2 - p1).Unit();

    double const dxyz[3] = {dir.X(), dir.Y(), dir.Z()};
    double const cp1[3] = {p1.X(), p1.Y(), p1.Z()};
    gGeoManager->InitTrack(cp1, dxyz);

    //might be helpful to have a point to a TGeoNode
    TGeoNode* node = gGeoManager->GetCurrentNode();

    //check that the points are not in the same volume already.
    //if they are in different volumes, keep stepping until you
    //are in the same volume as the second point
    while (!gGeoManager->IsSameLocation(p2.X(), p2.Y(), p2.Z())) {
      gGeoManager->FindNextBoundary();
      columnD += gGeoManager->GetStep() * node->GetMedium()->GetMaterial()->GetDensity();

      //the act of stepping puts you in the next node and returns that node
      node = gGeoManager->Step();
    } //end loop to get to volume of second point

    //now you are in the same volume as the last point, but not at that point.
    //get the distance between the current point and the last one
    Point_t const last = vect::makePointFromCoords(gGeoManager->GetCurrentPoint());
    double const lastStep = (p2 - last).R();
    columnD += lastStep * node->GetMedium()->GetMaterial()->GetDensity();

    return columnD;
  }

  //......................................................................
  std::string GeometryCore::Info(std::string indent /* = "" */) const
  {
    std::ostringstream sstr;
    Print(sstr, indent);
    return sstr.str();
  }

  //......................................................................
  std::vector<WireID> GeometryCore::ChannelToWire(raw::ChannelID_t channel) const
  {
    return fChannelMapAlg->ChannelToWire(channel);
  }

  //--------------------------------------------------------------------
  readout::ROPID GeometryCore::ChannelToROP(raw::ChannelID_t channel) const
  {
    return fChannelMapAlg->ChannelToROP(channel);
  }

  //----------------------------------------------------------------------------
  Length_t GeometryCore::WireCoordinate(Point_t const& pos, PlaneID const& planeid) const
  {
    return Plane(planeid).WireCoordinate(pos);
  }

  //----------------------------------------------------------------------------
  WireID GeometryCore::NearestWireID(Point_t const& worldPos, PlaneID const& planeid) const
  {
    return Plane(planeid).NearestWireID(worldPos);
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t GeometryCore::NearestChannel(Point_t const& worldPos,
                                                PlaneID const& planeid) const
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
    WireID const wireID = NearestWireID(worldPos, planeid);
    return wireID ? PlaneWireToChannel(wireID) : raw::InvalidChannelID;
  }

  //--------------------------------------
  raw::ChannelID_t GeometryCore::PlaneWireToChannel(WireID const& wireid) const
  {
    return fChannelMapAlg->PlaneWireToChannel(wireid);
  }

  //......................................................................
  void GeometryCore::WireEndPoints(WireID const& wireid, double* xyzStart, double* xyzEnd) const
  {
    Segment_t result = WireEndPoints(wireid);

    xyzStart[0] = result.start().X();
    xyzStart[1] = result.start().Y();
    xyzStart[2] = result.start().Z();
    xyzEnd[0] = result.end().X();
    xyzEnd[1] = result.end().Y();
    xyzEnd[2] = result.end().Z();

    if (xyzEnd[2] < xyzStart[2]) {
      //ensure that "End" has higher z-value than "Start"
      std::swap(xyzStart[0], xyzEnd[0]);
      std::swap(xyzStart[1], xyzEnd[1]);
      std::swap(xyzStart[2], xyzEnd[2]);
    }
    if (xyzEnd[1] < xyzStart[1] && std::abs(xyzEnd[2] - xyzStart[2]) < 0.01) {
      // if wire is vertical ensure that "End" has higher y-value than "Start"
      std::swap(xyzStart[0], xyzEnd[0]);
      std::swap(xyzStart[1], xyzEnd[1]);
      std::swap(xyzStart[2], xyzEnd[2]);
    }
  }

  //Changed to use WireIDsIntersect(). Apr, 2015 T.Yang
  //......................................................................
  bool GeometryCore::ChannelsIntersect(raw::ChannelID_t c1,
                                       raw::ChannelID_t c2,
                                       double& y,
                                       double& z) const
  {
    // [GP] these errors should be exceptions, and this function is deprecated
    // because it violates interoperability
    std::vector<WireID> chan1wires = ChannelToWire(c1);
    if (chan1wires.empty()) {
      mf::LogError("ChannelsIntersect")
        << "1st channel " << c1 << " maps to no wire (is it a real one?)";
      return false;
    }
    std::vector<WireID> chan2wires = ChannelToWire(c2);
    if (chan2wires.empty()) {
      mf::LogError("ChannelsIntersect")
        << "2nd channel " << c2 << " maps to no wire (is it a real one?)";
      return false;
    }

    if (chan1wires.size() > 1) {
      mf::LogWarning("ChannelsIntersect")
        << "1st channel " << c1 << " maps to " << chan2wires.size() << " wires; using the first!";
      return false;
    }
    if (chan2wires.size() > 1) {
      mf::LogError("ChannelsIntersect")
        << "2nd channel " << c2 << " maps to " << chan2wires.size() << " wires; using the first!";
      return false;
    }

    WireIDIntersection widIntersect;
    if (this->WireIDsIntersect(chan1wires[0], chan2wires[0], widIntersect)) {
      y = widIntersect.y;
      z = widIntersect.z;
      return true;
    }
    else {
      y = widIntersect.y;
      z = widIntersect.z;
      return false;
    }
  }

  //......................................................................
  bool GeometryCore::WireIDsIntersect(const WireID& wid1,
                                      const WireID& wid2,
                                      WireIDIntersection& widIntersect) const
  {
    static_assert(std::numeric_limits<decltype(widIntersect.y)>::has_infinity,
                  "the vector coordinate type can't represent infinity!");
    constexpr auto infinity = std::numeric_limits<decltype(widIntersect.y)>::infinity();

    if (!WireIDIntersectionCheck(wid1, wid2)) {
      widIntersect.y = widIntersect.z = infinity;
      widIntersect.TPC = TPCID::InvalidID;
      return false;
    }

    // get the endpoints to see if wires intersect
    Segment_t const w1 = WireEndPoints(wid1);
    Segment_t const w2 = WireEndPoints(wid2);

    // TODO extract the coordinates in the right way;
    // is it any worth, since then the result is in (y, z), whatever it means?
    bool const cross = IntersectLines(w1.start().Y(),
                                      w1.start().Z(),
                                      w1.end().Y(),
                                      w1.end().Z(),
                                      w2.start().Y(),
                                      w2.start().Z(),
                                      w2.end().Y(),
                                      w2.end().Z(),
                                      widIntersect.y,
                                      widIntersect.z);
    if (!cross) {
      widIntersect.y = widIntersect.z = infinity;
      widIntersect.TPC = TPCID::InvalidID;
      return false;
    }
    bool const within = lar::util::PointWithinSegments(w1.start().Y(),
                                                       w1.start().Z(),
                                                       w1.end().Y(),
                                                       w1.end().Z(),
                                                       w2.start().Y(),
                                                       w2.start().Z(),
                                                       w2.end().Y(),
                                                       w2.end().Z(),
                                                       widIntersect.y,
                                                       widIntersect.z);

    widIntersect.TPC = (within ? wid1.TPC : TPCID::InvalidID);

    // return whether the intersection is within the length of both wires
    return within;
  }

  //......................................................................
  bool GeometryCore::WireIDsIntersect(const WireID& wid1,
                                      const WireID& wid2,
                                      Point_t& intersection) const
  {
    //
    // This is not a real 3D intersection: the wires do not cross, since they
    // are required to belong to two different planes.
    //
    // After Christopher Backhouse suggestion, we take the point on the first
    // wire which is closest to the other one.
    //
    //
    static_assert(std::numeric_limits<decltype(intersection.X())>::has_infinity,
                  "the vector coordinate type can't represent infinity!");
    constexpr auto infinity = std::numeric_limits<decltype(intersection.X())>::infinity();

    if (!WireIDIntersectionCheck(wid1, wid2)) {
      intersection = {infinity, infinity, infinity};
      return false;
    }

    WireGeo const& wire1 = Wire(wid1);
    WireGeo const& wire2 = Wire(wid2);

    // distance of the intersection point from the center of the two wires:
    IntersectionPointAndOffsets<Point_t> intersectionAndOffset =
      WiresIntersectionAndOffsets(wire1, wire2);
    intersection = intersectionAndOffset.point;

    bool const within = ((std::abs(intersectionAndOffset.offset1) <= wire1.HalfL()) &&
                         (std::abs(intersectionAndOffset.offset2) <= wire2.HalfL()));

    // return whether the intersection is within the length of both wires
    return within;
  }

  //----------------------------------------------------------------------------
  PlaneID GeometryCore::ThirdPlane(PlaneID const& pid1, PlaneID const& pid2) const
  {
    // how many planes in the TPC pid1 belongs to:
    const unsigned int nPlanes = Nplanes(pid1);
    if (nPlanes != 3) {
      throw cet::exception("GeometryCore")
        << "ThirdPlane() supports only TPCs with 3 planes, and I see " << nPlanes << " instead\n";
    }

    PlaneID::PlaneID_t target_plane = nPlanes;
    for (PlaneID::PlaneID_t iPlane = 0; iPlane < nPlanes; ++iPlane) {
      if ((iPlane == pid1.Plane) || (iPlane == pid2.Plane)) continue;
      if (target_plane != nPlanes) {
        throw cet::exception("GeometryCore")
          << "ThirdPlane() found too many planes that are not " << std::string(pid1) << " nor "
          << std::string(pid2) << "! (first " << target_plane << ", then " << iPlane << ")\n";
      } // if we had a target already
      target_plane = iPlane;
    } // for
    if (target_plane == nPlanes) {
      throw cet::exception("GeometryCore")
        << "ThirdPlane() can't find a plane that is not " << std::string(pid1) << " nor "
        << std::string(pid2) << "!\n";
    }

    return PlaneID(pid1, target_plane);
  }

  //----------------------------------------------------------------------------
  double GeometryCore::ThirdPlaneSlope(PlaneID const& pid1,
                                       double slope1,
                                       PlaneID const& pid2,
                                       double slope2,
                                       PlaneID const& output_plane) const
  {
    CheckIndependentPlanesOnSameTPC(pid1, pid2, "ThirdPlaneSlope()");

    TPCGeo const& TPC = this->TPC(pid1);

    // We need the "wire coordinate direction" for each plane.
    // This is perpendicular to the wire orientation.
    // PlaneGeo::PhiZ() defines the right orientation too.
    return ComputeThirdPlaneSlope(TPC.Plane(pid1).PhiZ(),
                                  slope1,
                                  TPC.Plane(pid2).PhiZ(),
                                  slope2,
                                  TPC.Plane(output_plane).PhiZ());
  }

  //----------------------------------------------------------------------------
  double GeometryCore::ThirdPlaneSlope(PlaneID const& pid1,
                                       double slope1,
                                       PlaneID const& pid2,
                                       double slope2) const
  {
    PlaneID target_plane = ThirdPlane(pid1, pid2);
    return ThirdPlaneSlope(pid1, slope1, pid2, slope2, target_plane);
  }

  //----------------------------------------------------------------------------
  double GeometryCore::ThirdPlane_dTdW(PlaneID const& pid1,
                                       double slope1,
                                       PlaneID const& pid2,
                                       double slope2,
                                       PlaneID const& output_plane) const
  {
    CheckIndependentPlanesOnSameTPC(pid1, pid2, "ThirdPlane_dTdW()");

    TPCGeo const& TPC = this->TPC(pid1);

    double angle[3], pitch[3];
    PlaneGeo const* const planes[3] = {
      &TPC.Plane(pid1), &TPC.Plane(pid2), &TPC.Plane(output_plane)};

    // We need wire pitch and "wire coordinate direction" for each plane.
    // The latter is perpendicular to the wire orientation.
    // PlaneGeo::PhiZ() defines the right orientation too.
    for (size_t i = 0; i < 3; ++i) {
      angle[i] = planes[i]->PhiZ();
      pitch[i] = planes[i]->WirePitch();
    }

    return ComputeThirdPlane_dTdW(
      angle[0], pitch[0], slope1, angle[1], pitch[1], slope2, angle[2], pitch[2]);
  }

  //----------------------------------------------------------------------------
  double GeometryCore::ThirdPlane_dTdW(PlaneID const& pid1,
                                       double slope1,
                                       PlaneID const& pid2,
                                       double slope2) const
  {
    PlaneID target_plane = ThirdPlane(pid1, pid2);
    return ThirdPlane_dTdW(pid1, slope1, pid2, slope2, target_plane);
  }

  //----------------------------------------------------------------------------
  // Given slopes dTime/dWire in two planes, return with the slope in the 3rd plane.
  // Requires slopes to be in the same metrics,
  // e.g. converted in a distances ratio.
  // B. Baller August 2014
  // Rewritten by T. Yang Apr 2015 using the equation in H. Greenlee's talk:
  // https://cdcvs.fnal.gov/redmine/attachments/download/1821/larsoft_apr20_2011.pdf
  // slide 2
  double GeometryCore::ComputeThirdPlaneSlope(double angle1,
                                              double slope1,
                                              double angle2,
                                              double slope2,
                                              double angle3)
  {
    // note that, if needed, the trigonometric functions can be pre-calculated.

    // Can't resolve very small slopes
    if ((std::abs(slope1) < 0.001) && (std::abs(slope2)) < 0.001) return 0.001;

    // We need the "wire coordinate direction" for each plane.
    // This is perpendicular to the wire orientation.
    double slope3 = 0.001;
    if (std::abs(slope1) > 0.001 && std::abs(slope2) > 0.001) {
      slope3 =
        (+(1. / slope1) * std::sin(angle3 - angle2) - (1. / slope2) * std::sin(angle3 - angle1)) /
        std::sin(angle1 - angle2);
    }
    if (slope3 != 0.)
      slope3 = 1. / slope3;
    else
      slope3 = 999.;

    return slope3;
  }

  //----------------------------------------------------------------------------
  double GeometryCore::ComputeThirdPlane_dTdW(double angle1,
                                              double pitch1,
                                              double dTdW1,
                                              double angle2,
                                              double pitch2,
                                              double dTdW2,
                                              double angle_target,
                                              double pitch_target)
  {
    // we need to convert dt/dw into homogeneous coordinates, and then back;
    // slope = [dT * (TDCperiod / driftVelocity)] / [dW * wirePitch]
    // The coefficient of dT is assumed to be the same for all the planes,
    // and it finally cancels out. Pitches cancel out only if they are all
    // the same.
    return pitch_target *
           ComputeThirdPlaneSlope(angle1, dTdW1 / pitch1, angle2, dTdW2 / pitch2, angle_target);
  }

  //......................................................................
  // This function is called if it is determined that two wires in a single TPC must overlap.
  // To determine the yz coordinate of the wire intersection, we need to know the
  // endpoints of both wires in xyz-space, and also their orientation (angle), and the
  // inner dimensions of the TPC frame.
  // Note: This calculation is entirely dependent  on an accurate GDML description of the TPC!
  // Mitch - Feb., 2011
  // Changed to use WireIDsIntersect(). It does not check whether the intersection is on both wires (the same as the old behavior). T. Yang - Apr, 2015
  //--------------------------------------------------------------------
  bool GeometryCore::IntersectionPoint(WireID const& wid1,
                                       WireID const& wid2,
                                       double& y,
                                       double& z) const
  {
    WireIDIntersection widIntersect;
    bool const found = WireIDsIntersect(wid1, wid2, widIntersect);
    y = widIntersect.y;
    z = widIntersect.z;
    return found;
  }

  //============================================================================
  //===  TPC set information
  //===
  //--------------------------------------------------------------------
  unsigned int GeometryCore::NTPCsets(readout::CryostatID const& cryoid) const
  {
    return fChannelMapAlg->NTPCsets(cryoid);
  }

  //--------------------------------------------------------------------
  unsigned int GeometryCore::MaxTPCsets() const
  {
    return fChannelMapAlg->MaxTPCsets();
  }

  //--------------------------------------------------------------------
  bool GeometryCore::HasTPCset(readout::TPCsetID const& tpcsetid) const
  {
    return fChannelMapAlg->HasTPCset(tpcsetid);
  }

  //--------------------------------------------------------------------
  readout::TPCsetID GeometryCore::FindTPCsetAtPosition(Point_t const& worldLoc) const
  {
    return TPCtoTPCset(FindTPCAtPosition(worldLoc));
  }

  //--------------------------------------------------------------------
  readout::TPCsetID GeometryCore::TPCtoTPCset(TPCID const& tpcid) const
  {
    return fChannelMapAlg->TPCtoTPCset(tpcid);
  }

  //--------------------------------------------------------------------
  std::vector<TPCID> GeometryCore::TPCsetToTPCs(readout::TPCsetID const& tpcsetid) const
  {
    return fChannelMapAlg->TPCsetToTPCs(tpcsetid);
  }

  //============================================================================
  //===  Readout plane information
  //===
  //--------------------------------------------------------------------
  unsigned int GeometryCore::NROPs(readout::TPCsetID const& tpcsetid) const
  {
    return fChannelMapAlg->NROPs(tpcsetid);
  }

  //--------------------------------------------------------------------
  unsigned int GeometryCore::MaxROPs() const
  {
    return fChannelMapAlg->MaxROPs();
  }

  //--------------------------------------------------------------------
  bool GeometryCore::HasROP(readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->HasROP(ropid);
  }

  //--------------------------------------------------------------------
  readout::ROPID GeometryCore::WirePlaneToROP(PlaneID const& planeid) const
  {
    return fChannelMapAlg->WirePlaneToROP(planeid);
  }

  //--------------------------------------------------------------------
  std::vector<PlaneID> GeometryCore::ROPtoWirePlanes(readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->ROPtoWirePlanes(ropid);
  }

  //--------------------------------------------------------------------
  std::vector<TPCID> GeometryCore::ROPtoTPCs(readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->ROPtoTPCs(ropid);
  }

  //--------------------------------------------------------------------
  raw::ChannelID_t GeometryCore::FirstChannelInROP(readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->FirstChannelInROP(ropid);
  }

  //--------------------------------------------------------------------
  View_t GeometryCore::View(readout::ROPID const& ropid) const
  {
    return View(fChannelMapAlg->FirstWirePlaneInROP(ropid));
  }

  //--------------------------------------------------------------------
  SigType_t GeometryCore::SignalType(readout::ROPID const& ropid) const
  {
    return fChannelMapAlg->SignalTypeForROPID(ropid);
  }

  //============================================================================
  //--------------------------------------------------------------------
  // Return gdml string which gives sensitive opdet name
  std::string GeometryCore::OpDetGeoName(CryostatID const& cid) const
  {
    return Cryostat(cid).OpDetGeoName();
  }

  //--------------------------------------------------------------------
  // Convert OpDet, Cryo into unique OpDet number
  unsigned int GeometryCore::OpDetFromCryo(unsigned int o, unsigned int c) const
  {
    static bool Loaded = false;
    static std::vector<unsigned int> LowestID;
    static unsigned int NCryo;

    CryostatID const cid{c};
    // If not yet loaded static parameters, do it
    if (Loaded == false) {

      Loaded = true;

      // Store the lowest ID for each cryostat
      NCryo = Ncryostats();
      LowestID.resize(NCryo + 1);
      LowestID.at(0) = 0;
      for (size_t cryo = 0; cryo != NCryo; ++cryo) {
        LowestID.at(cryo + 1) = LowestID.at(cryo) + Cryostat(cid).NOpDet();
      }
    }

    if ((c < NCryo) && (o < Cryostat(cid).NOpDet())) { return LowestID.at(c) + o; }

    throw cet::exception("OpDetCryoToOpID Error")
      << "Coordinates c=" << c << ", o=" << o << " out of range. Abort\n";
  }

  //--------------------------------------------------------------------
  const OpDetGeo& GeometryCore::OpDetGeoFromOpChannel(unsigned int OpChannel) const
  {
    return this->OpDetGeoFromOpDet(this->OpDetFromOpChannel(OpChannel));
  }

  //--------------------------------------------------------------------
  const OpDetGeo& GeometryCore::OpDetGeoFromOpDet(unsigned int OpDet) const
  {
    static bool Loaded = false;
    static std::vector<unsigned int> LowestID;
    static size_t NCryo;
    // If not yet loaded static parameters, do it
    if (Loaded == false) {

      Loaded = true;

      // Store the lowest ID for each cryostat
      NCryo = Ncryostats();
      LowestID.resize(NCryo + 1);
      LowestID[0] = 0;
      for (size_t cryo = 0; cryo != NCryo; ++cryo) {
        LowestID[cryo + 1] = LowestID[cryo] + Cryostat(CryostatID(cryo)).NOpDet();
      }
    }

    for (size_t i = 0; i != NCryo; ++i) {
      if ((OpDet >= LowestID[i]) && (OpDet < LowestID[i + 1])) {
        int c = i;
        int o = OpDet - LowestID[i];
        return this->Cryostat(CryostatID(c)).OpDet(o);
      }
    }
    // If we made it here, we didn't find the right combination. abort
    throw cet::exception("OpID To OpDetCryo error") << "OpID out of range, " << OpDet << "\n";
  }

  //--------------------------------------------------------------------
  // Find the closest OpChannel to this point, in the appropriate cryostat
  unsigned int GeometryCore::GetClosestOpDet(Point_t const& point) const
  {
    CryostatGeo const* cryo = PositionToCryostatPtr(point);
    if (!cryo) return std::numeric_limits<unsigned int>::max();
    int o = cryo->GetClosestOpDet(point);
    return OpDetFromCryo(o, cryo->ID().Cryostat);
  }

  //--------------------------------------------------------------------
  bool GeometryCore::WireIDIntersectionCheck(const WireID& wid1, const WireID& wid2) const
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
        << "1st wire " << wid1 << " does not exist (max wire number: " << Nwires(wid1.planeID())
        << ")";
      return false;
    }
    if (!HasWire(wid2)) {
      mf::LogError("WireIDIntersectionCheck")
        << "2nd wire " << wid2 << " does not exist (max wire number: " << Nwires(wid2.planeID())
        << ")";
      return false;
    }
    return true;
  }

  //--------------------------------------------------------------------
  //--- ROOTGeoNodeForwardIterator
  //---

  ROOTGeoNodeForwardIterator::ROOTGeoNodeForwardIterator(TGeoNode const* start_node)
  {
    if (start_node) {
      current_path.push_back({start_node, 0U});
      reach_deepest_descendant();
    }
  }

  ROOTGeoNodeForwardIterator& ROOTGeoNodeForwardIterator::operator++()
  {
    if (current_path.empty()) return *this;
    if (current_path.size() == 1) {
      current_path.pop_back();
      return *this;
    }

    // I am done; all my descendants were also done already;
    // first look at my younger siblings
    NodeInfo_t& current = current_path.back();
    NodeInfo_t const& parent = current_path[current_path.size() - 2];
    if (++(current.sibling) < parent.self->GetNdaughters()) {
      // my next sibling exists, let's parse his descendents
      current.self = parent.self->GetDaughter(current.sibling);
      reach_deepest_descendant();
    }
    else
      current_path.pop_back(); // no sibling, it's time for mum
    return *this;
  }

  //--------------------------------------------------------------------
  std::vector<TGeoNode const*> ROOTGeoNodeForwardIterator::get_path() const
  {
    std::vector<TGeoNode const*> node_path(current_path.size());
    std::transform(current_path.begin(),
                   current_path.end(),
                   node_path.begin(),
                   [](NodeInfo_t const& node_info) { return node_info.self; });
    return node_path;
  }

  //--------------------------------------------------------------------
  void ROOTGeoNodeForwardIterator::reach_deepest_descendant()
  {
    TGeoNode const* descendent = current_path.back().self;
    while (descendent->GetNdaughters() > 0) {
      descendent = descendent->GetDaughter(0);
      current_path.push_back({descendent, 0U});
    }
  }

} // namespace geo
