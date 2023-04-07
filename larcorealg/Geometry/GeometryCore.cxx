/**
 * @file   larcorealg/Geometry/GeometryCore.cxx
 * @brief  Access the description of detector geometry - implementation file
 * @author brebel@fnal.gov
 * @see    larcorealg/Geometry/GeometryCore.h
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
#include "larcorealg/Geometry/geo_vectors_utils.h"         // geo::vect
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h" // geo::vect

// Framework includes
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>

// C/C++ includes
#include <algorithm> // std::transform()
#include <cassert>
#include <cctype>   // ::tolower()
#include <cmath>    // std::abs() ...
#include <cstddef>  // size_t
#include <iterator> // std::back_inserter()
#include <limits>   // std::numeric_limits<>
#include <numeric>  // std::accumulate
#include <sstream>  // std::ostringstream
#include <utility>  // std::swap()
#include <vector>

namespace geo {

  template <typename T, typename Arg>
  auto bind(bool (T::*ft)(Arg const&, Arg const&) const, T* t)
  {
    return [t, ft](auto const& a, auto const& b) { return (t->*ft)(a, b); };
  }

  //......................................................................
  // Constructor.
  GeometryCore::GeometryCore(fhicl::ParameterSet const& pset,
                             std::unique_ptr<GeometryBuilder> builder,
                             std::unique_ptr<GeoObjectSorter> sorter)
    : Iteration{details::GeometryIterationPolicy{this}, details::ToGeometryElement{this}}
    , fSorter{std::move(sorter)}
    , fCompareAuxDets{bind(&GeoObjectSorter::compareAuxDets, fSorter.get())}
    , fCompareCryostats{bind(&GeoObjectSorter::compareCryostats, fSorter.get())}
    , fCompareTPCs{bind(&GeoObjectSorter::compareTPCs, fSorter.get())}
    , fCompareOpDets{bind(&GeoObjectSorter::compareOpDets, fSorter.get())}
    , fSurfaceY(pset.get<double>("SurfaceY"))
    , fDetectorName(pset.get<std::string>("Name"))
    , fPositionWiggle(pset.get<double>("PositionEpsilon", 1.e-4))
    , fBuilder{std::move(builder)}
  {
    std::transform(fDetectorName.begin(), fDetectorName.end(), fDetectorName.begin(), ::tolower);
  }

  //......................................................................
  void GeometryCore::LoadGeometryFile(std::string gdmlfile, std::string rootfile)
  {
    if (fManager) {
      throw cet::exception("GeometryCore") << "Reloading a geometry is not supported.\n";
    }

    if (gdmlfile.empty()) {
      throw cet::exception("GeometryCore") << "No GDML Geometry file specified!\n";
    }

    if (rootfile.empty()) {
      throw cet::exception("GeometryCore") << "No ROOT Geometry file specified!\n";
    }

    // [20210630, petrillo@slac.stanford.edu]
    // ROOT 6.22.08 allows us to choose the representation of lengths in the geometry
    // objects parsed from GDML.  In LArSoft we want them to be centimeters (ROOT
    // standard).  This was tracked as Redmine issue #25990, and I leave this mark
    // because I feel that we'll be back to it not too far in the future.  Despite the
    // documentation (ROOT 6.22/08), it seems the units are locked from the beginning,
    // so we unlock without prejudice.
    TGeoManager::LockDefaultUnits(false);
    TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);
    TGeoManager::LockDefaultUnits(true);

    TGeoManager::Import(rootfile.c_str());
    gGeoManager->LockGeometry();

    BuildGeometry();

    fManager = gGeoManager;
    fGDMLfile = std::move(gdmlfile);
    fROOTfile = std::move(rootfile);

    SortGeometry();

    mf::LogInfo("GeometryCore") << "New detector geometry loaded from "
                                << "\n\t" << fROOTfile << "\n\t" << fGDMLfile << "\n";
  }

  //......................................................................
  void GeometryCore::SortGeometry()
  {
    mf::LogInfo("GeometryCore") << "Sorting volumes...";

    std::sort(fAuxDets.begin(), fAuxDets.end(), fCompareAuxDets);
    std::sort(fCryostats.begin(), fCryostats.end(), fCompareCryostats);

    // Renumber cryostats according to sorted order
    CryostatID::CryostatID_t c = 0;
    for (CryostatGeo& cryo : fCryostats) {
      cryo.SortSubVolumes(fCompareTPCs, fCompareOpDets);
      cryo.UpdateAfterSorting(CryostatID{c});
      ++c;
    }
  }

  //......................................................................
  TGeoManager* GeometryCore::ROOTGeoManager() const
  {
    assert(fManager == gGeoManager);
    return fManager;
  }

  //......................................................................
  unsigned int GeometryCore::NOpDets() const
  {
    int N = 0;
    for (auto const& cryo : Iterate<CryostatGeo>())
      N += cryo.NOpDet();
    return N;
  }

  //......................................................................
  unsigned int GeometryCore::NAuxDetSensitive(size_t const& aid) const
  {
    if (aid < NAuxDets()) { return AuxDets()[aid].NSensitiveVolume(); }
    throw cet::exception("Geometry")
      << "Requested AuxDet index " << aid << " is out of range: " << NAuxDets();
  }

  //......................................................................
  //
  // Return the geometry description of the ith cryostat in the detector.
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
  AuxDetGeo const& GeometryCore::AuxDet(unsigned int const ad) const
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
    return cryo->PositionToTPCID(point, 1. + fPositionWiggle);
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
  CryostatGeo const& GeometryCore::PositionToCryostat(Point_t const& point) const
  {
    if (auto cstat = PositionToCryostatPtr(point)) { return *cstat; }
    throw cet::exception("GeometryCore") << "Can't find any cryostat at position " << point << "\n";
  }

  //......................................................................
  std::string const& GeometryCore::GetWorldVolumeName() const
  {
    // For now, and possibly forever, this is a constant (given the definition of
    // "nodeNames" above).
    static std::string const worldVolumeName{"volWorld"};
    return worldVolumeName;
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
    double const halfwidth = pBox->GetDX();
    double const halfheight = pBox->GetDY();
    double const halflength = pBox->GetDZ();

    return {trans.LocalToWorld(Point_t{-halfwidth, -halfheight, -halflength}),
            trans.LocalToWorld(Point_t{+halfwidth, +halfheight, +halflength})};
  }

  //......................................................................
  /** **************************************************************************
   * @brief Iterator to navigate through all the nodes
   *
   * Note that this is not a fully standard forward iterator in that it lacks of the
   * postfix operator. The reason is that it's too expensive and it should be avoided.
   *
   * Also I did not bother declaring the standard type definitions (that's just laziness).
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
   * These iterators are one use only, and they can't be reset after a loop is completed.
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
  };

  struct CollectNodesByName {
    std::vector<TGeoNode const*> nodes;

    CollectNodesByName(std::set<std::string> const& names) : matcher(names) {}

    /// If the name of the node matches, records the end node
    void operator()(TGeoNode const& node)
    {
      if (matcher(node)) nodes.push_back(&node);
    }

    void operator()(ROOTGeoNodeForwardIterator const& iter) { operator()(**iter); }

  private:
    NodeNameMatcherClass matcher;
  };

  struct CollectPathsByName {
    std::vector<std::vector<TGeoNode const*>> paths;

    CollectPathsByName(std::set<std::string> const& names) : matcher(names) {}

    /// If the name of the node matches, records the node full path
    void operator()(ROOTGeoNodeForwardIterator const& iter)
    {
      if (matcher(**iter)) paths.push_back(iter.get_path());
    }

  private:
    NodeNameMatcherClass matcher;
  };

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
  unsigned int GeometryCore::MaxTPCs() const
  {
    unsigned int maxTPCs = 0;
    for (CryostatGeo const& cryo : fCryostats) {
      unsigned int maxTPCsInCryo = cryo.NTPC();
      if (maxTPCsInCryo > maxTPCs) maxTPCs = maxTPCsInCryo;
    }
    return maxTPCs;
  }

  //......................................................................
  unsigned int GeometryCore::TotalNTPC() const
  {
    // it looks like C++11 lambdas have made STL algorithms easier to use, but only so
    // much:
    return std::accumulate(
      fCryostats.begin(), fCryostats.end(), 0U, [](unsigned int sum, CryostatGeo const& cryo) {
        return sum + cryo.NTPC();
      });
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
  void GeometryCore::BuildGeometry()
  {
    std::cerr << "GeometryCore top node:" << gGeoManager->GetTopNode() << '\n';
    GeoNodePath path{gGeoManager->GetTopNode()};
    fCryostats = fBuilder->extractCryostats(path);
    fAuxDets = fBuilder->extractAuxiliaryDetectors(path);
  }

  //......................................................................
  //
  // Return the total mass of the detector
  //
  //
  double GeometryCore::TotalMass(std::string vol) const
  {
    //the TGeoNode::GetVolume() returns the TGeoVolume of the detector outline and ROOT
    //calculates the mass in kg for you
    TGeoVolume* gvol = gGeoManager->FindVolumeFast(vol.c_str());
    if (gvol) return gvol->Weight();

    throw cet::exception("GeometryCore")
      << "could not find specified volume '" << vol << " 'to determine total mass\n";
  }

  //......................................................................
  double GeometryCore::MassBetweenPoints(Point_t const& p1, Point_t const& p2) const
  {
    //The purpose of this method is to determine the column density between the two points
    //given.  Do that by starting at p1 and stepping until you get to the node of p2.
    //calculate the distance between the point just inside that node and p2 to get the
    //last bit of column density
    double columnD = 0.;

    //first initialize a track - get the direction cosines
    Vector_t const dir = (p2 - p1).Unit();

    double const dxyz[3] = {dir.X(), dir.Y(), dir.Z()};
    double const cp1[3] = {p1.X(), p1.Y(), p1.Z()};
    gGeoManager->InitTrack(cp1, dxyz);

    //might be helpful to have a point to a TGeoNode
    TGeoNode* node = gGeoManager->GetCurrentNode();

    //check that the points are not in the same volume already.  if they are in different
    //volumes, keep stepping until you are in the same volume as the second point
    while (!gGeoManager->IsSameLocation(p2.X(), p2.Y(), p2.Z())) {
      gGeoManager->FindNextBoundary();
      columnD += gGeoManager->GetStep() * node->GetMedium()->GetMaterial()->GetDensity();

      //the act of stepping puts you in the next node and returns that node
      node = gGeoManager->Step();
    } //end loop to get to volume of second point

    //now you are in the same volume as the last point, but not at that point.  get the
    //distance between the current point and the last one
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

  //============================================================================
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
  OpDetGeo const& GeometryCore::OpDetGeoFromOpDet(unsigned int OpDet) const
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
        return Cryostat(CryostatID(c)).OpDet(o);
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

    // I am done; all my descendants were also done already; first look at my younger
    // siblings
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
