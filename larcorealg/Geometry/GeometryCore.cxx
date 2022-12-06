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
#include <cctype>    // ::tolower()
#include <cmath>     // std::abs() ...
#include <cstddef>   // size_t
#include <iterator>  // std::back_inserter()
#include <limits>    // std::numeric_limits<>
#include <numeric>   // std::accumulate
#include <sstream>   // std::ostringstream
#include <utility>   // std::swap()
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
  GeometryCore::GeometryCore(fhicl::ParameterSet const& pset,
                             std::unique_ptr<GeoObjectSorter const> sorter)
    : Iteration{details::GeometryIterationPolicy{this}, ToGeometryElement{this}}
    , fSorter{std::move(sorter)}
    , fSurfaceY(pset.get<double>("SurfaceY"))
    , fDetectorName(pset.get<std::string>("Name"))
    , fMinWireZDist(pset.get<double>("MinWireZDist", 3.0))
    , fPositionWiggle(pset.get<double>("PositionEpsilon", 1.e-4))
    , fBuilderParameters(pset.get<fhicl::ParameterSet>("Builder", {}))
  {
    std::transform(fDetectorName.begin(), fDetectorName.end(), fDetectorName.begin(), ::tolower);
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

    fGDMLfile = std::move(gdmlfile);
    fROOTfile = std::move(rootfile);

    SortGeometry();

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
    fCryostats = {};
    fAuxDets = {};
  }

  //......................................................................
  void GeometryCore::SortGeometry()
  {
    mf::LogInfo("GeometryCore") << "Sorting volumes...";

    fSorter->SortAuxDets(fAuxDets);
    fSorter->SortCryostats(fCryostats);

    // Renumber cryostats according to sorted order
    CryostatID::CryostatID_t c = 0;
    for (CryostatGeo& cryo : fCryostats) {
      cryo.SortSubVolumes(*fSorter);
      cryo.UpdateAfterSorting(CryostatID{c});
      ++c;
    }

    // Update views
    std::set<View_t> updatedViews;
    for (auto const& tpc : Iterate<TPCGeo>()) {
      auto const& TPCviews = tpc.Views();
      updatedViews.insert(TPCviews.cbegin(), TPCviews.cend());
    }
    allViews = move(updatedViews);
  }

  //......................................................................
  TGeoManager* GeometryCore::ROOTGeoManager() const { return gGeoManager; }

  //......................................................................
  unsigned int GeometryCore::NOpDets() const
  {
    int N = 0;
    for (auto const& cryo : Iterate<geo::CryostatGeo>())
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
  // Number of different views, or wire orientations
  unsigned int GeometryCore::Nviews() const { return MaxPlanes(); }

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
    // For now, and possibly forever, this is a constant (given the
    // definition of "nodeNames" above).
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

  private:
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

  private:
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
  // This method returns the distance between wires in the specified view
  // it assumes all planes of a given view have the same pitch
  double GeometryCore::WireAngleToVertical(View_t view, TPCID const& tpcid) const
  {
    // loop over the planes in cryostat 0, tpc 0 to find the plane with the
    // specified view
    TPCGeo const& tpc = TPC(tpcid);
    for (unsigned int p = 0; p < tpc.Nplanes(); ++p) {
      PlaneGeo const& plane = tpc.Plane(p);
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
    fCryostats = builder.extractCryostats(path);
    fAuxDets = builder.extractAuxiliaryDetectors(path);
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

  //......................................................................
  std::optional<WireIDIntersection> GeometryCore::WireIDsIntersect(const WireID& wid1,
                                                                   const WireID& wid2) const
  {
    if (!WireIDIntersectionCheck(wid1, wid2)) { return std::nullopt; }

    // get the endpoints to see if wires intersect
    Segment_t const w1 = WireEndPoints(wid1);
    Segment_t const w2 = WireEndPoints(wid2);

    // TODO extract the coordinates in the right way;
    // is it any worth, since then the result is in (y, z), whatever it means?
    WireIDIntersection result;
    bool const cross = IntersectLines(w1.start().Y(),
                                      w1.start().Z(),
                                      w1.end().Y(),
                                      w1.end().Z(),
                                      w2.start().Y(),
                                      w2.start().Z(),
                                      w2.end().Y(),
                                      w2.end().Z(),
                                      result.y,
                                      result.z);
    if (!cross) { return std::nullopt; }
    bool const within = lar::util::PointWithinSegments(w1.start().Y(),
                                                       w1.start().Z(),
                                                       w1.end().Y(),
                                                       w1.end().Z(),
                                                       w2.start().Y(),
                                                       w2.start().Z(),
                                                       w2.end().Y(),
                                                       w2.end().Z(),
                                                       result.y,
                                                       result.z);

    result.TPC = (within ? wid1.TPC : TPCID::InvalidID);
    return result;
  }

  //......................................................................
  bool GeometryCore::WireIDsIntersect(const WireID& wid1,
                                      const WireID& wid2,
                                      Point_t& intersection) const
  {
    // This is not a real 3D intersection: the wires do not cross, since they
    // are required to belong to two different planes.
    //
    // We take the point on the first wire which is closest to the
    // other one.
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

    return std::abs(intersectionAndOffset.offset1) <= wire1.HalfL() &&
           std::abs(intersectionAndOffset.offset2) <= wire2.HalfL();
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

    TPCGeo const& tpc = TPC(pid1);

    // We need the "wire coordinate direction" for each plane.
    // This is perpendicular to the wire orientation.
    // PlaneGeo::PhiZ() defines the right orientation too.
    return ComputeThirdPlaneSlope(tpc.Plane(pid1).PhiZ(),
                                  slope1,
                                  tpc.Plane(pid2).PhiZ(),
                                  slope2,
                                  tpc.Plane(output_plane).PhiZ());
  }

  //----------------------------------------------------------------------------
  double GeometryCore::ThirdPlaneSlope(PlaneID const& pid1,
                                       double slope1,
                                       PlaneID const& pid2,
                                       double slope2) const
  {
    return ThirdPlaneSlope(pid1, slope1, pid2, slope2, ThirdPlane(pid1, pid2));
  }

  //----------------------------------------------------------------------------
  double GeometryCore::ThirdPlane_dTdW(PlaneID const& pid1,
                                       double slope1,
                                       PlaneID const& pid2,
                                       double slope2,
                                       PlaneID const& output_plane) const
  {
    CheckIndependentPlanesOnSameTPC(pid1, pid2, "ThirdPlane_dTdW()");

    TPCGeo const& tpc = TPC(pid1);

    double angle[3], pitch[3];
    PlaneGeo const* const planes[3] = {
      &tpc.Plane(pid1), &tpc.Plane(pid2), &tpc.Plane(output_plane)};

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
    return ThirdPlane_dTdW(pid1, slope1, pid2, slope2, ThirdPlane(pid1, pid2));
  }

  //----------------------------------------------------------------------------
  // Given slopes dTime/dWire in two planes, return with the slope in the 3rd plane.
  // Requires slopes to be in the same metrics,
  // e.g. converted in a distances ratio.
  // Note: Uses equation in H. Greenlee's talk:
  //       https://cdcvs.fnal.gov/redmine/attachments/download/1821/larsoft_apr20_2011.pdf
  //       slide 2
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
