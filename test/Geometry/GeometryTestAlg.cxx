/**
 * @file   GeometryTestAlg.cxx
 * @brief  Unit test for geometry functionalities: implementation file
 * @date   2011/02/17
 * @author brebel@fnal.gov
 * @see    GeometryTestAlg.h
 */

// our header
#include "GeometryTestAlg.h"

// LArSoft includes
#include "larcorealg/CoreUtils/DumpUtils.h" // lar::dump::vector3D(), ...
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/Decomposer.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/SimpleGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/geo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi<>
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"          // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h" // readout::ROPID

// Framework includes
#include "cetlib/pow.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "RtypesCore.h" // for kTRUE
#include "TClass.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"
#include "TObjArray.h"
#include "TStopwatch.h"

// C/C++ standard libraries
#include <algorithm> // std::copy()
#include <array>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <iterator> // std::inserter()
#include <limits>   // std::numeric_limits<>
#include <memory>
#include <set>
#include <sstream>
#include <stdint.h>
#include <stdlib.h> // for abort
#include <string>
#include <vector>

//------------------------------------------------------------------------------
namespace {

  /// Returns whether the CET exception e contains the specified category cat
  bool hasCategory(cet::exception const& e, std::string const& cat)
  {
    for (auto const& e_category : e.history())
      if (e_category == cat) return true;
    return false;
  } // hasCategory()

  /// Returns a convenience string for the specified direction
  template <typename Vector>
  std::string directionName(Vector const& v)
  {
    if (v == geo::Xaxis<Vector>()) return "x";
    if (v == geo::Yaxis<Vector>()) return "y";
    if (v == geo::Zaxis<Vector>()) return "z";
    if (v == -geo::Xaxis<Vector>()) return "-x";
    if (v == -geo::Yaxis<Vector>()) return "-y";
    if (v == -geo::Zaxis<Vector>()) return "-z";
    std::ostringstream sstr;
    sstr << v;
    return sstr.str();
  } // directionName()

} // local namespace

namespace geo {

  //......................................................................
  GeometryTestAlg::GeometryTestAlg(fhicl::ParameterSet const& pset)
    : geom(nullptr)
    , fDisableValidWireIDcheck(pset.get<bool>("DisableWireBoundaryCheck", false))
    , fExpectedWirePitches(pset.get<std::vector<double>>("ExpectedWirePitches", {}))
    , fExpectedPlanePitches(pset.get<std::vector<double>>("ExpectedPlanePitches", {}))
    , fComputeMass(pset.get("ComputeMass", true))
  {
    // initialize the list of non-fatal exceptions
    std::vector<std::string> NonFatalErrors(
      pset.get<std::vector<std::string>>("ForgiveExceptions", std::vector<std::string>()));
    std::copy(NonFatalErrors.begin(),
              NonFatalErrors.end(),
              std::inserter(fNonFatalExceptions, fNonFatalExceptions.end()));

    // initialize the list of tests to be run
    //
    // our name selector accepts everything by default;
    // the default set skips the following:
    fRunTests.AddToDefinition("default", "-CheckOverlaps", "-ThoroughCheck", "-PrintWires");
    fRunTests.ParseNames("@default"); // let's start from default

    // read and parse the test list from the configuration parameters (if any)
    fRunTests.Parse(pset.get<std::vector<std::string>>("RunTests", {}));

    std::ostringstream sstr;
    fRunTests.PrintConfiguration(sstr);

    mf::LogInfo("GeometryTestAlg") << "Test selection:" << sstr.str();

  } // GeometryTestAlg::GeometryTestAlg()

  //......................................................................
  unsigned int GeometryTestAlg::Run()
  {

    if (!geom) {
      throw cet::exception("GeometryTestAlg")
        << "GeometryTestAlg not configured: no valid geometry provided.\n";
    }

    unsigned int nErrors = 0; // currently unused

    // change the printed version number when changing the "GeometryTest" output
    //
    // Version 1.1:
    //   more TPC information when printing all geometry
    //
    mf::LogVerbatim("GeometryTest") << "GeometryTest version 1.1";

    mf::LogInfo("GeometryTestInfo") << "Running on detector: '" << geom->DetectorName() << "'";

    mf::LogVerbatim("GeometryTest") << "  Running on detector: '" << geom->DetectorName() << "'"
                                    << "\nGeometry file: " << geom->ROOTFile();

    try {
      if (shouldRunTests("DetectorIntro")) {
        MF_LOG_INFO("GeometryTest") << "detector introduction:";
        printDetectorIntro();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("CheckOverlaps")) {
        MF_LOG_INFO("GeometryTest") << "test for overlaps ...";
        gGeoManager->CheckOverlaps(1e-5);
        gGeoManager->PrintOverlaps();
        if (!gGeoManager->GetListOfOverlaps()->IsEmpty()) {
          mf::LogError("GeometryTest") << gGeoManager->GetListOfOverlaps()->GetSize()
                                       << " overlaps found in geometry during overlap test!";
          ++nErrors;
        }
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("ThoroughCheck")) {
        MF_LOG_INFO("GeometryTest") << "thorough geometry test ...";
        gGeoManager->CheckGeometryFull();
        if (!gGeoManager->GetListOfOverlaps()->IsEmpty()) {
          mf::LogError("GeometryTest") << gGeoManager->GetListOfOverlaps()->GetSize()
                                       << " overlaps found in geometry during thorough test!";
          ++nErrors;
        }
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("Cryostat")) {
        MF_LOG_INFO("GeometryTest") << "test Cryostat methods ...";
        testCryostat();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("FindVolumes")) {
        MF_LOG_INFO("GeometryTest") << "test FindAllVolumes method ...";
        testFindVolumes();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("PlaneDirections")) {
        MF_LOG_INFO("GeometryTest") << "test plane directions...";
        testPlaneDirections();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("WireOrientations")) {
        MF_LOG_INFO("GeometryTest") << "test wire orientations...";
        testWireOrientations();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("ChannelToROP")) {
        MF_LOG_INFO("GeometryTest") << "test channel to ROP and back ...";
        testChannelToROP();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("ChannelToWire")) {
        MF_LOG_INFO("GeometryTest") << "test channel to plane wire and back ...";
        testChannelToWire();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("FindPlaneCenters")) {
        MF_LOG_INFO("GeometryTest") << "test find plane centers...";
        testFindPlaneCenters();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("WireCoordFromPlane")) {
        MF_LOG_INFO("GeometryTest") << "test PlaneGeo::WireCoordinate...";
        testWireCoordFromPlane();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("ParallelWires")) {
        MF_LOG_INFO("GeometryTest") << "test wire parallelism...";
        testParallelWires();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("PlanePointDecomposition")) {
        MF_LOG_INFO("GeometryTest") << "test plane point decomposition...";
        testPlanePointDecomposition();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("PlaneProjections")) {
        MF_LOG_INFO("GeometryTest") << "test PlaneGeo::PointProjection...";
        testPlaneProjection();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("WireCoordAngle")) {
        MF_LOG_INFO("GeometryTest") << "testWireCoordAngle...";
        testWireCoordAngle();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("Projection")) {
        MF_LOG_INFO("GeometryTest") << "testProject...";
        testProject();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("WirePos")) {
        MF_LOG_INFO("GeometryTest") << "testWirePos...";
        // There is a contradiction here, and these must be tested differently
        // Testing based on detector ID should NOT become common practice
        MF_LOG_INFO("GeometryTest") << "disabled.";
      }

      if (shouldRunTests("NearestWire")) {
        MF_LOG_INFO("GeometryTest") << "testNearestWire...";
        testNearestWire();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("WireIntersection")) {
        MF_LOG_INFO("GeometryTest") << "testWireIntersection...";
        testWireIntersection();
        MF_LOG_INFO("GeometryTest") << "testWireIntersection complete";
      }

      if (shouldRunTests("ThirdPlane")) {
        MF_LOG_INFO("GeometryTest") << "testThirdPlane...";
        testThirdPlane();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("ThirdPlaneSlope")) {
        MF_LOG_INFO("GeometryTest") << "testThirdPlaneSlope...";
        testThirdPlane_dTdW();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("WirePitch")) {
        MF_LOG_INFO("GeometryTest") << "testWirePitch...";
        testWirePitch();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("InterWireProjectedDistance")) {
        MF_LOG_INFO("GeometryTest") << "testInterWireProjectedDistance...";
        testInterWireProjectedDistance();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("PlanePitch")) {
        MF_LOG_INFO("GeometryTest") << "testPlanePitch...";
        testPlanePitch();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("Stepping")) {
        MF_LOG_INFO("GeometryTest") << "testStepping...";
        testStepping();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("FindAuxDet")) {
        MF_LOG_INFO("GeometryTest") << "testFindAuxDet...";
        testFindAuxDet();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }

      if (shouldRunTests("PrintWires")) {
        MF_LOG_INFO("GeometryTest") << "printAllGeometry...";
        printAllGeometry();
        MF_LOG_INFO("GeometryTest") << "complete.";
      }
    }
    catch (cet::exception& e) {
      mf::LogWarning("GeometryTest") << "exception caught: \n" << e;
      if (fNonFatalExceptions.count(e.category()) == 0) throw;
    }

    std::ostringstream out;
    if (!fRunTests.CheckQueryRegistry(out)) {
      throw cet::exception("GeometryTest") << "(postumous) configuration error detected!\n"
                                           << out.str();
    }

    mf::LogInfo log("GeometryTest");
    log << "Tests completed:";
    auto tests_run = fRunTests.AcceptedNames();
    if (tests_run.empty()) { log << "\n  no test run"; }
    else {
      log << "\n  " << tests_run.size() << " tests run:\t ";
      for (std::string const& test_name : tests_run)
        log << " " << test_name;
    }
    auto tests_skipped = fRunTests.RejectedNames();
    if (!tests_skipped.empty()) {
      log << "\n  " << tests_skipped.size() << " tests skipped:\t ";
      for (std::string const& test_name : tests_skipped)
        log << " " << test_name;
    }

    return nErrors;
  } // GeometryTestAlg::Run()

  //......................................................................
  void GeometryTestAlg::printDetectorIntro() const
  {

    geo::WireGeo const& testWire = geom->Wire(geo::WireID(0, 0, 1, 10));
    mf::LogVerbatim log("GeometryTest");
    log << "Wire Rmax  " << testWire.RMax() << "\nWire length " << 2. * testWire.HalfL()
        << "\nWire Rmin  " << testWire.RMin();

    if (fComputeMass) { log << "\nTotal mass " << geom->TotalMass(); }

    log << "\nNumber of views " << geom->Nviews() << "\nNumber of channels " << geom->Nchannels()
        << "\nMaximum number of:"
        << "\n  TPC in a cryostat: " << geom->MaxTPCs()
        << "\n  planes in a TPC:   " << geom->MaxPlanes()
        << "\n  wires in a plane:  " << geom->MaxWires() << "\nTotal number of TPCs "
        << geom->TotalNTPC() << "\nAuxiliary detectors  " << geom->NAuxDets();

  } // GeometryTestAlg::printDetectorIntro()

  //......................................................................
  void GeometryTestAlg::printChannelSummary()
  {
    static unsigned int OneSeg = 0;
    static unsigned int TwoSegs = 0;
    static unsigned int ThreeSegs = 0;
    static unsigned int FourSegs = 0;
    uint32_t channels = geom->Nchannels();
    auto const ntpc = geom->NTPC(CryostatID{0});
    if (ntpc > 1) channels /= ntpc / 2;

    for (uint32_t c = 0; c < channels; c++) {

      unsigned int ChanSize = geom->ChannelToWire(c).size();

      if (ChanSize == 1)
        ++OneSeg;
      else if (ChanSize == 2)
        ++TwoSegs;
      else if (ChanSize == 3)
        ++ThreeSegs;
      else if (ChanSize == 4)
        ++FourSegs;
    }

    mf::LogVerbatim("GeometryTest") << "OneSeg: " << OneSeg << ",  TwoSegs: " << TwoSegs
                                    << ",  ThreeSegs: " << ThreeSegs << ",  FourSegs: " << FourSegs;
  }

  //......................................................................
  void GeometryTestAlg::printVolBounds()
  {
    geo::CryostatGeo::LocalPoint_t const cryo_origin{};
    for (auto const& cryostat : geom->Iterate<geo::CryostatGeo>()) {
      auto const world = cryostat.toWorldCoords(cryo_origin);

      mf::LogVerbatim("GeometryTest") << "Cryo " << cryostat.ID();
      mf::LogVerbatim("GeometryTest") << "    -x: " << world.X() - cryostat.HalfWidth();
      mf::LogVerbatim("GeometryTest") << "    +x: " << world.X() + cryostat.HalfWidth();
      mf::LogVerbatim("GeometryTest") << "    -y: " << world.Y() - cryostat.HalfHeight();
      mf::LogVerbatim("GeometryTest") << "    +y: " << world.Y() + cryostat.HalfHeight();
      mf::LogVerbatim("GeometryTest") << "    -z: " << world.Z() - cryostat.Length() / 2;
      mf::LogVerbatim("GeometryTest") << "    +z: " << world.Z() + cryostat.Length() / 2;

      geo::TPCGeo::LocalPoint_t const tpc_origin{};
      for (unsigned int t = 0; t < geom->NTPC(cryostat.ID()); ++t) {
        auto const world = cryostat.TPC(t).toWorldCoords(tpc_origin);

        mf::LogVerbatim("GeometryTest") << "  TPC " << t;
        mf::LogVerbatim("GeometryTest") << "    -x: " << world.X() - cryostat.TPC(t).HalfWidth();
        mf::LogVerbatim("GeometryTest") << "    +x: " << world.X() + cryostat.TPC(t).HalfWidth();
        mf::LogVerbatim("GeometryTest") << "    -y: " << world.Y() - cryostat.TPC(t).HalfHeight();
        mf::LogVerbatim("GeometryTest") << "    +y: " << world.Y() + cryostat.TPC(t).HalfHeight();
        mf::LogVerbatim("GeometryTest") << "    -z: " << world.Z() - cryostat.TPC(t).Length() / 2;
        mf::LogVerbatim("GeometryTest") << "    +z: " << world.Z() + cryostat.TPC(t).Length() / 2;
      }
    }
  }

  //......................................................................
  // great sanity check for geometry, only call in analyze when debugging
  void GeometryTestAlg::printDetDim()
  {
    for (auto const& cid : geom->Iterate<CryostatID>()) {
      mf::LogVerbatim("GeometryTest") << "Cryo " << cid;
      mf::LogVerbatim("GeometryTest") << "    width: " << geom->CryostatHalfWidth(cid);
      mf::LogVerbatim("GeometryTest") << "    height: " << geom->CryostatHalfHeight(cid);
      mf::LogVerbatim("GeometryTest") << "    length: " << geom->CryostatLength(cid);

      TPCID const tpcid{cid, 0};
      mf::LogVerbatim("GeometryTest") << "  TPC 0";
      mf::LogVerbatim("GeometryTest") << "    width: " << geom->DetHalfWidth(tpcid);
      mf::LogVerbatim("GeometryTest") << "    height: " << geom->DetHalfHeight(tpcid);
      mf::LogVerbatim("GeometryTest") << "    length: " << geom->DetLength(tpcid);
    }
  }

  //......................................................................
  // great sanity check for volume sorting, only call in analyze when debugging
  void GeometryTestAlg::printWirePos()
  {
    CryostatID const cs{0};
    for (unsigned int t = 0; t < std::floor(geom->NTPC(cs) / 12) + 1; ++t) {
      for (unsigned int p = 0; p < 3; ++p) {
        auto const& plane = geom->Cryostat(cs).TPC(t).Plane(p);
        for (unsigned int w = 0; w < plane.Nwires(); w++) {
          auto const xyz = plane.Wire(w).GetCenter();
          std::cout << "WireID (" << cs << ", " << t << ", " << p << ", " << w
                    << "):  x = " << xyz.X() << ", y = " << xyz.Y() << ", z = " << xyz.Z() << '\n';
        }
      }
    }
  }

  //......................................................................
  // great insanity: print all wires in a TPC
  void GeometryTestAlg::printWiresInTPC(const geo::TPCGeo& tpc, std::string indent /* = "" */) const
  {
    const unsigned int nPlanes = tpc.Nplanes();

    tpc.PrintTPCInfo(mf::LogVerbatim("GeometryTest") << indent, indent, geo::TPCGeo::MaxVerbosity);

    for (unsigned int p = 0; p < nPlanes; ++p) {
      const geo::PlaneGeo& plane = tpc.Plane(p);
      const unsigned int nWires = plane.Nwires();

      plane.PrintPlaneInfo(mf::LogVerbatim("GeometryTest") << indent << "  ",
                           indent + "      ",
                           8 /* large verbosity */
      );

      for (unsigned int w = 0; w < nWires; ++w) {
        const geo::WireGeo& wire = plane.Wire(w);
        // this additional check is preserved to test alternative transformation
        // code paths; center is expected to match wire.GetCenter()
        geo::WireGeo::LocalPoint_t const local{};
        auto const center = wire.toWorldCoords(local);

        // the wire should be aligned on z axis, half on each side of 0,
        // in its local frame
        mf::LogVerbatim("GeometryTest")
          << indent << "    wire #" << w << " at " << lar::dump::array<3>(center) << "\n"
          << indent << "       start at " << lar::dump::vector3D(wire.GetStart()) << "\n"
          << indent << "      middle at " << lar::dump::vector3D(wire.GetCenter()) << "\n"
          << indent << "         end at " << lar::dump::vector3D(wire.GetEnd());
      } // for wire
    }   // for plane
  }

  void GeometryTestAlg::printAllGeometry() const
  {
    const unsigned int nCryostats = geom->Ncryostats();
    mf::LogVerbatim("GeometryTest")
      << "Detector " << geom->DetectorName() << " has " << nCryostats << " cryostats:";
    for (auto const& cryostat : geom->Iterate<geo::CryostatGeo>()) {
      const unsigned int nTPCs = cryostat.NTPC();
      mf::LogVerbatim("GeometryTest")
        << "  cryostat #" << cryostat.ID() << " at " << lar::dump::vector3D(cryostat.GetCenter())
        << " cm has " << nTPCs << " TPC(s):";
      for (unsigned int t = 0; t < nTPCs; ++t) {
        const geo::TPCGeo& tpc = cryostat.TPC(t);
        printWiresInTPC(tpc, "    ");
      } // for TPC
    }   // for cryostat
    printAuxiliaryDetectors();
    mf::LogVerbatim("GeometryTest") << "End of detector " << geom->DetectorName() << " geometry.";
  }

  //......................................................................
  void GeometryTestAlg::printAuxiliaryDetectors() const
  {
    mf::LogVerbatim log("GeometryTest");

    unsigned int const nAuxDets = geom->NAuxDets();
    log << "There are " << nAuxDets << " auxiliary detectors:";
    for (unsigned int iDet = 0; iDet < nAuxDets; ++iDet) {
      log << "\n[#" << iDet << "] ";
      printAuxDetGeo(log, geom->AuxDet(iDet), "  ", "");
    } // for
  }

  //......................................................................
  template <typename Stream>
  void GeometryTestAlg::printAuxDetGeo(Stream&& out,
                                       geo::AuxDetGeo const& auxDet,
                                       std::string indent,
                                       std::string firstIndent) const
  {

    lar::util::RealComparisons<double> coordIs(1e-4);

    auto const center = auxDet.GetCenter();
    auto const normal = auxDet.GetNormalVector();

    out << firstIndent << "\"" << auxDet.Name() << "\" centered at " << lar::dump::array<3U>(center)
        << " cm, size ( " << (2.0 * auxDet.HalfWidth1());
    if (coordIs.nonEqual(auxDet.HalfWidth1(), auxDet.HalfWidth2()))
      out << "/" << (2.0 * auxDet.HalfWidth2());
    out << " x " << (2.0 * auxDet.HalfHeight()) << " x " << auxDet.Length() << " ) cm"
        << ", normal facing " << lar::dump::array<3U>(normal);
    unsigned int nSensitive = auxDet.NSensitiveVolume();
    switch (nSensitive) {
    case 0: break;
    case 1:
      out << "\n" << indent << "with sensitive volume ";
      printAuxDetSensitiveGeo(
        std::forward<Stream>(out), auxDet.SensitiveVolume(0U), indent + "  ", "");
      break;
    default:
      out << "\n" << indent << "with " << auxDet.NSensitiveVolume() << " sensitive detectors:";
      for (unsigned int iSens = 0; iSens < nSensitive; ++iSens) {
        out << "\n" << indent << "  [#" << iSens << "] ";
        printAuxDetSensitiveGeo(
          std::forward<Stream>(out), auxDet.SensitiveVolume(iSens), indent + "  ", "");
      } // for
      break;
    } // if sensitive detectors

  } // GeometryTestAlg::printAuxDetGeo()

  //......................................................................
  template <typename Stream>
  void GeometryTestAlg::printAuxDetSensitiveGeo(Stream&& out,
                                                geo::AuxDetSensitiveGeo const& auxDetSens,
                                                std::string /* indent */,
                                                std::string firstIndent) const
  {

    lar::util::RealComparisons<double> coordIs(1e-4);

    auto const center = auxDetSens.GetCenter();
    auto const normal = auxDetSens.GetNormalVector();

    out << firstIndent << "centered at " << lar::dump::array<3U>(center) << " cm, size ( "
        << (2.0 * auxDetSens.HalfWidth1());
    if (coordIs.nonEqual(auxDetSens.HalfWidth1(), auxDetSens.HalfWidth2()))
      out << "/" << (2.0 * auxDetSens.HalfWidth2());
    out << " x " << (2.0 * auxDetSens.HalfHeight()) << " x " << auxDetSens.Length() << " ) cm"
        << ", normal facing " << lar::dump::array<3U>(normal);

  } // GeometryTestAlg::printAuxDetSensitiveGeo()

  //......................................................................
  void GeometryTestAlg::testCryostat()
  {
    mf::LogVerbatim("GeometryTest")
      << "There are " << geom->Ncryostats() << " cryostats in the detector";

    for (auto const& cryo : geom->Iterate<geo::CryostatGeo>()) {

      {
        mf::LogVerbatim log("GeometryTest");

        log << "\n\tCryostat " << cryo.ID() << " " << cryo.Volume()->GetName()
            << " Dimensions [cm]: " << cryo.Width() << " x " << cryo.Height() << " x "
            << cryo.Length();
        if (fComputeMass) { log << "\n\t\tmass [kg]: " << cryo.Mass(); }
        log << "\n\t\tCryostat boundaries:"
            << "  -x:" << cryo.MinX() << " +x:" << cryo.MaxX() << "  -y:" << cryo.MinY()
            << " +y:" << cryo.MaxY() << "  -z:" << cryo.MinZ() << " +z:" << cryo.MaxZ();
      }

      MF_LOG_DEBUG("GeometryTest") << "\t testing GeometryCore::PoitionToCryostat....";
      geo::CryostatID cid;
      try {
        cid = geom->PositionToCryostatID(cryo.GetCenter());
      }
      catch (cet::exception& e) {
        mf::LogWarning("FailedToLocateCryostat") << "\n exception caught:" << e;
        if (fNonFatalExceptions.count(e.category()) == 0) throw;
      }
      if (cid != cryo.ID()) {
        throw cet::exception("CryostatTest")
          << "Geometry points the middle of cryostat " << std::string(cryo.ID())
          << " to a different one (" << std::string(cid) << ")\n";
      }
      MF_LOG_DEBUG("GeometryTest") << "done";

      MF_LOG_DEBUG("GeometryTest") << "\t Now test the TPCs associated with this cryostat";
      testTPC(cryo.ID());
    }

    return;
  }

  //......................................................................
  unsigned int GeometryTestAlg::testFindWorldVolumes()
  {

    unsigned int nErrors = 0;

    std::set<std::string> volume_names;
    std::vector<TGeoNode const*> nodes;

    // world
    volume_names.insert(geom->GetWorldVolumeName());
    nodes = geom->FindAllVolumes(volume_names);
    {
      mf::LogVerbatim log("GeometryTest");
      log << "Found " << nodes.size() << " world volumes '" << geom->GetWorldVolumeName() << "':";
      for (TGeoNode const* node : nodes) {
        TGeoVolume const* pVolume = node->GetVolume();
        log << "\n - '" << pVolume->GetName() << "' (a " << pVolume->GetShape()->GetName() << ")";
      } // for
    }   // anonymous block
    if (nodes.size() != 1) {
      ++nErrors;
      mf::LogError("GeometryTest") << "Found " << nodes.size() << " world volumes '"
                                   << geom->GetWorldVolumeName() << "! [expecting: one!!]";
    } // if nodes

    return nErrors;
  } // GeometryTestAlg::testFindWorldVolumes()

  unsigned int GeometryTestAlg::testFindCryostatVolumes()
  {

    unsigned int nErrors = 0;

    std::set<std::string> volume_names;
    volume_names.insert(geom->GetWorldVolumeName());
    volume_names.insert("volCryostat");

    std::vector<TGeoNode const*> nodes = geom->FindAllVolumes(volume_names);

    mf::LogVerbatim log("GeometryTest");
    log << "Found " << nodes.size() << " world and cryostat volumes:";
    for (TGeoNode const* node : nodes) {
      TGeoVolume const* pVolume = node->GetVolume();
      log << "\n - '" << pVolume->GetName() << "' (a " << pVolume->GetShape()->GetName() << ")";
    } // for
    if (nodes.size() != (1 + geom->Ncryostats())) {
      ++nErrors;
      mf::LogError("GeometryTest") << "Found " << nodes.size()
                                   << " world and cryostat volumes! "
                                      "[expecting: 1 world and "
                                   << geom->Ncryostats() << " cryostats]";
    } // if nodes

    return nErrors;
  } // GeometryTestAlg::testFindCryostatVolumes()

  unsigned int GeometryTestAlg::testFindTPCvolumePaths()
  {
    unsigned int nErrors = 0;

    // search the full path of all TPCs
    std::set<std::string> volume_names;
    for (auto const& TPC : geom->Iterate<geo::TPCGeo>())
      volume_names.insert(TPC.TotalVolume()->GetName());

    // get the right answer: how many TPCs?
    const unsigned int NTPCs = geom->TotalNTPC();

    std::vector<std::vector<TGeoNode const*>> node_paths = geom->FindAllVolumePaths(volume_names);

    mf::LogVerbatim log("GeometryTest");
    log << "Found " << node_paths.size() << " TPC volumes:";
    for (auto const& path : node_paths) {
      TGeoNode const* node = path.back();
      TGeoVolume const* pVolume = node->GetVolume();
      log << "\n - '" << pVolume->GetName() << "' (a " << pVolume->GetShape()->GetName()
          << ") with " << (path.size() - 1) << " ancestors";
      for (TGeoNode const* pNode : path) {
        TGeoVolume const* pVolume = pNode->GetVolume();
        log << "\n      * '" << pVolume->GetName() << "' (a " << pVolume->GetShape()->GetName()
            << ") with a " << pNode->GetMatrix()->IsA()->GetName() << " that "
            << (pNode->GetMatrix()->IsTranslation() ? "is" : "is not") << " a simple translation";
      } // for node
    }   // for path
    if (node_paths.size() != NTPCs) {
      ++nErrors;
      mf::LogError("GeometryTest") << "Found " << node_paths.size()
                                   << " TPC volumes! "
                                      "[expecting: "
                                   << NTPCs << " TPCs]";
    } // if missed some

    return nErrors;
  } // GeometryTestAlg::testFindTPCvolumePaths()

  void GeometryTestAlg::testFindVolumes()
  {
    /*
     * Finds and checks a selected number of volumes by name:
     * - world volume
     * - cryostat volumes
     * - TPCs (returns full paths)
     */

    unsigned int nErrors = 0;

    nErrors += testFindWorldVolumes();
    nErrors += testFindCryostatVolumes();
    nErrors += testFindTPCvolumePaths();

    if (nErrors != 0) {
      throw cet::exception("FindVolumes")
        << "Collected " << nErrors << " errors during FindAllVolumes() test!\n";
    }

  } // GeometryTestAlg::testFindVolumes()

  //......................................................................
  void GeometryTestAlg::testTPC(geo::CryostatID const& cid)
  {
    geo::CryostatGeo const& cryo = geom->Cryostat(cid);

    mf::LogVerbatim("GeometryTest") << "\tThere are " << cryo.NTPC() << " TPCs in the detector";

    for (size_t t = 0; t < cryo.NTPC(); ++t) {
      geo::TPCID const tpcid(cid, t);
      geo::TPCGeo const& tpc = cryo.TPC(tpcid);
      decltype(auto) activeCenter = tpc.GetActiveVolumeCenter();

      {
        mf::LogVerbatim log{"GeometryTest"};
        log << "\n\t\tTPC " << tpcid << " " << geom->GetLArTPCVolumeName(tpcid) << " has "
            << tpc.Nplanes() << " planes."
            << "\n\t\tTPC location: ( " << tpc.MinX() << " ; " << tpc.MinY() << " ; " << tpc.MinZ()
            << " ) =>  ( " << tpc.MaxX() << " ; " << tpc.MaxY() << " ; " << tpc.MaxZ() << " ) [cm]"
            << "\n\t\tTPC Dimensions (W x H x L, cm): " << tpc.Width() << " ("
            << directionName(tpc.WidthDir()) << ")"
            << " x " << tpc.Height() << " (" << directionName(tpc.HeightDir()) << ")"
            << " x " << tpc.Length() << " (" << directionName(tpc.LengthDir()) << ")"
            << "\n\t\tTPC Active Dimensions: " << 2. * tpc.ActiveHalfWidth() << " x "
            << 2. * tpc.ActiveHalfHeight() << " x " << tpc.ActiveLength() << " around ( "
            << activeCenter.X() << " ; " << activeCenter.Y() << " ; " << activeCenter.Z()
            << " ) cm";
        if (fComputeMass) log << "\n\t\tTPC mass: " << tpc.ActiveMass();
        log << "\n\t\tTPC drift distance: " << tpc.DriftDistance()
            << ", direction: " << tpc.DriftDir();
      }

      for (size_t p = 0; p < tpc.Nplanes(); ++p) {
        geo::PlaneGeo const& plane = tpc.Plane(p);

        // first line indented with two tabs, the others with two more spaces;
        // very verbose (8)
        plane.PrintPlaneInfo(mf::LogVerbatim("GeometryTest") << "\t\t", "\t\t  ", 8);

        mf::LogVerbatim("GeometryTest") << "\t\t  pitch from plane 0 is " << tpc.Plane0Pitch(p);

      } // for plane
      geo::DriftDirection_t dir = tpc.DriftDirection();
      if (dir == geo::kNegX) {
        mf::LogVerbatim("GeometryTest")
          << "\t\tdrift direction is towards negative values: " << tpc.DriftDir();
      }
      else if (dir == geo::kPosX) {
        mf::LogVerbatim("GeometryTest")
          << "\t\tdrift direction is towards positive values: " << tpc.DriftDir();
      }
      else {
        throw cet::exception("UnknownDriftDirection") << "\t\tdrift direction is unknown\n";
      }

      MF_LOG_DEBUG("GeometryTest") << "\t testing PositionToTPC...";
      // pick a position in the middle of the TPC in the world coordinates
      geo::TPCGeo::LocalPoint_t const localLoc{};
      auto const worldLoc = tpc.toWorldCoords(localLoc);

      const unsigned int tpcNo = cryo.PositionToTPCID(worldLoc, 1 + 1.e-4).TPC;

      if (tpcNo != t)
        throw cet::exception("BadTPCLookupFromPosition")
          << "TPC look up returned tpc = " << tpcNo << " should be " << t << "\n";

      MF_LOG_DEBUG("GeometryTest") << "done.";
    } // for TPC
  }

  //......................................................................
  void GeometryTestAlg::testPlaneDirections() const
  {
    /*
     * The test verifies that all the planes in the geometry respect the
     * orientation requirements:
     *
     *   { (wire direction) , (wire coordinate increase), (plane normal) }
     *
     *   { (width direction) , (length direction), (plane normal) }
     *
     * both be a positively defined base.
     */

    lar::util::RealComparisons<double> coordIs(1e-5);

    unsigned int nErrors = 0;
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {
      //
      // check the ( wire ; wire coordinate ; normal) base
      //

      // this funny way declares a reference or not, depending on return type
      decltype(auto) planeNormal = plane.GetNormalDirection();
      decltype(auto) wireCoordDir = plane.GetIncreasingWireDirection();
      decltype(auto) wireDir = plane.GetWireDirection();

      double const wireFrame = wireDir.Cross(wireCoordDir).Dot(planeNormal);

      if (coordIs.nonEqual(wireFrame, 1.0)) {
        ++nErrors;

        mf::LogProblem("GeometryTestAlg")
          << "Plane " << plane.ID() << " has wire direction " << wireDir
          << " wire coordinate direction " << wireCoordDir << " and normal " << planeNormal
          << " , yielding to a non-positive plane-coordinate definition"
          << " (l x w . n = " << wireFrame << ")";
      } // if error

      //
      // check the ( width ; depth ; normal) base
      //

      decltype(auto) widthDir = plane.WidthDir();
      decltype(auto) depthDir = plane.DepthDir();
      double const geoFrame = widthDir.Cross(depthDir).Dot(planeNormal);

      if (coordIs.nonEqual(geoFrame, 1.0)) {
        ++nErrors;

        mf::LogProblem("GeometryTestAlg")
          << "Plane " << plane.ID() << " has width direction " << widthDir << " depth direction "
          << depthDir << " and normal " << planeNormal
          << " , yielding to a non-positive plane-coordinate definition"
          << " (w x d . n = " << geoFrame << ")";
      } // if error

    } // for plane

    if (nErrors > 0) {
      throw cet::exception("GeometryTestAlg")
        << "testPlaneDirections() accumulated " << nErrors << " errors (see messages above)\n";
    }

  } // GeometryTestAlg::testPlaneDirections()

  //......................................................................
  void GeometryTestAlg::testWireOrientations() const
  {
    /*
     * The test verifies that all the wires in the geometry respect the
     * orientation requirement described in geo::WireGeo documentation:
     *
     *   { (wire direction) , (wire coordinate increase), (plane normal) }
     *
     * be a positively defined base.
     *
     */

    unsigned int nErrors = 0;
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      // this funny way declares a reference or not, depending on return type
      decltype(auto) planeNormal = plane.GetNormalDirection();
      decltype(auto) wireCoordDir = plane.GetIncreasingWireDirection();

      unsigned int nWires = plane.Nwires();
      for (unsigned int wireNo = 0; wireNo < nWires; ++wireNo) {

        geo::WireGeo const& wire = plane.Wire(wireNo);

        double positive = wire.Direction().Cross(wireCoordDir).Dot(planeNormal);

        if (positive < 0.5) {
          ++nErrors;

          // detail the problem; details of the plane should be read in the
          // output from other tests
          decltype(auto) wireDir = wire.Direction();
          mf::LogProblem("GeometryTestAlg")
            << "Wire " << plane.ID() << " W: " << wireNo << " has direction ( " << wireDir
            << " , yielding to a non-positive plane-coordinate definition"
            << " (l x w . n = " << positive << ")";
        } // if error

      } // for wire

    } // for plane

    if (nErrors > 0) {
      throw cet::exception("GeometryTestAlg")
        << "testWireOrientations() accumulated " << nErrors << " errors (see messages above)\n";
    }

  } // GeometryTestAlg::testWireOrientations()

  //......................................................................
  void GeometryTestAlg::testWireCoordFromPlane() const
  {

    //
    // For each wire:
    //
    // * picks points lying on the planes including a wire and the normal to the
    //   wire plane (which have all the same wire coordinate)
    //
    // * tests that the coordinates are as expected (wire number times pitch)
    //

    unsigned int nErrors = 0;
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      auto const nWires = plane.Nwires();
      auto const wirePitch = plane.WirePitch();

      double const driftDistance = geom->TPC(plane.ID()).DriftDistance();

      decltype(auto) planeNormal = plane.GetNormalDirection();

      for (geo::WireID::WireID_t wireNo = 0; wireNo < nWires; ++wireNo) {

        geo::WireGeo const& wire = plane.Wire(wireNo);

        double const expected = wireNo * wirePitch;

        // sample 7 points on wire
        constexpr int shifts = 3;
        double const step = wire.HalfL() / (std::abs(shifts) + 1);
        for (int iOfs = -shifts; iOfs <= shifts; ++iOfs) {

          double const offset = iOfs * step;

          auto const basePoint = wire.GetPositionFromCenter(offset);

          // at 4 different distances from the plane
          constexpr int quotas = 4;
          double const jump = driftDistance / (std::abs(quotas) + 1);
          for (int iQuota = 0; iQuota < quotas; ++iQuota) {

            // translate the point along the normal to the plane;
            // this should not change the result
            auto const point = basePoint + iQuota * jump * planeNormal;

            double const distance = plane.PlaneCoordinate(geo::vect::toPoint(point));

            if (std::abs(distance - expected) > 1e-4) {
              mf::LogProblem("GeometryTestAlg")
                << "Point " << point << "  (offset: " << iOfs << "x" << step << ", at " << iQuota
                << "x" << jump << " from plane) is reported to be " << distance
                << " cm far from wire " << plane.ID() << " W: " << wireNo << " (" << expected
                << " expected)";
              ++nErrors;
            } // if unexpected

          } // for quotas

        } // for iOfs

      } // for wires

    } // for planes

    if (nErrors > 0) {
      throw cet::exception("GeometryTestAlg")
        << "testWireCoordFromPlane() accumulated " << nErrors << " errors (see messages above)\n";
    }

  } // GeometryTestAlg::testWireCoordFromPlane()

  //......................................................................
  void GeometryTestAlg::testParallelWires() const
  {

    //
    // checks that all the wires in the same plane are parallel
    //
    auto const vectorIs = lar::util::makeVector3DComparison(1e-8);

    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      decltype(auto) genDir = plane.GetWireDirection();

      geo::WireID::WireID_t wireNo = 0;
      for (geo::WireGeo const& wire : plane.IterateWires()) {

        geo::WireID const wireID(plane.ID(), wireNo++);

        decltype(auto) wireDir = wire.Direction();

        if (vectorIs.nonEqual(wireDir, genDir)) {
          throw cet::exception("ParallelWires")
            << "Wire " << std::string(wireID) << " has direction " << wireDir
            << ", not parallel to wire direction " << genDir << " according to the plane "
            << std::string(plane.ID()) << "\n";
        }

      } // for wires in plane

    } // for planes

  } // GeometryTestAlg::testParallelWires()

  //......................................................................
  void GeometryTestAlg::testPlanePointDecomposition() const
  {

    //
    // For each plane:
    //
    // 1) create a plane point with arbitrary distance from the plane,
    //    wire coordinate multiple (N) of wire pitch, and wire direction
    //    coordinate 0 or half a wire length in either directions
    //
    // 2) compose into a 3D vector, and verify that the nearest wire is the one
    //    expected (N)
    //
    // 3) decompose back the 3D vector, and verify that the result matches the
    //    starting decomposition
    //
    // 4) also verify singly PointProjection() and DistanceFromPlane()
    //
    // 5) verify DriftPoint()
    //
    //

    lar::util::RealComparisons<double> coordIs(1e-5);
    auto vectorIs = lar::util::makeVector3DComparison(coordIs);

    unsigned int nErrors = 0;
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      auto const& planeNorm = plane.GetNormalDirection();
      auto const& wirePitch = plane.WirePitch();
      auto const& refPoint = plane.ProjectionReferencePoint();

      unsigned int const lastWire = plane.Nwires() - 1;
      geo::WireID::WireID_t wireNo = 0;
      for (geo::WireGeo const& wire : plane.IterateWires()) {

        geo::WireID const wireID(plane.ID(), wireNo++);

        constexpr double distance = 5.0; // 5 cm from the plane

        auto const wireDirStep = wire.HalfL() / 2.0; // quarter of the length

        for (int iWDStep = -1; iWDStep <= 1; ++iWDStep) {

          //
          // prepare expectation
          //
          auto const wireDirOffset = iWDStep * wireDirStep;

          auto const expectedPoint =
            wire.GetCenter() + wireDirOffset * wire.Direction() + distance * planeNorm;

          auto const expectedWireCoord = wirePitch * wireID.Wire;
          auto const expectedWireDirCoord =
            wireDirOffset + wire.Direction().Dot(wire.GetCenter() - refPoint);

          geo::PlaneGeo::WireCoordProjection_t const expectedProj(expectedWireDirCoord,
                                                                  expectedWireCoord);

          //
          // composition
          //
          auto point = plane.ComposePoint(distance, expectedProj);

          if (vectorIs.nonEqual(point, expectedPoint)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] ComposePoint(): "
              << "Point with projection " << expectedProj << " and distance from plane " << distance
              << " was reported as " << point << " while it is expected to be at " << expectedPoint
              << " (on wire " << std::string(wireID) << ")";
          } // if wrong point

          //
          // decomposition
          //
          auto const decomp = plane.DecomposePoint(point);
          if (coordIs.nonEqual(decomp.distance, distance)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] DecomposePoint(): "
              << "Point " << point << " (on wire " << std::string(wireID) << ")"
              << " is reported to have distance from plane " << decomp.distance << " cm, while "
              << distance << " is expected";
          } // if wrong distance
          if (coordIs.nonEqual(decomp.projection.X(), expectedWireDirCoord)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] DecomposePoint(): "
              << "Point " << point << " (on wire " << std::string(wireID) << ")"
              << " is reported to have wire direction coordinate " << decomp.projection.X()
              << " cm, while " << expectedWireDirCoord << " is expected";
          } // if wrong coordinate along the wire
          if (coordIs.nonEqual(decomp.projection.Y(), expectedWireCoord)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] DecomposePoint(): "
              << "Point " << point << " (on wire " << std::string(wireID) << ")"
              << " is reported to have wire coordinate " << decomp.projection.Y() << " cm, while "
              << expectedWireCoord << " is expected";
          } // if wrong wire coordinate

          //
          // projection
          //
          auto const proj = plane.PointProjection(point);
          if (coordIs.nonEqual(proj.X(), expectedWireDirCoord)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] PointProjection(): "
              << "Point " << point << " (on wire " << std::string(wireID) << ")"
              << " is reported to have wire direction coordinate " << proj.X() << " cm, while "
              << expectedWireDirCoord << " is expected";
          } // if wrong wire coordinate
          if (coordIs.nonEqual(proj.Y(), expectedWireCoord)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] PointProjection(): "
              << "Point " << point << " (on wire " << std::string(wireID) << ")"
              << " is reported to have wire coordinate " << proj.Y() << " cm, while "
              << expectedWireCoord << " is expected";
          } // if wrong wire coordinate

          //
          // distance
          //
          auto const dist = plane.DistanceFromPlane(point);
          if (coordIs.nonEqual(dist, distance)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] DistanceFromPlane(): "
              << "Point " << point << " (on wire " << std::string(wireID) << ")"
              << " is reported to have distance " << dist << " cm from plane " << plane.ID()
              << ", while " << distance << " is expected";
          } // if wrong distance

          //
          // drift
          //
          // BUG the double brace syntax is required to work around clang bug 21629
          // (https://bugs.llvm.org/show_bug.cgi?id=21629)
          std::array<double, 3> drifts{{-distance, distance, 2. * distance}};
          for (double drift : drifts) {

            // DriftPoint() moves the point in the drift direction,
            // which is opposite to the plane normal:
            auto const expectedDistance = distance - drift;

            //
            // drift it by a known value
            //
            auto point = expectedPoint;
            plane.DriftPoint(point, drift);

            //
            // check the new distance
            //
            auto dist = plane.DistanceFromPlane(point);
            if (coordIs.nonEqual(dist, expectedDistance)) {
              ++nErrors;
              mf::LogProblem("GeometryTestAlg")
                << "[testPlanePointDecomposition] DriftPoint(): "
                << "Point " << expectedPoint << " (distant " << distance << " cm from the plane)"
                << " (on wire " << std::string(wireID) << ")"
                << " drifted by " << drift << " became " << point
                << " which has distance from plane " << dist << " cm, while " << expectedDistance
                << " was expected";
            } // if wrong distance

          } // for drifts

          //
          // containment
          //
          // skip this test for the first and last wire, since the containment
          // area might well exclude the whole wire if it's short enough
          bool const isExtremeWire = (wireID.Wire == 0) || (wireID.Wire == lastWire);
          if (!isExtremeWire && !plane.isProjectionOnPlane(expectedPoint)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] isProjectionOnPlane(): "
              << "Point " << expectedPoint << " (center of " << wireID
              << ") is not believed to be on the plane, but it should.";
          }

        } // for different wire direction steps

      } // for wires in plane

    } // for planes

    if (nErrors > 0) {
      throw cet::exception("GeometryTestAlg") << "testPlanePointDecomposition() accumulated "
                                              << nErrors << " errors (see messages above)\n";
    }

  } // GeometryTestAlg::testPlanePointDecomposition()

  //......................................................................
  void GeometryTestAlg::testWireCoordAngle() const
  {
    /*
     * Tests that the angle PhiZ() actually points to the next wire.
     *
     * The test, for each plane, performs the following:
     * - pick the middle wire, verify that we can get the expected wire
     *   coordinate for its centre
     * - move one wire pitch away from the centre in the direction determined
     *   by PhiZ(), verify that the coordinate increases by 1
     */

    for (auto const& planeid : geom->Iterate<geo::PlaneID>()) {

      geo::PlaneGeo const& plane = geom->Plane(planeid);

      // define the wires to work with
      const unsigned int nWires = plane.Nwires();

      geo::WireID middle_wire_id(planeid, nWires / 2);
      geo::WireID next_wire_id(planeid, nWires / 2 + 1);

      if (next_wire_id.Wire >= nWires) {
        throw cet::exception("WeirdGeometry")
          << "Plane " << std::string(planeid) << " has only " << nWires << " wires?!?\n";
      }

      geo::WireGeo const& middle_wire = geom->Wire(middle_wire_id);
      decltype(auto) middle_wire_center = middle_wire.GetCenter();
      MF_LOG_TRACE("GeometryTest")
        << "Center of " << middle_wire_id << " at " << middle_wire_center;

      // cross check: we can find the middle wire
      const double middle_coord = geom->WireCoordinate(middle_wire_center, planeid);

      if (std::abs(middle_coord - double(middle_wire_id.Wire)) > 2e-3) {
        throw cet::exception("WireCoordAngle")
          << "Center of " << std::string(middle_wire_id) << " at (" << middle_wire_center.X()
          << "; " << middle_wire_center.Y() << "; " << middle_wire_center.Z()
          << ") has wire coordinate " << middle_coord << " (" << middle_wire_id.Wire
          << " expected)\n";
      } // if

      // the check: this coordinate should lie on the next wire
      const double pitch = plane.WirePitch();
      decltype(auto) wireCoordDir = plane.GetIncreasingWireDirection();

      MF_LOG_TRACE("GeometryTest")
        << "  pitch: " << pitch << " wire coord dir: cos(phi_z): " << wireCoordDir;

      auto on_next_wire = middle_wire_center + pitch * wireCoordDir;

      const double next_coord = geom->WireCoordinate(on_next_wire, planeid);

      if (std::abs(next_coord - double(next_wire_id.Wire)) > 2e-3) {
        std::cerr << "  pitch: " << pitch << " wire coord dir: " << wireCoordDir
                  << "\n  start wire: " << middle_wire_id << " centered at " << middle_wire_center
                  << "\n  querying point " << on_next_wire << std::endl;
        throw cet::exception("WireCoordAngle")
          << "Position " << on_next_wire << " is expected to be on wire "
          << std::string(next_wire_id) << " but it has wire coordinate " << next_coord << "\n";
      } // if

    } // for

  } // GeometryTestAlg::testWireCoordAngle()

  //......................................................................
  void GeometryTestAlg::testChannelToROP() const
  {

    // test that an invalid channel yields an invalid ROP
    try {
      readout::ROPID invalidROP = geom->ChannelToROP(raw::InvalidChannelID);
      if (invalidROP.isValid) {
        throw cet::exception("testChannelToROP")
          << "ROP from an invalid channel (" << raw::InvalidChannelID << ") is "
          << std::string(invalidROP) << " !?\n";
      } // if invalid rop is not invalid
    }
    catch (cet::exception const& e) {
      mf::LogWarning("testChannelToROP")
        << "Non-compilant ChannelToROP() throws on invalid channel.";
    }

    // for each channel, test that its ROP contains it;
    // we assume each ROP contains contiguous channel IDs
    for (raw::ChannelID_t channel = 0; channel < geom->Nchannels(); ++channel) {

      readout::ROPID const ropid = geom->ChannelToROP(channel);

      auto const firstChannel = geom->FirstChannelInROP(ropid);
      auto const lastChannel = firstChannel + geom->Nchannels(ropid);

      if ((channel < firstChannel) || (channel >= lastChannel)) {
        throw cet::exception("testChannelToROP")
          << "Channel " << channel << " comes from ROP " << std::string(ropid)
          << ", which contains only channels from " << firstChannel << " to " << lastChannel
          << " (excluded)\n";
      } // if

    } // for channel

  } // GeometryTestAlg::testChannelToROP()

  //......................................................................
  void GeometryTestAlg::testChannelToWire() const
  {
    using std::begin;
    using std::end;

    geo::PlaneID lastPlane; // starts invalid
    geo::View_t planeView = geo::kUnknown;
    geo::SigType_t planeSigType = geo::kMysteryType;

    for (auto const& testWireID : geom->Iterate<geo::WireID>()) {

      raw::ChannelID_t channel = geom->PlaneWireToChannel(testWireID);

      if (!raw::isValidChannelID(channel)) {
        throw cet::exception("BadChannelLookup")
          << "Invalid channel returned for wire " << std::string(testWireID) << "\n";
      }

      auto const wireIDs = geom->ChannelToWire(channel);

      if (wireIDs.empty()) {
        throw cet::exception("BadChannelLookup")
          << "List of wires for channel #" << channel << " is empty;"
          << " it should have contained at least " << std::string(testWireID) << "\n";
      }

      if (std::count(begin(wireIDs), end(wireIDs), testWireID) != 1) {
        mf::LogError log("GeometryTest");
        log << wireIDs.size() << " wire IDs associated with channel #" << channel << ":";
        for (auto const& wid : wireIDs)
          log << "\n  " << wid;
        throw cet::exception("BadChannelLookup")
          << "Wire " << std::string(testWireID) << " is on channel #" << channel
          << " but ChannelToWire() does not map the channel to that wire\n";
      }

      // currently (LArSoft 6.12) signal type from channel and from plane use
      // the same underlying code, so the following test is not very valuable
      auto const channelSigType = geom->SignalType(channel);
      if (channelSigType != geom->SignalType(testWireID.planeID())) {
        throw cet::exception("BadChannelLookup")
          << "Geometry service claims channel #" << channel << " to be of type " << channelSigType
          << " but that the plane of " << std::string(testWireID) << " is of type "
          << geom->SignalType(testWireID.planeID()) << "\n";
      }

      auto const channelView = geom->View(channel);
      if (channelView != geom->Plane(testWireID).View()) {
        throw cet::exception("BadChannelLookup")
          << "Geometry service claims channel #" << channel << " should be on view "
          << geom->View(channel) << " but the plane of " << std::string(testWireID)
          << " claims to be on view " << geom->Plane(testWireID).View() << "\n";
      }

      // check that all the channels on the same plane are consistent
      // (sort of: it's more like "all contiguous wires in the same plane")
      if (testWireID.planeID() == lastPlane) {
        if (planeView != channelView) {
          throw cet::exception("BadChannelLookup")
            << "Geometry service claims channel #" << channel << " is on view " << channelView
            << " but its plane " << std::string(lastPlane) << " claims to be on view " << planeView
            << "\n";
        }
        if (planeSigType != channelSigType) {
          throw cet::exception("BadChannelLookup")
            << "Geometry service claims channel #" << channel << " is if type " << channelSigType
            << " but its plane " << std::string(lastPlane) << " to be of type " << planeSigType
            << "\n";
        }
      }
      else {
        lastPlane = testWireID.planeID();
        planeView = channelView;
        planeSigType = channelSigType;
      }

    } // for wires in detector

  } // GeometryTestAlg::testChannelToWire()

  //......................................................................
  void GeometryTestAlg::testFindPlaneCenters()
  {
    geo::PlaneGeo::LocalPoint_t const origin{};
    for (auto const& plane : geom->Iterate<PlaneGeo>(CryostatID{0})) {
      auto const xyzW = plane.toWorldCoords(origin);
      mf::LogVerbatim("GeometryTest") << "\n\tplane " << plane.ID() << " is centered at (x,y,z) = ("
                                      << xyzW.X() << "," << xyzW.Y() << "," << xyzW.Z() << ")";
    }
  }

  //......................................................................
  void GeometryTestAlg::testPlaneProjectionReference() const
  {

    //
    // Check the definition of the projection reference
    //
    lar::util::RealComparisons<double> coordIs(1e-8);
    unsigned int nErrors = 0;
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      auto const reference = plane.ProjectionReferencePoint();

      auto decomp = plane.DecomposePoint(reference);

      if (coordIs.nonZero(decomp.distance)) {
        MF_LOG_ERROR("GeometryTest") << "Plane " << plane.ID() << " reference point " << reference
                                     << " has distance " << decomp.distance << " cm (should be 0)";
        ++nErrors;
      }

      if (coordIs.nonZero(decomp.projection.X()) || coordIs.nonZero(decomp.projection.Y())) {
        MF_LOG_ERROR("GeometryTest")
          << "Plane " << plane.ID() << " reference point " << reference << " has projection ( "
          << decomp.projection.X() << " ; " << decomp.projection.Y() << " ) cm (should be (0;0) )";
        ++nErrors;
      }

    } // for planes
    if (nErrors > 0) {
      throw cet::exception("GeoTestPlaneProjectionReference")
        << "Accumulated " << nErrors << " errors (see messages above)\n";
    }
  } // GeometryTestAlg::testPlaneProjectionReference()

  //......................................................................
  void GeometryTestAlg::testPlanePointDecompositionFrame() const
  {

    //
    // For each plane:
    //
    // 1) create a plane point with arbitrary distance from the plane,
    //    width and depth coordinates all across the plane
    //
    // 2) compose into a 3D vector
    //
    // 3) decompose back the 3D vector, and verify that the result matches the
    //    starting decomposition
    //
    // 4) also verify singly PointProjection() and DistanceFromPlane()
    //
    // 5) verify DriftPoint()
    //
    //

    lar::util::RealComparisons<double> coordIs(1e-5);
    auto vectorIs = lar::util::makeVector3DComparison(coordIs);

    // steps on each side of the center, within the plane:
    constexpr int steps = 5;
    static_assert(steps > 0, "Steps must be a positive number.");
    // how many width units we go far from the center (1 means stay inside)
    constexpr int nOutsides = 1;

    unsigned int nErrors = 0;
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      auto const& planeNorm = plane.GetNormalDirection();
      auto const& widthDir = plane.WidthDir();
      auto const& depthDir = plane.DepthDir();
      auto const& refPoint = plane.GetCenter();

      double const halfWidth = plane.Width() / 2;
      double const halfDepth = plane.Depth() / 2;
      double const dW_2 = halfWidth / (steps * 2); // half width step
      double const dD_2 = halfDepth / (steps * 2); // half depth step

      for (int iW = -nOutsides * steps; iW <= +nOutsides * steps; ++iW) {

        double const expected_w = dW_2 * (iW * 2 + 1); // width coordinate
        bool const inWidth = std::abs(expected_w) <= halfWidth;

        for (int iD = -nOutsides * steps; iD <= +nOutsides * steps; ++iD) {

          double const expected_d = dD_2 * (iD * 2 + 1); // depth coordinate
          bool const inDepth = std::abs(expected_d) <= halfDepth;

          constexpr double distance = 5.0; // we might test this too...

          //
          // prepare expectation
          //
          auto const expectedPoint =
            refPoint + expected_w * widthDir + expected_d * depthDir + distance * planeNorm;

          geo::PlaneGeo::WidthDepthProjection_t const expectedProj(expected_w, expected_d);

          //
          // composition
          //
          auto point = plane.ComposePoint(distance, expectedProj);

          if (vectorIs.nonEqual(point, expectedPoint)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecompositionFrame] ComposePoint(): "
              << "Point with projection " << expectedProj << " (width: " << expected_w
              << ", depth: " << expected_d << ") and distance from plane " << distance
              << " was reported as " << point << " while it is expected to be at " << expectedPoint;
          } // if wrong point

          //
          // decomposition
          //
          auto const decomp = plane.DecomposePointWidthDepth(point);
          if (coordIs.nonEqual(decomp.distance, distance)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] DecomposePointWidthDepth(): "
              << "Point " << point << " (width: " << expected_w << ", depth: " << expected_d
              << ") is reported to have distance from plane " << decomp.distance << " cm, while "
              << distance << " is expected";
          } // if wrong distance
          if (coordIs.nonEqual(decomp.projection.X(), expected_w)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] DecomposePointWidthDepth(): "
              << "Point " << point << " (width: " << expected_w << ", depth: " << expected_d
              << ") is reported to have width direction coordinate " << decomp.projection.X()
              << " cm, while " << expected_w << " is expected";
          } // if wrong coordinate along the wire
          if (coordIs.nonEqual(decomp.projection.Y(), expected_d)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] DecomposePointWidthDepth(): "
              << "Point " << point << " (width: " << expected_w << ", depth: " << expected_d
              << ") is reported to have depth direction coordinate " << decomp.projection.Y()
              << " cm, while " << expected_d << " is expected";
          } // if wrong wire coordinate

          //
          // projection
          //
          auto const proj = plane.PointWidthDepthProjection(point);
          if (coordIs.nonEqual(proj.X(), expected_w)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] PointProjection(): "
              << "Point " << point << " (width: " << expected_w << ", depth: " << expected_d
              << ") is reported to have width direction coordinate " << proj.X() << " cm, while "
              << expected_w << " is expected";
          } // if wrong wire coordinate
          if (coordIs.nonEqual(proj.Y(), expected_d)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] PointProjectionWidthDepth(): "
              << "Point " << point << " (width: " << expected_w << ", depth: " << expected_d
              << ") is reported to have wire coordinate " << proj.Y() << " cm, while " << expected_d
              << " is expected";
          } // if wrong wire coordinate

          //
          // distance
          //
          auto const dist = plane.DistanceFromPlane(point);
          if (coordIs.nonEqual(dist, distance)) {
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecomposition] DistanceFromPlane(): "
              << "Point " << point << " (width: " << expected_w << ", depth: " << expected_d
              << ") is reported to have distance from plane " << dist << " cm, while " << distance
              << " is expected";
          } // if wrong distance

          //
          // drift
          //
          // BUG the double brace syntax is required to work around clang bug 21629
          // (https://bugs.llvm.org/show_bug.cgi?id=21629)
          std::array<double, 3> drifts{{-distance, distance, 2. * distance}};
          for (double drift : drifts) {

            // DriftPoint() moves the point in the drift direction,
            // which is opposite to the plane normal:
            auto const expectedDistance = distance - drift;

            //
            // drift it by a known value
            //
            auto point = expectedPoint;
            plane.DriftPoint(point, drift);

            //
            // check the new distance
            //
            auto dist = plane.DistanceFromPlane(point);
            if (coordIs.nonEqual(dist, expectedDistance)) {
              ++nErrors;
              mf::LogProblem("GeometryTestAlg")
                << "[testPlanePointDecomposition] DriftPoint(): "
                << "Point " << expectedPoint << " (distant " << distance << " cm from the plane)"
                << " (width: " << expected_w << ", depth: " << expected_d << ") drifted by "
                << drift << " became " << point << " which has distance from plane " << dist
                << " cm, while " << expectedDistance << " was expected";
            } // if wrong distance

          } // for drifts

          //
          // containment
          //
          const bool expected_onPlane = inWidth && inDepth;
          const bool onPlane = plane.isProjectionOnPlane(expectedPoint);
          if (expected_onPlane != onPlane) {
            // always
            ++nErrors;
            mf::LogProblem("GeometryTestAlg")
              << "[testPlanePointDecompositionFrame] isProjectionOnPlane(): "
              << "Point " << expectedPoint << " (width: " << expected_w << ", depth: " << expected_d
              << ") is" << (onPlane ? "" : " not") << " believed to be on plane " << plane.ID()
              << ", but it should" << (expected_onPlane ? "" : " not be") << ".";
          }

        } // for different wire direction steps

      } // for wires in plane

    } // for planes

    if (nErrors > 0) {
      throw cet::exception("GeometryTestAlg") << "testPlanePointDecomposition() accumulated "
                                              << nErrors << " errors (see messages above)\n";
    }

  } // GeometryTestAlg::testPlanePointDecompositionFrame()

  //......................................................................
  void GeometryTestAlg::testPlaneProjectionOnFrame() const
  {

    //
    // Tests:
    //
    // 1. containment (isProjectionOnPlane())
    //
    //
    // 2. capping by closest point
    //

    lar::util::RealComparisons<double> coordIs(1e-5);
    auto vectorIs = lar::util::makeVector3DComparison(coordIs);
    auto const& vector2Dis = vectorIs.comp2D();

    // steps on each side of the center, within the plane:
    constexpr int steps = 5;
    static_assert(steps > 0, "Steps must be a positive number.");
    // how many width units we go far from the center (1 means stay inside)
    constexpr int nOutsides = 2;

    unsigned int nErrors = 0;
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      double const halfWidth = plane.Width() / 2;
      double const halfDepth = plane.Depth() / 2;
      double const dW_2 = halfWidth / (steps * 2); // half width step
      double const dD_2 = halfDepth / (steps * 2); // half depth step

      for (int iW = -nOutsides * steps; iW <= +nOutsides * steps; ++iW) {

        double const expected_w = dW_2 * (iW * 2 + 1); // width coordinate
        bool const inWidth = std::abs(expected_w) <= halfWidth;

        for (int iD = -nOutsides * steps; iD <= +nOutsides * steps; ++iD) {

          double const expected_d = dD_2 * (iD * 2 + 1); // depth coordinate
          bool const inDepth = std::abs(expected_d) <= halfDepth;

          for (double distance : {-30., 0.0, +30.0}) {

            //
            // prepare expectation
            //
            // definition of the test point
            geo::PlaneGeo::WidthDepthProjection_t const expected_proj(expected_w, expected_d);

            auto const expected_point = plane.ComposePoint(distance, expected_proj);

            //
            // 1. Containment test
            //
            const bool expected_onPlane = inWidth && inDepth;
            const bool onPlane = plane.isProjectionOnPlane(expected_point);
            if (expected_onPlane != onPlane) {
              ++nErrors;
              mf::LogProblem("GeometryTestAlg")
                << "[testPlaneProjectionOnFrame] isProjectionOnPlane(): "
                << "Point " << expected_point << " (width: " << expected_w
                << ", depth: " << expected_d << ") is" << (onPlane ? "" : " not")
                << " believed to be on plane " << plane.ID() << ", but it should"
                << (expected_onPlane ? "" : " not be") << ".";
            }

            //
            // 2. capping by closest point
            //
            // 2.1. capping projection
            //
            geo::PlaneGeo::WidthDepthProjection_t expected_movedProjection(
              (inWidth          ? expected_w :
               (expected_w < 0) ? -halfWidth :
                                  +halfWidth),
              (inDepth          ? expected_d :
               (expected_d < 0) ? -halfDepth :
                                  +halfDepth));
            auto movedProjection = plane.MoveProjectionToPlane(expected_proj);
            if (vector2Dis.nonEqual(movedProjection, expected_movedProjection)) {
              ++nErrors;
              mf::LogProblem("GeometryTestAlg")
                << "[testPlaneProjectionOnFrame] moveProjectionOnPlane():"
                << "Projection of ooint " << expected_point << " (width: " << expected_w
                << ", depth: " << expected_d << ") (" << (onPlane ? "on" : "off") << "-plane "
                << plane.ID() << ") was moved to " << movedProjection
                << " while it should have moved to " << expected_movedProjection << ".";
            }

            //
            // 2.2. capping point
            //
            auto expected_movedPoint = plane.ComposePoint(distance, expected_movedProjection);

            auto movedPoint = plane.MovePointOverPlane(expected_point);
            if (vectorIs.nonEqual(movedPoint, expected_movedPoint)) {
              ++nErrors;
              mf::LogProblem("GeometryTestAlg")
                << "[testPlaneProjectionOnFrame] movePointOnPlane(): "
                << "Point " << expected_point << " (width: " << expected_w
                << ", depth: " << expected_d << ") (" << (onPlane ? "on" : "off") << "-plane "
                << plane.ID() << ") was moved to " << movedPoint
                << " while it should have moved to " << expected_movedPoint << ".";
            }

          } // for distance

        } // for depth step

      } // for width step

    } // for planes

    if (nErrors > 0) {
      throw cet::exception("GeoTestPlaneProjection")
        << "Accumulated " << nErrors << " errors (see messages above)\n";
    }

  } // testPlaneProjectionOnFrame()

  //......................................................................
  void GeometryTestAlg::testPlaneProjection() const
  {

    //
    // Check the definition of the reference
    //

    testPlaneProjectionReference();

    //
    // Check the projections and point composition in the plane frame reference
    //
    testPlanePointDecompositionFrame();

    //
    // Check containment
    //
    testPlaneProjectionOnFrame();

  } // GeometryTestAlg::testPlaneProjection()

  //......................................................................
  void GeometryTestAlg::testStandardWirePos()
  {
    for (auto const& cryostat : geom->Iterate<geo::CryostatGeo>()) {
      for (size_t t = 0; t < cryostat.NTPC(); ++t) {
        const geo::TPCGeo* tpc = &cryostat.TPC(t);

        for (size_t i = 0; i < tpc->Nplanes(); ++i) {
          const geo::PlaneGeo* plane = &tpc->Plane(i);

          for (size_t j = 1; j < plane->Nwires(); ++j) {

            geo::WireGeo const& wire = plane->Wire(j);
            geo::WireGeo const& wireprev = plane->Wire(j - 1);

            auto const xyz = wire.GetCenter();
            auto const xyzprev = wireprev.GetCenter();

            // wires increase in +z order
            if (xyz.Z() < xyzprev.Z())
              throw cet::exception("WireOrderProblem")
                << "\n\twires do not increase in +z order in"
                << "Cryostat " << cryostat.ID() << ", TPC " << t << ", Plane " << i << ";  at wire "
                << j << "\n";

          } // end loop over wires
        }   // end loop over planes
      }     // end loop over tpcs
    }       // end loop over cryostats
  }

  //......................................................................
  void GeometryTestAlg::testAPAWirePos()
  {
    geo::TPCGeo::LocalPoint_t const origin{};
    for (auto const& cryostat : geom->Iterate<geo::CryostatGeo>()) {
      for (size_t t = 0; t < cryostat.NTPC(); ++t) {
        const geo::TPCGeo& tpc = cryostat.TPC(t);
        auto const tpcworld = tpc.toWorldCoords(origin);

        for (size_t i = 0; i < tpc.Nplanes(); ++i) {
          const geo::PlaneGeo* plane = &tpc.Plane(i);

          for (size_t j = 1; j < plane->Nwires(); ++j) {
            geo::WireGeo const& wire = plane->Wire(j);
            geo::WireGeo const& wireprev = plane->Wire(j - 1);

            auto const xyz = wire.GetCenter();
            auto const xyzprev = wireprev.GetCenter();

            // top TPC wires increase in -y
            if (tpcworld.X() > 0 && xyz.Y() > xyzprev.Y())
              throw cet::exception("WireOrderProblem")
                << "\n\ttop TPC wires do not increase in -y order in"
                << "Cryostat " << cryostat.ID() << ", TPC " << t << ", Plane " << i << ";  at wire "
                << j << "\n";
            // bottom TPC wires increase in +y
            if (tpcworld.Y() < 0 && xyz.Y() < xyzprev.Y())
              throw cet::exception("WireOrderProblem")
                << "\n\tbottom TPC wires do not increase in +y order in"
                << "Cryostat " << cryostat.ID() << ", TPC " << t << ", Plane " << i << ";  at wire "
                << j << "\n";
          } // end loop over wires
        }   // end loop over planes
      }     // end loop over tpcs
    }       // end loop over cryostats
  }

  //......................................................................
  inline std::array<double, 3> GeometryTestAlg::GetIncreasingWireDirection(
    const geo::PlaneGeo& plane)
  {
    auto IncreasingWireDir = plane.GetIncreasingWireDirection();
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    return {{IncreasingWireDir.X(), IncreasingWireDir.Y(), IncreasingWireDir.Z()}};
  } // GeometryTestAlg::GetIncreasingWireDirection()

  //......................................................................
  void GeometryTestAlg::testNearestWire()
  {
    // Even if you comment it out, please leave the TStopWatch code
    // in this code for additional testing. The NearestChannel routine
    // is the most frequently called in the simulation, so its execution time
    // is an important component of LArSoft's speed.
    TStopwatch stopWatch;
    stopWatch.Start();

    bool bTestWireCoordinate = true;

    // get a wire and find its center
    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      geo::PlaneID const& planeID = plane.ID();
      const unsigned int NWires = plane.Nwires();

      decltype(auto) IncreasingWireDir = plane.GetIncreasingWireDirection();

      MF_LOG_DEBUG("GeoTestWireCoordinate")
        << "The direction of increasing wires for plane " << planeID
        << " (theta=" << plane.Wire(0).ThetaZ() << " pitch=" << plane.WirePitch()
        << " orientation=" << (plane.Orientation() == geo::kHorizontal ? "H" : "V")
        << (plane.WireIDincreasesWithZ() ? "+" : "-") << ") is " << IncreasingWireDir;

      geo::WireID::WireID_t w = 0;
      for (geo::WireGeo const& wire : plane.IterateWires()) {

        geo::WireID wireID(planeID, w++);

        decltype(auto) wire_center = wire.GetCenter();

        uint32_t nearest = 0;
        std::vector<geo::WireID> wireIDs;

        try {
          nearest = geom->NearestChannel(wire_center, planeID);
        }
        catch (cet::exception& e) {
          mf::LogWarning("GeoTestCaughtException") << e;
          if (fNonFatalExceptions.count(e.category()) == 0) throw;
        }

        try {
          wireIDs = geom->ChannelToWire(nearest);

          if (wireIDs.empty()) {
            throw cet::exception("BadPositionToChannel")
              << "test point is at " << wire_center << " nearest channel is " << nearest << " for "
              << std::string(wireID) << "\n";
          }
        }
        catch (cet::exception& e) {
          mf::LogWarning("GeoTestCaughtException") << e;
          if (fNonFatalExceptions.count(e.category()) == 0) throw;
        }

        if (std::find(wireIDs.begin(), wireIDs.end(), wireID) == wireIDs.end()) {
          throw cet::exception("BadPositionToChannel")
            << "Current wire " << std::string(wireID) << " has a world position at " << wire_center
            << "\nNearestWire for this position is "
            << geom->NearestWireID(wire_center, planeID).Wire << "\nNearestChannel is " << nearest
            << "\nShould be channel " << geom->PlaneWireToChannel(wireID);
        } // if good lookup fails

        // nearest wire, integral and floating point
        try {
          // The test consists in sampling NStep (=5) points between the current
          // wire and the previous/next, following the normal to the wire.
          // We expect WireCoordinate() to reflect the same shift.

          // using absolute value just in case (what happens if w1 > w2?)

          const double pitch = std::abs(geom->WirePitch(planeID));

          auto const step = pitch * IncreasingWireDir;

          constexpr int NSteps = 5; // odd value avoids testing half-way
          for (int i = -NSteps; i <= +NSteps; ++i) {
            // we move away by this fraction of wire:
            const double f = NSteps ? (double(i) / NSteps) : 0.0;

            // these are the actual shifts on the positive directions y and z
            auto const delta = f * step;
            auto const wire_shifted = wire_center + delta;

            // we expect this wire number
            const double expected_wire = wireID.Wire + f;

            double wire_from_wc = 0;
            if (bTestWireCoordinate) {
              try {
                wire_from_wc = geom->WireCoordinate(wire_shifted, planeID);
              }
              catch (cet::exception& e) {
                if (hasCategory(e, "NotImplemented")) {
                  MF_LOG_ERROR("GeoTestWireCoordinate")
                    << "WireCoordinate() is not implemented for your ChannelMap;"
                       " skipping the test";
                  bTestWireCoordinate = false;
                }
                else if (bTestWireCoordinate)
                  throw;
              }
            }
            if (bTestWireCoordinate) {
              if (std::abs(wire_from_wc - expected_wire) > 1e-3) {
                //  throw cet::exception("GeoTestErrorWireCoordinate")
                mf::LogError("GeoTestWireCoordinate")
                  << "wire " << wireID << " [center: " << wire_center << "] on step of " << i << "/"
                  << NSteps << " x" << step << " = " << delta << " cm shows " << wire_from_wc
                  << ", " << expected_wire << " expected.\n";
              } // if mismatch

            } // if testing WireCoordinate

            if ((expected_wire > -0.5) && (expected_wire < NWires - 0.5)) {
              const unsigned int expected_wire_number = std::round(expected_wire);
              unsigned int wire_number_from_wc;
              try {
                wire_number_from_wc = geom->NearestWireID(wire_shifted, planeID).Wire;
              }
              catch (cet::exception& e) {
                throw cet::exception("GeoTestErrorWireCoordinate", "", e)
                  //  MF_LOG_ERROR("GeoTestWireCoordinate")
                  << "wire " << std::string(wireID) << " [center: " << wire_center
                  << "] on step of " << i << "/" << NSteps << " x" << step << " = " << delta
                  << " cm failed NearestWire(), " << expected_wire_number
                  << " expected (more precisely, " << expected_wire << ").\n";
              }

              if (mf::isDebugEnabled()) {
                // In debug mode, we print a lot and we don't (fatally) complain
                std::stringstream e;
                e << "wire " << wireID << " [center: " << wire_center << "] on step of " << i << "/"
                  << NSteps << " x" << step << " = " << delta << " cm near to "
                  << wire_number_from_wc;
                if (wire_number_from_wc != expected_wire_number) {
                  e << ", " << expected_wire_number << " expected (more precisely, "
                    << expected_wire << ").";
                  // throw e;
                  MF_LOG_ERROR("GeoTestWireCoordinate") << e.str();
                }
                else {
                  mf::LogVerbatim("GeoTestWireCoordinate") << e.str();
                }
              }
              else if (wire_number_from_wc != expected_wire_number) {
                // In production mode, we don't print anything and throw on error
                throw cet::exception("GeoTestErrorWireCoordinate")
                  << "wire " << std::string(wireID) << " [center: " << wire_center
                  << "] on step of " << i << "/" << NSteps << " x" << step << " = " << delta
                  << " cm near to " << wire_number_from_wc << ", " << expected_wire_number
                  << " expected (more precisely, " << expected_wire << ").";
              }
            } // if shifted wire not outside boundaries

          } // for i

        } // try
        catch (cet::exception& e) {
          mf::LogWarning("GeoTestCaughtException") << e;
          if (fNonFatalExceptions.count(e.category()) == 0) throw;
        }

      } // for all wires in the plane
    }   // end loop over planes

    stopWatch.Stop();
    MF_LOG_DEBUG("GeoTestWireCoordinate") << "\tdone testing closest channel";
    stopWatch.Print();

    // trigger an exception with NearestChannel
    mf::LogVerbatim("GeoTestWireCoordinate") << "\tattempt to cause an exception to be caught "
                                             << "when looking for a nearest channel";

    // pick a position out of the world
    double posWorld[3];
    geom->WorldBox(nullptr, posWorld + 0, nullptr, posWorld + 1, nullptr, posWorld + 2);
    for (int i = 0; i < 3; ++i)
      posWorld[i] *= 2.;

    bool hasThrown = false;
    unsigned int nearest_to_what = 0;
    try {
      nearest_to_what = geom->NearestChannel(geo::vect::toPoint(posWorld), geom->plane_zero);
    }
    catch (const geo::InvalidWireError& e) {
      mf::LogWarning("GeoTestCaughtException")
        << e << "\nReturned wire would be: " << e.badWire() << ", suggested: " << e.suggestedWire();
      hasThrown = true;
    }
    catch (cet::exception& e) {
      mf::LogWarning("GeoTestCaughtException") << e;
      hasThrown = true;
    }
    if (!hasThrown) {
      if (fDisableValidWireIDcheck) {
        // ok, then why do we disable it?
        // an implementation might prefer to cap the wire number and go on
        // instead of throwing.
        MF_LOG_WARNING("GeoTestWireCoordinate")
          << "GeometryCore::NearestChannel() did not raise an exception"
             " on out-of-world position ("
          << posWorld[0] << "; " << posWorld[1] << "; " << posWorld[2] << "), and returned "
          << nearest_to_what
          << " instead.\n"
             "This is normally considered a failure.";
      }
      else {
        throw cet::exception("GeoTestErrorNearestChannel")
          << "GeometryCore::NearestChannel() did not raise an exception"
             " on out-of-world position ("
          << posWorld[0] << "; " << posWorld[1] << "; " << posWorld[2] << "), and returned "
          << nearest_to_what << " instead\n";
      }
    }
  }

  //......................................................................
  bool GeometryTestAlg::isWireAlignedToPlaneDirections(geo::PlaneGeo const& plane,
                                                       geo::Vector_t const& wireDir) const
  {
    /*
     * Returns `true` if `wireDir` is aligned with plane frame or wire direction
     */

    auto const isOrthogonalTo = [&wireDir](geo::Vector_t const& other) {
      return std::abs(geo::vect::dot(wireDir, other)) < 1.0e-5;
    };

    // we dislike wires aligned to plane frame:
    if (isOrthogonalTo(plane.WidthDir())) return true;
    if (isOrthogonalTo(plane.DepthDir())) return true;

    if (isOrthogonalTo(plane.GetIncreasingWireDirection())) return true;

    return false;

  } // isWireAlignedToPlaneDirections()

  void GeometryTestAlg::testWireIntersection() const
  {
    /*
     * This is a test for geo::GeometryCore::WireIDsIntersect() and
     * geo::WireGeo::IntersectionWith() methods, that return whether two wires
     * intersect, and where.
     *
     * The test strategy is to check all the TPC one by one:
     * - if a query for wires on different cryostats fails
     * - if a query for wires on different TPCs fails
     * - if a query for wires on the same plane fails
     * - for points at the centre of a grid SplitY x SplitZ on the wire planes,
     *   test these point by testWireIntersectionAt() function (see)
     * All tests are performed; at the end, the test is considered a failure
     * if any of the single tests failed.
     */

    unsigned int nErrors = 0;
    for (auto const& TPC : geom->Iterate<geo::TPCGeo>()) {

      MF_LOG_DEBUG("GeometryTest") << "Wire intersection test on " << TPC.ID();

      // sanity: wires on different cryostats
      if (TPC.ID().Cryostat < geom->Ncryostats() - 1) {

        geo::WireID const w1{geo::PlaneID{TPC.ID(), 0}, 0};
        geo::WireGeo const& wire1 = geom->Wire(w1);
        geo::Vector_t const& wireDir = wire1.Direction();

        geo::CryostatID const otherCryoID{TPC.ID().Cryostat + 1};
        geo::CryostatGeo const& otherCryo = geom->Cryostat(otherCryoID);
        geo::PlaneGeo const* otherPlane = nullptr;
        for (auto const& plane : geom->Iterate<geo::PlaneGeo>(otherCryo.ID())) {
          if (isWireAlignedToPlaneDirections(plane, wireDir)) continue;
          otherPlane = &plane;
          break;
        }
        if (otherPlane) {
          geo::WireID const w2{otherPlane->ID(), 1};
          MF_LOG_TRACE("GeometryTest") << "Off cryostat test (" << w1 << "): chosen wire " << w2;
          geo::Point_t xingPoint;
          if (geom->WireIDsIntersect(w1, w2, xingPoint)) {
            MF_LOG_ERROR("GeometryTest")
              << "WireIDsIntersect() on " << w1 << " and " << w2 << " returned " << xingPoint
              << " cm, while should have reported no intersection at all";
            ++nErrors;
          } // if intersect
          try {
            // the value of result is not checked here
            wire1.IntersectionWith(geom->Wire(w2));
          }
          catch (...) {
            MF_LOG_ERROR("GeometryTest") << "WiresIntersect() on " << w1 << " and " << w2
                                         << " threw an exception, which should not have";
            ++nErrors;
          }
        }
        else {
          MF_LOG_WARNING("GeometryTest")
            << "No wire plane found in " << otherCryo.ID() << " with wires not aligned with "
            << w1.asPlaneID() << ", " << wireDir << "; off-cryostat sanity check skipped.";
        }

      } // if not the last cryostat

      // sanity: wires on different TPC
      if (TPC.ID().TPC < geom->NTPC(TPC.ID().asCryostatID()) - 1) {

        geo::PlaneID const refPlaneID{TPC.ID(), 0};
        geo::WireID const w1{refPlaneID, 0};
        geo::WireGeo const& wire1 = geom->Wire(w1);
        geo::Vector_t const& wireDir = wire1.Direction();

        geo::CryostatGeo const& cryo = geom->Cryostat(TPC.ID());
        geo::PlaneGeo const* otherPlane = nullptr;
        for (auto const& plane : geom->Iterate<geo::PlaneGeo>(cryo.ID())) {
          if (plane.ID().asTPCID() == TPC.ID()) continue; // on the same TPC
          if (isWireAlignedToPlaneDirections(plane, wireDir)) continue;
          otherPlane = &plane;
          break;
        }
        if (otherPlane) {
          geo::WireID const w2{otherPlane->ID(), 1};
          MF_LOG_TRACE("GeometryTest") << "Off TPC test (" << w1 << "): chosen wire " << w2;
          geo::Point_t xingPoint;
          if (geom->WireIDsIntersect(w1, w2, xingPoint)) {
            MF_LOG_ERROR("GeometryTest")
              << "WireIDsIntersect() on " << w1 << " and " << w2 << " returned " << xingPoint
              << ", while should have reported no intersection at all";
            ++nErrors;
          } // if intersect
          try {
            // the value of result is not checked here
            wire1.IntersectionWith(geom->Wire(w2));
          }
          catch (...) {
            MF_LOG_ERROR("GeometryTest") << "WiresIntersect() on " << w1 << " and " << w2
                                         << " threw an exception, which should not have";
            ++nErrors;
          }
        }
        else {
          MF_LOG_WARNING("GeometryTest")
            << "No wire plane found in any TPC of " << cryo.ID() << " with wires not aligned with "
            << refPlaneID << ", " << wireDir << "; off-TPC sanity check skipped.";
        }
      } // if not the last TPC

      // sanity: wires on same plane
      for (geo::PlaneGeo const& plane : TPC.IteratePlanes()) {
        geo::WireID const w1{plane.ID(), 0}, w2{plane.ID(), 1};
        geo::Point_t xingPoint;
        if (geom->WireIDsIntersect(w1, w2, xingPoint)) {
          MF_LOG_ERROR("GeometryTest")
            << "WireIDsIntersect() on " << w1 << " and " << w2 << " returned " << xingPoint
            << ", while should have reported no intersection at all";
          ++nErrors;
        } // if intersect
        // prerequisites of WiresIntersect() are not met here, no test possible
      } // for all planes

      // sample the area covered by all planes, split into SplitW and SplitD
      // rectangles; drift coordinate is chosen roughly in the middle of the TPC
      geo::PlaneGeo const& refPlane = TPC.SmallestPlane();
      constexpr unsigned int SplitW = 19, SplitD = 17;

      auto const driftOffset = -TPC.DriftDistance() / 2.0 * TPC.DriftDir();
      auto const refPoint = refPlane.GetCenter() + driftOffset;

      decltype(auto) coverage = refPlane.ActiveArea();
      const double stepW = coverage.width.length() / SplitW;
      const double stepD = coverage.depth.length() / SplitD;
      const int stepsW = SplitW / 2;
      const int stepsD = SplitD / 2;

      // let's pick a point:
      for (int iW = -stepsW; iW <= +stepsW; ++iW) {

        auto const widthOffset = (iW * stepW) * refPlane.WidthDir();

        for (int iD = -stepsD; iD < +stepsD; ++iD) {

          auto const depthOffset = (iD * stepD) * refPlane.DepthDir();

          auto const point = refPoint + widthOffset + depthOffset;

          // finally, let's test this point...
          nErrors += testWireIntersectionAt(TPC, point);
        } // for y
      }   // for z
    }     // for TPC

    if (nErrors > 0) {
      throw cet::exception("GeoTestWireIntersection")
        << "Accumulated " << nErrors << " errors (see messages above)\n";
    }

  } // GeometryTestAlg::testWireIntersection()

  unsigned int GeometryTestAlg::testWireIntersectionAt(const geo::TPCGeo& TPC,
                                                       geo::Point_t const& point) const
  {
    /* Tests WireIDsIntersect() and WiresIntersect() on the specified point on
     * the wire planes of a given TPC.
     *
     * The test follows this strategy:
     * - find the ID of the wires closest to the point on each plane
     * - for all wire plane pairing, ask for the intersection between the wires
     * - fail if the returned point is farther than half a pitch from the
     *   original point
     *
     * This function returns the number of accumulated failures.
     */

    using geo::vect::dot;

    unsigned int nErrors = 0;

    const unsigned int NPlanes = TPC.Nplanes();

    bool const bDriftOnX = (TPC.DriftDir() == geo::Xaxis()) || (TPC.DriftDir() == -geo::Xaxis());

    // collect information per plane:
    std::vector<double> WirePitch(NPlanes); // for convenience
    std::vector<geo::Vector_t> WireCoordDirs(NPlanes);
    std::vector<geo::WireID> WireIDs; // ID of the closest wire
    WireIDs.reserve(NPlanes);
    std::vector<double> WireDistances(NPlanes); // distance from the closest wire
    for (unsigned int iPlane = 0; iPlane < NPlanes; ++iPlane) {
      const geo::PlaneGeo& plane = TPC.Plane(iPlane);
      WireCoordDirs[iPlane] = plane.GetIncreasingWireDirection();
      WirePitch[iPlane] = plane.WirePitch();

      const double WireDistance = geom->WireCoordinate(point, plane.ID());
      WireIDs.emplace_back(plane.ID(), (unsigned int)std::round(WireDistance));
      WireDistances[iPlane] = (WireDistance - std::round(WireDistance)) * WirePitch[iPlane];

      MF_LOG_DEBUG("GeometryTest")
        << "Nearest wire to " << point << " on plane " << std::string(plane.ID())
        << " (pitch: " << WirePitch[iPlane] << ", coord.dir.=" << WireCoordDirs[iPlane] << ") is "
        << WireIDs[iPlane] << " (position: " << WireDistance << ")";
    } // for planes

    // test all the combinations

    lar::util::RealComparisons<double> coordIs(1e-3);
    auto vectorIs = lar::util::makeVector3DComparison(coordIs);
    for (unsigned int iPlane1 = 0; iPlane1 < NPlanes; ++iPlane1) {

      const geo::WireID& w1 = WireIDs[iPlane1];
      geo::PlaneGeo const& plane1 = TPC.Plane(w1);
      geo::WireGeo const& w1obj = plane1.Wire(w1);

      for (unsigned int iPlane2 = iPlane1 + 1; iPlane2 < NPlanes; ++iPlane2) {
        const geo::WireID& w2 = WireIDs[iPlane2];
        geo::PlaneGeo const& plane2 = TPC.Plane(w2);
        geo::WireGeo const& w2obj = plane2.Wire(w2);

        geo::Point_t xingPoint;
        if (!geom->WireIDsIntersect(w1, w2, xingPoint)) {
          MF_LOG_ERROR("GeometryTest")
            << "Wires " << w1 << " and " << w2 << " should intersect around " << point << " of TPC "
            << TPC.ID() << ", but they seem not to intersect at all!";
          ++nErrors;
          continue;
        }
        geo::Point_t xingPoint2 = xingPoint; // matching point on plane 2
        plane2.DriftPoint(xingPoint2);

        if (bDriftOnX) { // legacy code check

          geo::WireIDIntersection widIntersect;
          if (!geom->WireIDsIntersect(w1, w2, widIntersect)) {
            MF_LOG_ERROR("GeometryTest")
              << "Legacy check: wires " << w1 << " and " << w2 << " should intersect around "
              << point << " of TPC " << TPC.ID() << ", but they seem not to intersect at all!";
            ++nErrors;
          }

          if (coordIs.nonEqual(widIntersect.y, xingPoint.Y()) ||
              coordIs.nonEqual(widIntersect.z, xingPoint.Z())) {
            MF_LOG_ERROR("GeometryTest")
              << "Legacy check: wires " << w1 << " and " << w2 << " should intersect around "
              << point << " of TPC " << TPC.ID() << ", but legacy code says (?, " << widIntersect.y
              << ", " << widIntersect.z << ")!";
            ++nErrors;
          }

        } // bDriftOnX

        geo::Point_t xingPointInv;
        if (!geom->WireIDsIntersect(w2, w1, xingPointInv)) {
          MF_LOG_ERROR("GeometryTest")
            << "Wires " << w2 << " and " << w1 << " (reversed test) should intersect around "
            << point << " of TPC " << TPC.ID() << ", but they seem not to intersect at all!";
          ++nErrors;
          continue;
        }
        if (vectorIs.nonEqual(xingPointInv, xingPoint2)) {
          MF_LOG_ERROR("GeometryTest")
            << "WireIDsIntersect() gives different intersections for " << w1 << " and " << w2
            << ": " << xingPoint2 << " (direct+shift) and " << xingPointInv << " (reversed)";
          ++nErrors;
          continue;
        }

        // the expected distance between the probe point and the
        // intersection point is geometrically determined, given the distances
        // of the probe point from the two wires and the angle between the wires
        // the formula is a mix between the Carnot theorem and sine definition;
        // the definition of the angles is tricky though, and we rely on the
        // strong vector definition enforced in the geometry to get the proper
        // "cosine" and corresponding sine
        const double d1 = WireDistances[iPlane1], d2 = WireDistances[iPlane2],
                     cosAlpha = dot(WireCoordDirs[iPlane1], WireCoordDirs[iPlane2]);
        const double expected_d =
          std::sqrt((cet::square(d1) + cet::square(d2) - 2.0 * d1 * d2 * cosAlpha) /
                    (1 - cet::square(cosAlpha)));
        // the actual distance we have found:
        double const d = plane1.Projection(xingPoint - point).R();
        MF_LOG_DEBUG("GeometryTest")
          << " - wires " << w1 << " and " << w2 << " intersect at " << xingPoint << ", " << d
          << " cm far from starting point (expected: " << expected_d << ")";

        // precision of the test is an issue; the 10^-3 x pitch threshold
        // is roughly tuned so that we don't get errors
        lar::util::RealComparisons<double> wireCoordIs(
          std::max(WirePitch[iPlane1], WirePitch[iPlane2]) * 1e-3); // cm
        if (wireCoordIs.nonEqual(d, expected_d)) {
          MF_LOG_ERROR("GeometryTest")
            << "wires " << w1 << " and " << w2 << " intersect at " << xingPoint << ", " << d
            << " cm far from starting point: too far from the expected " << expected_d << " cm!";
          ++nErrors;
          continue;
        } // if too far

        geo::Point_t objXingPoint;

        // test that geo::WiresIntersection() gives the same result as
        // the already validated one from geom->WireIDsIntersect()
        objXingPoint = geo::WiresIntersection(w1obj, w2obj);
        if (vectorIs.nonEqual(objXingPoint, xingPoint)) {
          MF_LOG_ERROR("GeometryTest")
            << "geo::WiresIntersection() gives wrong intersection for " << w1 << " and " << w2
            << ": " << objXingPoint << " vs. " << xingPoint << " (expected)";
          ++nErrors;
          continue;
        }

        // test that geo::WireGeo::IntersectionWith() gives the same result as
        // the already validated one from geom->WireIDsIntersect()
        objXingPoint = w1obj.IntersectionWith(w2obj);
        if (vectorIs.nonEqual(objXingPoint, xingPoint)) {
          MF_LOG_ERROR("GeometryTest")
            << "geo::WireGeo[" << w1 << "]::IntersectionWith() gives wrong intersection with " << w2
            << ": " << objXingPoint << " vs. " << xingPoint << " (expected)";
          ++nErrors;
          continue;
        }

        objXingPoint = w2obj.IntersectionWith(w1obj);
        if (vectorIs.nonEqual(objXingPoint, xingPointInv)) {
          MF_LOG_ERROR("GeometryTest")
            << "geo::WireGeo[" << w2 << "]::IntersectionWith() gives wrong intersection with " << w1
            << ": " << objXingPoint << " vs. " << xingPoint << " (expected)";
          ++nErrors;
          continue;
        }

      } // for iPlane2
    }   // for iPlane1

    return nErrors;
  } // GeometryTestAlg::testWireIntersectionAt()

  //......................................................................
  void GeometryTestAlg::testThirdPlane() const
  {
    /*
     * This is a test for ThirdPlane() function, that returns the plane that is
     * not specified in the input.
     * Currently, the only implemented signature is designed for TPCs with 3
     * planes.
     *
     * The test strategy is to check all the TPC one by one:
     * - for all combinations of two planes, if the result is the expected one
     *
     * All tests are performed; at the end, the test is considered a failure
     * if any of the single tests failed.
     */

    unsigned int nErrors = 0;
    for (auto const& tpcid : geom->Iterate<geo::TPCID>()) {
      const geo::TPCGeo& TPC = geom->TPC(tpcid);

      const unsigned int nPlanes = TPC.Nplanes();
      MF_LOG_DEBUG("GeometryTest") << tpcid << " (" << nPlanes << " planes)";

      for (geo::PlaneID::PlaneID_t iPlane1 = 0; iPlane1 < nPlanes; ++iPlane1) {
        geo::PlaneID pid1(tpcid, iPlane1);

        for (geo::PlaneID::PlaneID_t iPlane2 = 0; iPlane2 < nPlanes; ++iPlane2) {
          geo::PlaneID pid2(tpcid, iPlane2);

          const bool isValidInput = (nPlanes == 3) && (iPlane1 != iPlane2);
          bool bError = false;
          geo::PlaneID third_plane;
          try {
            third_plane = geom->ThirdPlane(pid1, pid2);
          }
          catch (cet::exception const& e) {
            if (isValidInput) throw;
            // we have gotten the error we were looking for
            // if "GeometryCore" is included in the categories of the exception
            bError = hasCategory(e, "GeometryCore");
          } // try...catch

          MF_LOG_TRACE("GeometryTest") << "  [" << pid1 << "], [" << pid2 << "] => "
                                       << (bError ? "error" : std::string(third_plane));
          if (bError) continue; // we got what we wanted

          if (!bError && !isValidInput) {
            MF_LOG_ERROR("GeometryTest")
              << "ThirdPlane() on " << pid1 << " and " << pid2 << " returned " << third_plane
              << ", while should have thrown an exception";
            ++nErrors;
            continue;
          } // if no error

          if (third_plane != tpcid) {
            MF_LOG_ERROR("GeometryTest")
              << "ThirdPlane() on " << pid1 << " and " << pid2 << " returned " << third_plane
              << ", on a different TPC!!!";
            ++nErrors;
          }
          else if (!third_plane.isValid) {
            MF_LOG_ERROR("GeometryTest") << "ThirdPlane() on " << pid1 << " and " << pid2
                                         << " returned an invalid " << third_plane;
            ++nErrors;
          }
          else if (third_plane.Plane >= nPlanes) {
            MF_LOG_ERROR("GeometryTest")
              << "ThirdPlane() on " << pid1 << " and " << pid2 << " returned " << third_plane
              << " with plane out of range";
            ++nErrors;
          }
          else if (third_plane == pid1) {
            MF_LOG_ERROR("GeometryTest")
              << "ThirdPlane() on " << pid1 << " and " << pid2 << " returned " << third_plane
              << ", same as the first input";
            ++nErrors;
          }
          else if (third_plane == pid2) {
            MF_LOG_ERROR("GeometryTest")
              << "ThirdPlane() on " << pid1 << " and " << pid2 << " returned " << third_plane
              << ", same as the second input";
            ++nErrors;
          }

        } // for plane 2

      } // for plane 1

    } // for TPC

    if (nErrors > 0) {
      throw cet::exception("GeoTestThirdPlane")
        << "Accumulated " << nErrors << " errors (see messages above)\n";
    }

  } // GeometryTestAlg::testThirdPlane()

  //......................................................................
  void GeometryTestAlg::testThirdPlane_dTdW() const
  {
    /*
     * This is a test for ThirdPlane_dTdW() function, that returns the apparent
     * slope on a wire plane, given the ones observed on other two planes.
     *
     * The test strategy is to check all the TPC one by one:
     * - if a query for planes on different cryostats fails
     * - if a query for planes on different TPCs fails
     * - for selected 3D points, compute the three dT/dW and verify them by
     *   test these slopes by testThirdPlane__dTdW_at() function (see)
     *
     * All tests are performed; at the end, the test is considered a failure
     * if any of the single tests failed.
     */

    unsigned int nErrors = 0;
    for (geo::TPCID const& tpcid : geom->Iterate<geo::TPCID>()) {
      const geo::TPCGeo& TPC = geom->TPC(tpcid);

      const double driftVelocity = 0.1 * ((TPC.DriftDirection() == geo::kNegX) ? -1. : +1);

      const unsigned int NPlanes = TPC.Nplanes();
      MF_LOG_DEBUG("GeometryTest") << tpcid << " (" << NPlanes << " planes)";

      // sanity: planes on different cryostats
      if (tpcid.Cryostat < geom->Ncryostats() - 1) {
        geo::PlaneID p1{tpcid, 0}, p2{tpcid.Cryostat + 1, tpcid.TPC, 1};
        bool bError = false;
        double slope;
        try {
          slope = geom->ThirdPlane_dTdW(p1, +1.0, p2, -1.0);
        }
        catch (cet::exception const& e) {
          // we have gotten the error we were looking for
          // if "GeometryCore" is included in the categories of the exception
          bError = hasCategory(e, "GeometryCore");
        } // try...catch
        if (!bError) {
          MF_LOG_ERROR("GeometryTest")
            << "ThirdPlane_dTdW() on " << p1 << " and " << p2 << " returned " << slope
            << ", while should have thrown an exception";
          ++nErrors;
        } // if no error
      }   // if not the last cryostat

      // sanity: planes on different TPC
      if (tpcid.TPC < geom->NTPC(tpcid.asCryostatID()) - 1) {
        geo::PlaneID p1{tpcid, 0}, p2{tpcid.Cryostat, tpcid.TPC + 1, 1};
        bool bError = false;
        double slope;
        try {
          slope = geom->ThirdPlane_dTdW(p1, +1.0, p2, -1.0);
        }
        catch (cet::exception const& e) {
          // we have gotten the error we were looking for
          // if "GeometryCore" is included in the categories of the exception
          bError = hasCategory(e, "GeometryCore");
        } // try...catch
        if (!bError) {
          MF_LOG_ERROR("GeometryTest")
            << "ThirdPlane_dTdW() on " << p1 << " and " << p2 << " returned " << slope
            << ", while should have thrown an exception";
          ++nErrors;
        } // if no error
      }   // if not the last TPC in its cryostat

      // pick a point in the very middle of the TPC:
      const geo::Point_t A{TPC.CenterX(), TPC.CenterY(), TPC.CenterZ()};
      // pick a radius half the way to the closest border
      const double radius = std::min({TPC.HalfWidth(), TPC.HalfHeight(), TPC.Length() / 2.}) / 2.;

      // I arbitrary decide that the second point will have dX equal size as
      // the radius, and on the positive x direction (may be negative dT)
      const double dX = radius;
      const double dT = driftVelocity * dX;

      // sample a circle of SplitAngles directions around A
      constexpr unsigned int NAngles = 19;
      const double start_angle = 0.05;                             // radians
      const double step_angle = 2. * util::pi<double>() / NAngles; // radians

      for (unsigned int iAngle = 0; iAngle < NAngles; ++iAngle) {
        const double angle = start_angle + iAngle * step_angle;

        // define B as a point "radius" far from A in the angle direction,
        // with some arbitrary and fixed dx offset
        geo::Point_t const B{
          A.X() + dX, A.Y() + radius * std::sin(angle), A.Z() + radius * std::cos(angle)};

        // get the expectation; this function assumes a drift velocity of
        // 1 mm per tick by default; for the test, it does not matter
        std::vector<std::pair<geo::PlaneID, double>> dTdWs =
          ExpectedPlane_dTdW(A, B, driftVelocity);

        if (mf::isDebugEnabled()) {
          mf::LogTrace log("GeometryTest");
          log << "Expected dT/dW for a segment with " << radius << " cm long projection at "
              << angle << " rad, and dT " << dT << " cm:";
          for (auto const& dTdW_info : dTdWs)
            log << "  " << dTdW_info.first << " slope:" << dTdW_info.second;
        } // if debug

        // run the test
        nErrors += testThirdPlane_dTdW_at(dTdWs);

      } // for angle

    } // for TPC

    if (nErrors > 0) {
      throw cet::exception("GeoTestThirdPlane_dTdW")
        << "Accumulated " << nErrors << " errors (see messages above)\n";
    }

  } // GeometryTestAlg::testThirdPlane_dTdW()

  std::vector<std::pair<geo::PlaneID, double>> GeometryTestAlg::ExpectedPlane_dTdW(
    geo::Point_t const& A,
    geo::Point_t const& B,
    const double driftVelocity /* = 0.1 */
  ) const
  {
    // This function returns a list of entries, one for each plane:
    // - plane ID
    // - slope of the projection of AB from the plane, in dt/dw units

    // find which TPC we are taking about
    geo::TPCID tpcid = geom->FindTPCAtPosition(A);

    if (!tpcid.isValid) {
      throw cet::exception("GeometryTestAlg")
        << "ExpectedPlane_dTdW(): can't find any TPC containing point A (" << A.X() << "; " << A.Y()
        << "; " << A.Z() << ")";
    }

    if (geom->FindTPCAtPosition(B) != tpcid) {
      throw cet::exception("GeometryTestAlg")
        << "ExpectedPlane_dTdW(): point A (" << A.X() << "; " << A.Y() << "; " << A.Z()
        << ") is in " << std::string(tpcid) << " while point B (" << B.X() << "; " << B.Y() << "; "
        << B.Z() << ") is in " << std::string(geom->FindTPCAtPosition(B));
    }

    geo::TPCGeo const& TPC = geom->TPC(tpcid);

    // conversion from X coordinate to tick coordinate
    double dT_over_dX = 1. / driftVelocity;
    switch (TPC.DriftDirection()) {
    case geo::kPosX:
      // if the drift direction is toward positive x, higher x will reach the
      // plane earlier and have smaller t, hence the flip in sign
      dT_over_dX = -dT_over_dX;
      break;
    case geo::kNegX: break;
    default:
      throw cet::exception("InternalError")
        << "GeometryTestAlg::ExpectedPlane_dTdW(): drift direction #" << ((int)TPC.DriftDirection())
        << " of " << std::string(tpcid) << " not supported.\n";
    } // switch drift direction

    const unsigned int nPlanes = TPC.Nplanes();
    std::vector<std::pair<geo::PlaneID, double>> slopes(nPlanes);

    for (geo::PlaneID::PlaneID_t iPlane = 0; iPlane < nPlanes; ++iPlane) {
      geo::PlaneID pid(tpcid, iPlane);
      const double wA = geom->WireCoordinate(A, pid);
      const double wB = geom->WireCoordinate(B, pid);

      slopes[iPlane] = std::make_pair(pid, ((B.X() - A.X()) * dT_over_dX) / (wB - wA));

    } // for iPlane

    return slopes;
  } // GeometryTestAlg::ExpectedPlane_dTdW()

  unsigned int GeometryTestAlg::testThirdPlane_dTdW_at(
    std::vector<std::pair<geo::PlaneID, double>> const& plane_dTdW) const
  {
    /*
     * This function tests that for every combination of planes, the expected
     * slope is returned within some tolerance.
     * It returns the number of mistakes found.
     *
     * The parameter is a list if pair of expected slope on the paired plane.
     */

    unsigned int nErrors = 0;
    for (std::pair<geo::PlaneID, double> const& input1 : plane_dTdW) {
      for (std::pair<geo::PlaneID, double> const& input2 : plane_dTdW) {

        const bool bValidInput = input1.first != input2.first;

        for (std::pair<geo::PlaneID, double> const& output : plane_dTdW) {
          bool bError = false;
          double output_slope = 0.;
          try {
            output_slope = geom->ThirdPlane_dTdW(
              input1.first, input1.second, input2.first, input2.second, output.first);
          }
          catch (cet::exception const& e) {
            // catch only "GeometryCore" category, and only if input is faulty;
            // otherwise, rethrow
            if (bValidInput) throw;
            bError = hasCategory(e, "GeometryCore");
            if (!bError) throw;
            MF_LOG_TRACE("GeometryTest")
              << input1.first << " slope:" << input1.second << "  " << input2.first
              << " slope:" << input2.second << "  => exception";
            continue;
          }

          if (!bValidInput && !bError) {
            MF_LOG_ERROR("GeometryTest")
              << "GeometryCore::ThirdPlane_dTdW() on " << input1.first << " and " << input2.first
              << " should have thrown an exception, it returned " << output_slope << " instead";
            ++nErrors;
            continue;
          } // if faulty input and no error

          MF_LOG_TRACE("GeometryTest")
            << input1.first << " slope:" << input1.second << "  " << input2.first
            << " slope:" << input2.second << "  => " << output.first << " slope:" << output_slope;
          if (((output.second == 0.) && (output_slope > 1e-3)) ||
              std::abs(output_slope / output.second - 1.) > 1e-3) {
            MF_LOG_ERROR("testThirdPlane_dTdW_at")
              << "GeometryCore::ThirdPlane_dTdW(): " << input1.first << " slope:" << input1.second
              << "  " << input2.first << " slope:" << input2.second << "  => " << output.first
              << " slope:" << output_slope << "  (expected: " << output.second << ")";
            ++nErrors;
          } // if too far

          // now test the automatic detection of the other plane

        } // for output
      }   // for second input
    }     // for first input
    return nErrors;
  } // GeometryTestAlg::testThirdPlane_dTdW_at()

  //......................................................................
  void GeometryTestAlg::testWirePitch()
  {
    // loop over all planes and wires to be sure the pitch is consistent
    unsigned int nPitchErrors = 0;

    if (fExpectedWirePitches.empty()) {
      // hard code the value we think it should be for each detector;
      // this is legacy and you should not add anything:
      // add the expectation to the FHiCL configuration of the test instead
      if (geom->DetectorName() == "bo") { fExpectedWirePitches = {0.46977, 0.46977, 0.46977}; }
      if (!fExpectedWirePitches.empty()) {
        mf::LogInfo("WirePitch")
          << "Using legacy wire pitch parameters hard-coded for the detector '"
          << geom->DetectorName() << "'";
      }
    }
    if (fExpectedWirePitches.empty()) {
      mf::LogWarning("WirePitch")
        << "no expected wire pitch; I'll just check that they are all the same";
    }
    else {
      mf::LogInfo log("WirePitch");
      log << "Expected wire pitch per plane, in centimetres:";
      for (double pitch : fExpectedWirePitches)
        log << " " << pitch;
      log << " [...]";
    }

    for (geo::PlaneID const& planeid : geom->Iterate<geo::PlaneID>()) {

      geo::PlaneGeo const& plane = geom->Plane(planeid);
      const unsigned int nWires = plane.Nwires();
      if (nWires < 2) continue;

      geo::WireGeo const* pWire = &(plane.Wire(0));

      // which pitch to expect:
      // - if they did not tell us anything:
      //     get the one from the first two wires
      // - if they did tell something, but not for this plane:
      //     get the last pitch they told us
      // - if they told us about this plane: well, then use it!
      double expectedPitch = 0.;
      if (fExpectedWirePitches.empty()) {
        geo::WireGeo const& wire1 = plane.Wire(1); // pWire now points to wire0
        expectedPitch = geo::WireGeo::WirePitch(*pWire, wire1);
        MF_LOG_DEBUG("WirePitch") << "Wire pitch on " << planeid << ": " << expectedPitch << " cm";
      }
      else if (planeid.Plane < fExpectedWirePitches.size())
        expectedPitch = fExpectedWirePitches[planeid.Plane];
      else
        expectedPitch = fExpectedWirePitches.back();

      geo::WireID::WireID_t w = 0; // wire number
      while (++w < nWires) {
        geo::WireGeo const* pPrevWire = pWire;
        pWire = &(plane.Wire(w));

        const double thisPitch = std::abs(pWire->DistanceFrom(*pPrevWire));
        if (std::abs(thisPitch - expectedPitch) > 1e-5) {
          mf::LogProblem("WirePitch")
            << "ERROR: on plane " << planeid << " pitch between wires W:" << (w - 1)
            << " and W:" << w << " is " << thisPitch << " cm, not " << expectedPitch
            << " as expected!";
          ++nPitchErrors;
        } // if unexpected pitch
      }   // while

    } // for

    if (nPitchErrors > 0) {
      throw cet::exception("UnexpectedWirePitch")
        << "unexpected pitches between " << nPitchErrors << " wires!";
    } // end loop over planes

  } // GeometryTestAlg::testWirePitch()

  //......................................................................
  void GeometryTestAlg::testInterWireProjectedDistance() const
  {

    /*
     * For each wire plane:
     *  * we pick some projected directions; for each one:
     *     * we manipulate it in some way that does not affect the result
     *     * check that the projected distance is as expected
     *     * add some arbitrary component on the drift direction; for each one:
     *       * check that the projected distance is as expected
     *       * check that the 3D distance is as expected
     *
     * We do not test directions parallel to the wires because they get
     * numerically unstable and the expectation may potentially differ a lot
     * being calculated with a different procedure.
     */

    constexpr lar::util::RealComparisons cmp{1e-4};

    static double const V3 = std::sqrt(3.0);

    // BUG the deduction guide for std::array seems not to be implemented yet
    //     in Clang 5.0.0
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    //     std::array const testProjections = {
    std::array<geo::PlaneGeo::WireCoordProjection_t, 5U> const testProjections{{
      geo::PlaneGeo::WireCoordProjection_t{0.00, 1.00},
      geo::PlaneGeo::WireCoordProjection_t{0.75, 1.00},
      geo::PlaneGeo::WireCoordProjection_t{1.00, 1.00},
      geo::PlaneGeo::WireCoordProjection_t{V3, 1.00},
      geo::PlaneGeo::WireCoordProjection_t{2.40, 1.00}
      //     };
    }};
    // BUG the deduction guide for std::array seems not to be implemented yet
    //     in Clang 5.0.0
    // BUG the double brace syntax is required to work around clang bug 21629
    // (https://bugs.llvm.org/show_bug.cgi?id=21629)
    //     constexpr std::array testDriftOffsets { -20.0, -10.0, 0.0, +10.0, +20.0 };
    constexpr std::array<double, 5U> testDriftOffsets{{-20.0, -10.0, 0.0, +10.0, +20.0}};

    unsigned int nErrors = 0; // error count for the final report

    for (auto const& plane : geom->Iterate<geo::PlaneGeo>()) {

      double const pitch = plane.WirePitch();
      auto const normalDir = plane.GetNormalDirection();

      for (auto const& testProjBase : testProjections) {

        //
        // expected result is kind of encoded in the chosen projections
        //
        double const expected = testProjBase.R() * pitch;
        double const expectedSqr = cet::square(expected);

        // we flip the projection around: result should not change
        for (double dirL : {-1.0, 1.0})
          for (double dirW : {-1.0, 1.0})
            for (double scale : {0.5, 1.0, 3.0}) {

              //
              // test the projection directly
              //
              geo::PlaneGeo::WireCoordProjection_t const testProj{scale * dirL * testProjBase.X(),
                                                                  scale * dirW * testProjBase.Y()};

              double const interWireFromProj = plane.InterWireProjectedDistance(testProj);
              if (cmp.nonEqual(interWireFromProj, expected)) {
                mf::LogProblem("")
                  << "ERROR: on plane " << plane.ID()
                  << " distance between wires on projected direction " << testProj << " is "
                  << interWireFromProj << " cm (expected: " << expected << ")";
                ++nErrors;
              } // if unexpected result

              // this is how much we needed to expand the test direction vector
              // (happens to work for the special case expected = 0 too)
              double const dScale = expected / testProj.R();

              //
              // test a 3D direction
              //

              auto const baseDir = plane.ComposeVector(0.0, testProj);

              for (double const driftOffset : testDriftOffsets) {
                auto const testDir = baseDir + driftOffset * normalDir;

                double const interWire = plane.InterWireProjectedDistance(testDir);
                if (cmp.nonEqual(interWire, expected)) {
                  mf::LogProblem("")
                    << "ERROR: on plane " << plane.ID()
                    << " projected distance between wires on direction " << testDir
                    << " (from projection " << testProj << " and offset " << driftOffset << ") is "
                    << interWire << " cm (expected: " << expected << ")";
                  ++nErrors;
                } // if unexpected result

                // build the expectation adding the drift component to the projected
                double const expected3D =
                  std::sqrt(expectedSqr + cet::square(driftOffset * dScale));

                double const interWire3D = plane.InterWireDistance(testDir);
                if (cmp.nonEqual(interWire3D, expected3D)) {
                  mf::LogProblem("")
                    << "ERROR: on plane " << plane.ID() << " distance between wires on direction "
                    << testDir << " (from projection " << testProj << " and offset " << driftOffset
                    << ") is " << interWire3D << " cm (expected: " << expected3D << ")";
                  ++nErrors;
                } // if unexpected result

              } // for drifts

            } // for dirL, dirW and scale

      } // for all test projection directions

    } // for all planes

    if (nErrors > 0U) {
      throw cet::exception("InterWireProjectedDistance")
        << "unexpected distances in " << nErrors << " tests!";
    } // end loop over planes

  } // GeometryTestAlg::testInterWireProjectedDistance()

  //......................................................................
  void GeometryTestAlg::testPlanePitch()
  {
    // loop over all planes to be sure the pitch is consistent

    if (fExpectedPlanePitches.empty()) {
      // hard code the value we think it should be for each detector;
      // this is legacy and you should not add anything:
      // add the expectation to the FHiCL configuration of the test instead
      if (geom->DetectorName() == "bo") { fExpectedPlanePitches = {0.65}; }
      if (!fExpectedPlanePitches.empty()) {
        mf::LogInfo("PlanePitch")
          << "Using legacy plane pitch parameters hard-coded for the detector '"
          << geom->DetectorName() << "'";
      }
    }
    if (fExpectedPlanePitches.empty()) {
      mf::LogWarning("PlanePitch") << "no expected plane pitch;"
                                      " I'll just check that they are all the same";
    }
    else {
      mf::LogInfo log("PlanePitch");
      log << "Expected plane pitch per plane pair, in centimetres:";
      for (double pitch : fExpectedPlanePitches)
        log << " " << pitch;
      log << " [...]";
    }

    unsigned int nPitchErrors = 0;
    for (geo::TPCID const& tpcid : geom->Iterate<geo::TPCID>()) {

      geo::TPCGeo const& TPC = geom->TPC(tpcid);
      const unsigned int nPlanes = TPC.Nplanes();
      if (nPlanes < 2) continue;

      double expectedPitch = 0.;
      if (fExpectedPlanePitches.empty()) {
        expectedPitch = TPC.PlanePitch(0, 1);
        MF_LOG_DEBUG("PlanePitch") << "Plane pitch between the first two planes of " << tpcid
                                   << ": " << expectedPitch << " cm";
      }

      geo::PlaneID::PlaneID_t p = 0; // plane number
      while (++p < nPlanes) {
        // which pitch to expect:
        // - if they did not tell us anything:
        //     use the one from the first two planes (already in expectedPitch)
        // - if they did tell something, but not for this plane:
        //     get the last pitch they told us
        // - if they told us about this plane: well, then use it!
        if (!fExpectedPlanePitches.empty()) {
          if (p - 1 < fExpectedPlanePitches.size())
            expectedPitch = fExpectedPlanePitches[p - 1];
          else
            expectedPitch = fExpectedPlanePitches.back();
        } // if we have directions about plane pitch

        const double thisPitch = std::abs(TPC.PlanePitch(p - 1, p));
        if (std::abs(thisPitch - expectedPitch) > 1e-5) {
          mf::LogProblem("PlanePitch")
            << "ERROR: pitch of planes P:" << (p - 1) << " and P: " << p << " in " << tpcid
            << " is " << thisPitch << " cm, not " << expectedPitch << " as expected!";
          ++nPitchErrors;
        } // if unexpected pitch
      }   // while planes

    } // for TPCs

    if (nPitchErrors > 0) {
      throw cet::exception("UnexpectedPlanePitch")
        << "unexpected pitches between " << nPitchErrors << " planes!";
    } // end loop over planes

  } // GeometryTestAlg::testPlanePitch()

  //......................................................................

  void GeometryTestAlg::testStepping()
  {
    //
    // Test stepping. Example is similar to what one would do for photon
    // transport. Rattles photons around inside the scintillator
    // bouncing them off walls.
    //
    geo::WireGeo::LocalPoint_t const xyzWire{};
    geo::WireGeo::LocalVector_t const dxyzWire{0, sin(0.1), cos(0.1)};

    geo::PlaneID const planeid{0, 0, 1};
    auto const xyz_orig = geom->Plane(planeid).Wire(0).toWorldCoords(xyzWire);
    auto const dxyz_orig = geom->Plane(planeid).Wire(0).toWorldCoords(dxyzWire);

    mf::LogVerbatim("GeometryTest")
      << "\n\t" << xyz_orig.X() << "\t" << xyz_orig.Y() << "\t" << xyz_orig.Z();
    mf::LogVerbatim("GeometryTest")
      << "\t" << dxyz_orig.X() << "\t" << dxyz_orig.Y() << "\t" << dxyz_orig.Z();

    double xyz[3];
    double dxyz[3];
    xyz_orig.GetCoordinates(xyz);
    dxyz_orig.GetCoordinates(dxyz);
    gGeoManager->InitTrack(xyz, dxyz);
    for (int i = 0; i < 10; ++i) {
      const double* pos = gGeoManager->GetCurrentPoint();
      const double* dir = gGeoManager->GetCurrentDirection();
      mf::LogVerbatim("GeometryTest")
        << "\tnode = " << gGeoManager->GetCurrentNode()->GetName() << "\n\t\tpos="
        << "\t" << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\n\t\tdir="
        << "\t" << dir[0] << "\t" << dir[1] << "\t" << dir[2]
        << "\n\t\tmat = " << gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->GetName();

      gGeoManager->FindNextBoundary();
      gGeoManager->FindNormal();
      gGeoManager->Step(kTRUE, kTRUE);
    }

    xyz[0] = 306.108;
    xyz[1] = -7.23775;
    xyz[2] = 856.757;
    gGeoManager->InitTrack(xyz, dxyz);
    mf::LogVerbatim("GeometryTest")
      << "\tnode = " << gGeoManager->GetCurrentNode()->GetName()
      << "\n\tmat = " << gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->GetName();

    gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->Print();
  }

  //......................................................................

  void GeometryTestAlg::testProject()
  {
    double xlo, xhi;
    double ylo, yhi;
    double zlo, zhi;
    geom->WorldBox(&xlo, &xhi, &ylo, &yhi, &zlo, &zhi);

    double xyz[3] = {0.0, 0.0, 0.0};
    double dxyz1[3] = {1.0, 0.0, 0.0};
    double dxyz2[3] = {-1.0, 0.0, 0.0};
    double dxyz3[3] = {0.0, 1.0, 0.0};
    double dxyz4[3] = {0.0, -1.0, 0.0};
    double dxyz5[3] = {0.0, 0.0, 1.0};
    double dxyz6[3] = {0.0, 0.0, -1.0};

    double xyzo[3];
    geo::ProjectToBoxEdge(xyz, dxyz1, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[0] - xhi) > 1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz2, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[0] - xlo) > 1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz3, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[1] - yhi) > 1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz4, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[1] - ylo) > 1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz5, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[2] - zhi) > 1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz6, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (std::abs(xyzo[2] - zlo) > 1.E-6) abort();
  }

  //......................................................................
  bool GeometryTestAlg::CheckAuxDetAtPosition(geo::Point_t const pos, unsigned int expected) const
  {
    unsigned int foundDet = std::numeric_limits<unsigned int>::max();
    try {
      foundDet = geom->FindAuxDetAtPosition(pos);
    }
    catch (cet::exception const& e) {
      mf::LogProblem("GeometryTestAlg")
        << "Caught an exception while looking for aux det around " << lar::dump::array<3U>(pos)
        << " (within aux det #" << expected << "); message:\n"
        << e.what();
      return false;
    }
    if (foundDet != expected) {
      mf::LogProblem("GeometryTestAlg")
        << "Auxiliary detector at position " << lar::dump::array<3U>(pos)
        << ", expected within aux det #" << expected << ", was returned to be " << foundDet
        << " instead";
      return false;
    }
    return true;
  } // GeometryTestAlg::CheckAuxDetAtPosition()

  //......................................................................
  bool GeometryTestAlg::CheckAuxDetSensitiveAtPosition(geo::Point_t const pos,
                                                       unsigned int expectedDet,
                                                       unsigned int expectedSens) const
  {
    size_t foundDet = std::numeric_limits<unsigned int>::max();
    size_t foundSensDet = std::numeric_limits<unsigned int>::max();
    try {
      geom->FindAuxDetSensitiveAtPosition(pos, foundDet, foundSensDet);
    }
    catch (cet::exception const& e) {
      mf::LogProblem("GeometryTestAlg")
        << "Caught an exception while looking for aux det sensitive around "
        << lar::dump::array<3U>(pos) << " (within aux det #" << expectedDet << " sensitive volume #"
        << expectedSens << "); message:\n"
        << e.what();
      return false;
    }
    if ((foundDet != expectedDet) || (foundSensDet != expectedSens)) {
      mf::LogProblem("GeometryTestAlg")
        << "Auxiliary detector at position " << lar::dump::array<3U>(pos)
        << ", expected within aux det #" << expectedDet << ", sensitive volume #" << expectedSens
        << ", was returned to be in aux det #" << foundDet << " sensitive volume #" << foundSensDet
        << " instead";
      return false;
    }
    return true;
  } // GeometryTestAlg::CheckAuxDetAtPosition()

  //......................................................................
  void GeometryTestAlg::testFindAuxDet() const
  {

    /*
     *
     * Picks the center of each sensitive detector and verifies that the
     * correct sensitive detector and auxiliary detector are found.
     *
     */

    unsigned int nErrors = 0;

    unsigned int const nAuxDets = geom->NAuxDets();

    for (unsigned int iDet = 0; iDet < nAuxDets; ++iDet) {

      geo::AuxDetGeo const& auxDet = geom->AuxDet(iDet);
      unsigned int const nSensitive = auxDet.NSensitiveVolume();

      if (nSensitive == 0) {
        auto const center = auxDet.GetCenter();

        if (!CheckAuxDetAtPosition(center, iDet)) ++nErrors;
      }
      else { // if one or more sensitive detectors

        for (unsigned int iDetSens = 0; iDetSens < nSensitive; ++iDetSens) {

          geo::AuxDetSensitiveGeo const& auxDetSens = auxDet.SensitiveVolume(iDetSens);
          auto const center = auxDetSens.GetCenter();

          if (!CheckAuxDetAtPosition(center, iDet)) ++nErrors;
          if (!CheckAuxDetSensitiveAtPosition(center, iDet, iDetSens)) ++nErrors;

        } // for sensitive detectors

      } // if ... else

    } // for all auxiliary detectors

    if (nErrors != 0) {
      throw cet::exception("FindAuxDet")
        << "Collected " << nErrors << " errors during testFindAuxDet() test!\n";
    }

  } // GeometryTestAlg::testFindAuxDet()

  //......................................................................

  inline bool GeometryTestAlg::shouldRunTests(std::string test_name) const
  {
    return fRunTests(test_name);
  } // GeometryTestAlg::shouldRunTests()

} //end namespace
