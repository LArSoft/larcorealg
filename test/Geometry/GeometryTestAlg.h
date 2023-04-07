/**
 * @file   GeometryTestAlg.h
 * @brief  Unit test for geometry functionalities
 * @date   2011/02/17
 * @author brebel@fnal.gov
 * @see    GeometryTestAlg.cxx
 *
 * Refactored by Gianluca Petrillo on May 5th, 2015.
 */

#ifndef GEO_GEOMETRYTESTALG_H
#define GEO_GEOMETRYTESTALG_H

// LArSoft includes
#include "larcorealg/Geometry/fwd.h"
#include "larcorealg/TestUtils/NameSelector.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "fhiclcpp/fwd.h"

// C/C++ standard libraries
#include <array>
#include <set>
#include <string>
#include <vector>

namespace geo {

  // forward declarations
  struct CryostatID;
  struct PlaneID;

  /** **************************************************************************
   * @brief Performs tests on the geometry as seen by Geometry service
   *
   * Configuration parameters
   * =========================
   *
   * - **DisableWireBoundaryCheck** (boolean, default: false): the exceptions
   *   thrown when checking wire boundaries are not fatal
   * - **ExpectedWirePitches** (list of reals, default: empty):
   *   if specified, marks the expected uniform wire pitch, one entry for each
   *   plane; if there are more planes than elements, the missing pitches are
   *   assumed the same as the one in the last plane; if the parameter is
   *   not specified or empty, no expectation is taken, but the spacing is still
   *   checked to be uniform; exception: for legacy behaviour, default pitches
   *   are provided for some experiments
   * - **ExpectedPlanePitches** (list of reals, default: empty):
   *   if specified, marks the expected uniform plane pitch, one entry for each
   *   plane vs. the following. It works as /ExpectedWirePitches/
   * - **ForgiveExceptions** (list of strings, default: empty): the categories
   *   of exceptions in this list are "forgiven" (non-fatal)
   * - **RunTests** (string list): marks which tests to run; each entry can be
   *   specified with a "+" or "-" prepended, to indicate to add or to remove
   *   the test from the set to be run, respectively. The directives act against
   *   a pre-existing default set that executes all tests, except the ones
   *   explicitly marked as excluded by default in the list below:
   *   + `CheckOverlaps` (not in default) perform overlap checks
   *   + `ThoroughCheck` (not in default) makes ROOT perform full geometry check
   *   + `DetectorIntro`: prints some information about the detector
   *   + `FindVolumes`: checks it can find the volumes corresponding to world
   *     and all cryostats
   *   + `Cryostat`:
   *   + `WireOrientations`: checks that the definition of wire coordinates is
   *     matching the prescription
   *   + `WireCoordFromPlane`: checks `PlaneGeo::WireCoordinateFrom()`
   *   + `ChannelToWire`:
   *   + `FindPlaneCenters`:
   *   + `Projection`:
   *   + `WirePos`: currently disabled
   *   + `PlanePointDecomposition`: methods for projections and decompositions
   *     on the wire coordinate reference system
   *   + `PlaneProjections`: methods for projections on the wire planes in the
   *     reference system of the frame of the plane
   *   + `WireCoordAngle`: tests geo::PlaneGeo::PhiZ()
   *   + `NearestWire`: tests `WireCoordinate()` and `NearestWire()`
   *   + `WireIntersection`: tests `WireIDsIntersect()`
   *   + `ThirdPlane`: tests `ThirdPlane()`
   *   + `ThirdPlaneSlope`: tests `ThirdPlaneSlope()`
   *   + `WirePitch`:
   *   + `PlanePitch`:
   *   + `InterWireProjectedDistance`: tests `geo::PlaneGeo::InterWireProjectedDistance()`
   *   + `Stepping`:
   *   + `FindAuxDet`: test on location of nearest auxiliary detector
   *   + `PrintWires`: (not in default) prints *all* the wires in the geometry
   *   + `default`: represents the default set (optionally prepended by '@')
   *   + `!` (special): means to forget the tests configured so far; used as the
   *     first test name, removes the default list but leaves unchanged the
   *     default behaviour (the one specified with "+*" or "-*")
   * - **CheckForOverlaps** (boolean, default: false): equivalent to enabling
   *   `+CheckOverlaps` in `RunTests`
   * - **PrintWires**: (boolean, default: false): equivalent to enabling
   *   `+PrintWires` in `RunTests`
   * - **ComputeMass**: (boolean, default: `true`): prints the mass of the
   *     detector and cryostats, which may take a lot of time.
   */
  class GeometryTestAlg {
  public:
    explicit GeometryTestAlg(GeometryCore const* geom,
                             WireReadoutGeom const* wireReadoutGeom,
                             fhicl::ParameterSet const& pset);

    /// Virtual destructor
    virtual ~GeometryTestAlg() = default;

    /// Runs the test, returns a number of errors (very unlikely!)
    virtual unsigned int Run();

    /// Returns the direction on plane orthogonal to wires where wire number increases
    static std::array<double, 3> GetIncreasingWireDirection(PlaneGeo const& plane);

  private:
    GeometryCore const* geom;
    WireReadoutGeom const* wireReadoutGeom;

    bool fDisableValidWireIDcheck; ///< disable test on out-of-world NearestWire()
    std::set<std::string> fNonFatalExceptions;
    std::vector<double> fExpectedWirePitches;  ///< wire pitch on each plane
    std::vector<double> fExpectedPlanePitches; ///< plane pitch on each plane

    bool fComputeMass = true; ///< Whether to print the detector mass.

    // using as pointer just not to have to write the declaration in the header
    testing::NameSelector fRunTests; ///< test filter

    void printDetectorIntro() const;
    void printVolBounds();
    void printDetDim();
    void printWirePos();
    void printWiresInTPC(TPCGeo const& tpc, std::string indent = "") const;
    void printAllGeometry() const;
    void testFindVolumes();
    void testCryostat();
    void testTPC(CryostatID const& cid);
    void testPlaneDirections() const;
    void testWireOrientations() const;
    void testChannelToROP() const;
    void testChannelToWire() const;
    void testFindPlaneCenters();
    void testProject();
    void testPlaneProjectionReference() const;
    void testPlanePointDecompositionFrame() const;
    void testPlaneProjectionOnFrame() const;
    void testPlaneProjection() const;
    void testWireCoordFromPlane() const;
    void testParallelWires() const;
    void testPlanePointDecomposition() const;
    void testWireCoordAngle() const;
    void testWirePitch();
    void testInterWireProjectedDistance() const;
    void testPlanePitch();
    void testStandardWirePos();
    void testAPAWirePos();
    void testNearestWire();
    void testWireIntersection() const;
    void testThirdPlane() const;
    void testThirdPlane_dTdW() const;
    void testStepping();
    void testFindAuxDet() const;

    bool shouldRunTests(std::string test_name) const;

    // single tests for testFindVolumes
    unsigned int testFindWorldVolumes();
    unsigned int testFindCryostatVolumes();
    unsigned int testFindTPCvolumePaths();

    /// Method to print the auxiliary detectors on screen.
    void printAuxiliaryDetectors() const;

    /// Prints information of an auxiliary detector into the specified stream.
    template <typename Stream>
    void printAuxDetGeo(Stream&& out,
                        AuxDetGeo const& auxDet,
                        std::string indent,
                        std::string firstIndent) const;

    /// Prints information of an auxiliary detector into the specified stream.
    template <typename Stream>
    void printAuxDetGeo(Stream&& out, AuxDetGeo const& auxDet, std::string indent = "") const
    {
      printAuxDetGeo(std::forward<Stream>(out), auxDet, indent, indent);
    }

    /// Prints information of the sensitive auxiliary detector into a stream.
    template <typename Stream>
    void printAuxDetSensitiveGeo(Stream&& out,
                                 AuxDetSensitiveGeo const& auxDetSens,
                                 std::string indent,
                                 std::string firstIndent) const;

    /// Prints information of the sensitive auxiliary detector into a stream.
    template <typename Stream>
    void printAuxDetSensitiveGeo(Stream&& out,
                                 AuxDetSensitiveGeo const& auxDetSens,
                                 std::string indent = "") const
    {
      printAuxDetSensitiveGeo(std::forward<Stream>(out), auxDetSens, indent, indent);
    }

    /// Returns whether the auxiliary detector at `pos` is the `expected` one.
    bool CheckAuxDetAtPosition(Point_t const pos, unsigned int expected) const;

    /// Returns whether the auxiliary sensitive detector at `pos` is expected.
    bool CheckAuxDetSensitiveAtPosition(Point_t const pos,
                                        unsigned int expectedDet,
                                        unsigned int expectedSens) const;

    /// Helper function for `testWireIntersection()`.
    bool isWireAlignedToPlaneDirections(PlaneGeo const& plane, Vector_t const& wireDir) const;

    /// Performs the wire intersection test at a single point
    unsigned int testWireIntersectionAt(TPCGeo const& TPC, Point_t const& point) const;

    /// Returns dT/dW expected from the specified segment A-to-B
    std::vector<std::pair<PlaneID, double>>
    ExpectedPlane_dTdW(Point_t const& A, Point_t const& B, double const driftVelocity = -0.1) const;

    /// Performs the third plane slope test with a single configuration
    unsigned int testThirdPlane_dTdW_at(
      std::vector<std::pair<PlaneID, double>> const& plane_dTdW) const;
  };

  namespace details {
    /// Class telling whether a test needs to be run
    class TestTrackerClassBase {
    public:
      using TestList_t = std::set<std::string>;

      virtual ~TestTrackerClassBase() = default;

      /// Returns whether the specified test should run
      virtual bool ShouldRun(std::string test_name) const = 0;

      /// Checks the test and records the request
      bool operator()(std::string test_name);

      /// Allow the specified test to run
      virtual void PleaseRunAlso(std::string test_name) = 0;

      /// Returns the tests that have been run
      TestList_t const& RunTests() const { return run; }

      /// Returns the tests that have been skipped
      TestList_t const& SkippedTests() const { return skipped; }

      /// Returns the tests that have been queried
      TestList_t QueriedTests() const;

      /// Checks that the validity of the configuration (after the fact)
      virtual bool CheckQueriesRegistry() const;

      /// Prints information about the configuration of the filter
      virtual void PrintConfiguration(std::ostream&) const;

    protected:
      TestList_t run;     ///< requested tests that should be run
      TestList_t skipped; ///< requested tests that should be skipped

      virtual void RecordRequest(std::string test_name, bool bRun);

      /// Checks the test and records the request
      virtual bool Query(std::string test_name);

      /// Adds a vector of tests into a test set
      static void CopyList(TestList_t& dest, std::vector<std::string> const& from);
    }; // class TestTrackerClassBase

  } // namespace details

} // namespace geo

#endif // GEO_GEOMETRYTESTALG_H
