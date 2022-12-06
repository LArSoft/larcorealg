/**
 * @file   larcorealg/Geometry/GeometryCore.h
 * @brief  Access the description of detector geometry
 * @author brebel@fnal.gov
 * @see    larcorealg/Geometry/GeometryCore.cxx
 * @ingroup Geometry
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYCORE_H
#define LARCOREALG_GEOMETRY_GEOMETRYCORE_H

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/Iterable.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/details/GeometryIterationPolicy.h"
#include "larcorealg/Geometry/fwd.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"       // geo::vect namespace
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

// Framework and infrastructure libraries
#include "fhiclcpp/ParameterSet.h"

// External libraries
#include "range/v3/view.hpp"

// C/C++ standard libraries
#include <cstddef>  // size_t
#include <iterator> // std::forward_iterator_tag
#include <memory>   // std::shared_ptr<>
#include <optional>
#include <set>
#include <string>
#include <type_traits> // std::is_base_of<>
#include <utility>
#include <vector>

// ROOT class prototypes
class TGeoManager;
class TGeoNode;
class TGeoVolume;
class TGeoMaterial;

/// Namespace collecting geometry-related classes utilities
namespace geo {

  // BEGIN Geometry group ------------------------------------------------------
  /// @ingroup Geometry
  /// @{

  //
  // GeometryCore
  //

  /** **************************************************************************
   * @brief Description of geometry of one entire detector
   *
   * @note All lengths are specified in centimeters
   *
   *
   * How to correctly instantiate a GeometryCore object
   * ---------------------------------------------------
   *
   * Instantiation is a multi-step procedure:
   *
   * 1. construct a GeometryCore object (the "service provider"), with the full
   *    configuration; at this step, configuration is just stored
   *
   * 2. load a geometry with GeometryCore::LoadGeometryFile(); this loads the detector
   *    geometry information
   *
   * 3. prepare a channel map algorithm object (might use for example
   *    GeometryCore::DetectorName() or the detector geometry from the newly created
   *    object, but any use of channel mapping related functions is forbidden and it would
   *    yield undefined behaviour (expected to be catastrophic)
   *
   * 4. sort geometry according to the sorter provided by the channel map
   *
   * Step 3 (creation of the channel mapping algorithm object) can be performed at any
   * time before step 4, provided that no GeometryCore instance is needed for it.
   *
   *
   * Configuration parameters
   * ------------------------
   *
   * - *Name* (string; mandatory): string identifying the detector; it can be different
   *   from the base name of the file used to initialize the geometry; standard names are
   *   recommended by each experiment.  This name can be used, for example, to select
   *   which channel mapping algorithm to use.
   * - *SurfaceY* (real; mandatory): depth of the detector, in centimeters; see SurfaceY()
   *   for details
   * - *MinWireZDist* (real; default: 3)
   * - *PositionEpsilon* (real; default: 0.01%) set the default tolerance (see
   *   DefaultWiggle())
   *
   */
  class ToGeometryElement {
  public:
    explicit ToGeometryElement(GeometryCore const* geom) : fGeom{geom} {}
    template <typename T, typename Iterator>
    auto transform(Iterator const& iterator) const
    {
      return details::geometry_element_iterator<T, Iterator>{fGeom, iterator};
    }

  private:
    GeometryCore const* fGeom;
  };

  class GeometryCore : Iterable<details::GeometryIterationPolicy, ToGeometryElement> {
    using Iteration = Iterable<details::GeometryIterationPolicy, ToGeometryElement>;

  public:
    /// Simple class with two points (a pair with aliases).
    template <typename Point>
    struct Segment : public std::pair<Point, Point> {

      // use the base class constructors
      using std::pair<Point, Point>::pair;

      Point const& start() const { return this->first; }
      Point& start() { return this->first; }

      Point const& end() const { return this->second; }
      Point& end() { return this->second; }

    }; // struct Segment_t

    using Segment_t = Segment<Point_t>;

    /// Type of list of cryostats
    using CryostatList_t = std::vector<geo::CryostatGeo>;
    /// Type of list of auxiliary detectors
    using AuxDetList_t = std::vector<geo::AuxDetGeo>;

    /**
     * @brief Initialize geometry from a given configuration
     * @param pset configuration parameters
     *
     * This constructor does not load any geometry description.  The next step is to do
     * exactly that, by GeometryCore::LoadGeometryFile().
     */
    GeometryCore(fhicl::ParameterSet const& pset, std::unique_ptr<GeoObjectSorter const> sorter);
    GeometryCore(GeometryCore const&) = delete;
    GeometryCore(GeometryCore&&) = delete;
    GeometryCore& operator=(GeometryCore const&) = delete;
    GeometryCore& operator=(GeometryCore&&) = delete;

    /**
     * @brief Returns the tolerance used in looking for positions
     * @return the tolerance value
     *
     * This parameter is used as tolerance ("wiggle") for methods that require it
     * (e.g. `geo::CryostatGeo::FindTPCAtPosition()`).  Typically, it's a additional
     * fraction of tolerance: 0 means no tolerance, 0.1 means 10% tolerance.
     *
     * @todo Confirm the definition of wiggle: this one is taken from other doc
     */
    double DefaultWiggle() const { return fPositionWiggle; }

    /**
     * @brief Returns the full directory path to the geometry file source
     * @return the full directory path to the geometry file source
     *
     * This is the full path of the source of the detector geometry GeometryCore relies
     * on.
     */
    std::string const& ROOTFile() const { return fROOTfile; }

    /**
     * @brief Returns the full directory path to the GDML file source
     * @return the full directory path to the GDML file source
     *
     * This is the full path of the source of the detector geometry handed to the detector
     * simulation (GEANT).
     */
    std::string const& GDMLFile() const { return fGDMLfile; }

    // BEGIN Detector information
    /// @name Detector information
    /// @{

    //
    // global features
    //
    /// Returns a string with the name of the detector, as configured
    std::string const& DetectorName() const { return fDetectorName; }

    //
    // position
    //

    /// Returns a pointer to the world volume.
    TGeoVolume const* WorldVolume() const;

    /**
     * @brief Fills the arguments with the boundaries of the world
     * @param xlo (output) pointer to the lower x coordinate
     * @param xlo (output) pointer to the upper x coordinate
     * @param ylo (output) pointer to the lower y coordinate
     * @param ylo (output) pointer to the upper y coordinate
     * @param zlo (output) pointer to the lower z coordinate
     * @param zlo (output) pointer to the upper z coordinate
     * @throw cet::exception (`"GeometryCore"` category) if no world found
     * @see `GetWorldVolumeName()`
     *
     * This method fills the boundaries of the world volume
     * (`GetWorldVolumeName()`).
     *
     * If a pointer is null, its coordinate is skipped.
     *
     * @deprecated Use the version without arguments instead.
     */
    void WorldBox(double* xlo, double* xhi, double* ylo, double* yhi, double* zlo, double* zhi)
      const;

    /// Returns a box with the extremes of the world volume (from shape axes).
    /// @see `GetWorldVolumeName()`
    BoxBoundedGeo WorldBox() const;

    /**
     * @brief The position of the detector respect to earth surface
     * @return typical y position at surface in units of cm
     *
     * This is the depth (y) of the surface (where earth meets air) for this detector
     * site.  The number is expressed in world coordinates and in centimeters, and it
     * represents the y coordinate of earth surface.  A negative value means that the
     * origin of coordinates, typically matching the detector centre, is above surface.
     *
     * @todo check that this is actually how it is used
     */
    //
    Length_t SurfaceY() const { return fSurfaceY; }

    //
    // object description and information
    //

    /// Access to the ROOT geometry description manager
    TGeoManager* ROOTGeoManager() const;

    /// Return the name of the world volume (needed by Geant4 simulation)
    std::string const& GetWorldVolumeName() const;

    /// Returns the absolute  coordinates of the detector enclosure volume [cm].
    /// @param name name of the volume to be sought (default: `volDetEnclosure`)
    /// @throw cet::exception if the specified volume is not found
    BoxBoundedGeo DetectorEnclosureBox(std::string const& name = "volDetEnclosure") const;

    //@{
    /**
     * @brief Returns the name of the deepest volume containing specified point
     * @param point the location to query, in world coordinates
     * @return name of the volume containing the point
     *
     * @todo what happens if none?
     * @todo Unify the coordinates type
     */
    std::string VolumeName(Point_t const& point) const;
    //@}

    /**
     * @brief Returns all the nodes with volumes with any of the specified names
     * @param vol_names list of names of volumes
     * @return list of nodes found
     *
     * All the nodes in the geometry are checked, and all the ones that contain a volume
     * with a name among the ones specified in vol_names are saved in the collection and
     * returned.
     */
    std::vector<TGeoNode const*> FindAllVolumes(std::set<std::string> const& vol_names) const;

    /**
     * @brief Returns paths of all nodes with volumes with the specified names
     * @param vol_names list of names of volumes
     * @return list paths of the found nodes
     *
     * All the nodes in the geometry are checked, and the path of all the ones that
     * contain a volume with a name among the ones specified in vol_names is saved in the
     * collection and returned.
     *
     * A node path is a ordered list of all nodes leading to the final one, starting from
     * thetop level (root) down. The node at the `back()` of the path is the one with name
     * in vol_names.  No empty paths are returned.
     */
    std::vector<std::vector<TGeoNode const*>> FindAllVolumePaths(
      std::set<std::string> const& vol_names) const;

    /// Returns the material at the specified position
    TGeoMaterial const* Material(Point_t const& point) const;
    //@{
    /**
     * @brief Name of the deepest material containing the point xyz
     * @return material of the origin by default
     */
    std::string MaterialName(Point_t const& point) const;
    //@}

    //@{
    /// Returns the total mass [kg] of the specified volume (default: world).
    double TotalMass() const { return TotalMass(GetWorldVolumeName()); }
    double TotalMass(std::string vol) const;
    //@}

    //@{
    /**
     * @brief Returns the column density between two points.
     * @param p1 the first point
     * @param p2 the second point
     * @return the column density [kg / cm&sup2;]
     *
     * The column density is defined as
     * @f$ \int_{\vec{p}_{1}}^{\vec{p}_{2}} \rho(\vec{p}) d\vec{p} @f$
     * where @f$ \rho(\vec{p}) @f$ is the density at point @f$ \vec{p} @f$,
     * which the integral leads from `p1` to `p2` in a straight line.
     *
     * Both points are specified in world coordinates.
     */
    double MassBetweenPoints(Point_t const& p1, Point_t const& p2) const;
    //@}

    /// Prints geometry information with maximum verbosity.
    template <typename Stream>
    void Print(Stream&& out, std::string indent = "  ") const;

    /// @brief Returns a string with complete geometry information.
    /// @see `Print()`
    std::string Info(std::string indent = "  ") const;

    /// @}
    // END Detector information

    /// @name Cryostat access and information
    /// @{

    //
    // group features
    //

    //@{
    /**
     * @brief Returns the number of cryostats in the detector
     *
     * The NElements() and NSiblingElements() methods are overloaded and their return
     * depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Ncryostats() const { return Cryostats().size(); }
    unsigned int NElements() const { return Ncryostats(); }
    unsigned int NSiblingElements(CryostatID const&) const { return Ncryostats(); }
    //@}

    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified cryostat
     *
     * The HasElement() method is overloaded and its meaning depends on the type of ID.
     */
    bool HasCryostat(CryostatID const& cryoid) const { return cryoid.Cryostat < Ncryostats(); }
    bool HasElement(CryostatID const& cryoid) const { return HasCryostat(cryoid); }
    //@}

    //@{
    /**
     * @brief Returns the specified cryostat
     * @param cstat number of cryostat
     * @param cryoid cryostat ID
     * @return a constant reference to the specified cryostat
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     *
     * The GetElement() method is overloaded and its return depends on the type of ID.
     *
     * @todo Make the cryostat number mandatory (as CryostatID)
     */
    static constexpr CryostatID cryostat_zero{0};
    CryostatGeo const& Cryostat(CryostatID const& cryoid = cryostat_zero) const;
    CryostatGeo const& GetElement(CryostatID const& cryoid) const { return Cryostat(cryoid); }
    //@}

    //@{
    /**
     * @brief Returns the specified cryostat
     * @param cryoid cryostat ID
     * @return a constant pointer to the specified cryostat, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the type of ID.
     */
    CryostatGeo const* CryostatPtr(CryostatID const& cryoid) const
    {
      return HasCryostat(cryoid) ? &(Cryostats()[cryoid.Cryostat]) : nullptr;
    }
    CryostatGeo const* GetElementPtr(CryostatID const& cryoid) const { return CryostatPtr(cryoid); }
    //@}

    /**
     * @brief Returns the cryostat at specified location.
     * @param point the location [cm]
     * @return pointer to the `geo::CryostatGeo` including `point`, or `nullptr`
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    CryostatGeo const* PositionToCryostatPtr(Point_t const& point) const;

    /**
     * @brief Returns the ID of the cryostat at specified location.
     * @param point the location [cm]
     * @return ID of the cryostat including `point` (invalid if none)
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    CryostatID PositionToCryostatID(Point_t const& point) const;

    //@{
    /**
     * @brief Returns the cryostat at specified location.
     * @param point the location [cm]
     * @return a constant reference to the `geo::CryostatGeo` containing `point`
     * @throws cet::exception ("Geometry" category) if no cryostat matches
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    CryostatGeo const& PositionToCryostat(Point_t const& point) const;

    //
    // iterators
    //
    using Iteration::begin;
    using Iteration::end;
    using Iteration::Iterate;

    //
    // single object features
    //

    /// @name TPC access and information
    /// @{

    //
    // group features
    //

    /// Returns the largest number of TPCs a cryostat in the detector has
    unsigned int MaxTPCs() const;

    /// Returns the total number of TPCs in the detector
    unsigned int TotalNTPC() const;

    //@{
    /**
     * @brief Returns the total number of TPCs in the specified cryostat
     * @param cryoid cryostat number
     * @return number of TPCs in specified cryostat, or 0 if no cryostat found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their return
     * depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int NTPC(CryostatID const& cryoid = cryostat_zero) const
    {
      CryostatGeo const* pCryo = GetElementPtr(cryoid);
      return pCryo ? pCryo->NElements() : 0;
    }
    unsigned int NElements(CryostatID const& cryoid) const { return NTPC(cryoid); }
    unsigned int NSiblingElements(TPCID const& tpcid) const { return NTPC(tpcid); }
    //@}

    //
    // access
    //
    /// Returns whether we have the specified TPC
    bool HasTPC(TPCID const& tpcid) const
    {
      CryostatGeo const* pCryo = CryostatPtr(tpcid);
      return pCryo ? pCryo->HasTPC(tpcid) : false;
    }

    /// Returns whether we have the specified TPC
    bool HasElement(TPCID const& tpcid) const { return HasTPC(tpcid); }

    //@{
    /**
     * @brief Returns the specified TPC
     * @param tpcid ID of the tpc
     * @param tpc tpc number within the cryostat
     * @param cstat number of cryostat
     * @return a constant reference to the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     *
     * The GetElement() method is overloaded and its return depends on the type of ID.
     */
    static constexpr TPCID tpc_zero{cryostat_zero, 0};
    TPCGeo const& TPC(TPCID const& tpcid = tpc_zero) const { return Cryostat(tpcid).TPC(tpcid); }
    TPCGeo const& GetElement(TPCID const& tpcid) const { return TPC(tpcid); }
    //@}

    //@{
    /**
     * @brief Returns the specified TPC
     * @param tpcid TPC ID
     * @return a constant pointer to the specified TPC, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the type of ID.
     */
    TPCGeo const* TPCPtr(TPCID const& tpcid) const
    {
      CryostatGeo const* pCryo = CryostatPtr(tpcid);
      return pCryo ? pCryo->TPCPtr(tpcid) : nullptr;
    }
    TPCGeo const* GetElementPtr(TPCID const& tpcid) const { return TPCPtr(tpcid); }
    //@}

    //@{
    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param worldLoc 3D point (world reference frame, centimeters)
     * @return the TPC ID, or an invalid one if no TPC is there
     */
    TPCID FindTPCAtPosition(Point_t const& point) const;
    //@}

    /**
     * @brief Returns the TPC at specified location.
     * @param point the location [cm]
     * @return the `geo::TPCGeo` including `point`, or `nullptr` if none
     */
    TPCGeo const* PositionToTPCptr(Point_t const& point) const;

    //@{
    /**
     * @brief Returns the TPC at specified location.
     * @param point the location [cm]
     * @return a constant reference to the `geo::TPCGeo` including `point`
     * @throws cet::exception ("Geometry" category) if no TPC matches
     */
    TPCGeo const& PositionToTPC(Point_t const& point) const;
    //@}

    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param point the location [cm]
     * @return ID of the TPC at specified location, invalid if none
     * @see `PositionToTPC()`
     */
    TPCID PositionToTPCID(Point_t const& point) const;

    //
    // object description
    //

    /// @name Plane access and information
    /// @{

    //
    // group features
    //

    /// Returns the largest number of planes among all TPCs in this detector
    unsigned int MaxPlanes() const;

    //@{
    /**
     * @brief Returns the total number of planes in the specified TPC
     * @param tpcid TPC ID
     * @return number of planes in specified TPC, or 0 if no TPC found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their return
     * depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nplanes(TPCID const& tpcid = tpc_zero) const
    {
      TPCGeo const* pTPC = GetElementPtr(tpcid);
      return pTPC ? pTPC->NElements() : 0;
    }
    unsigned int NElements(TPCID const& tpcid) const { return Nplanes(tpcid); }
    unsigned int NSiblingElements(PlaneID const& planeid) const { return Nplanes(planeid); }
    //@}

    /**
     * @brief Returns the number of views (different wire orientations)
     *
     * Returns the number of different views, or wire orientations, in the detector.
     *
     * The function assumes that all TPCs in all cryostats of a detector have the same
     * number of planes, which should be a safe assumption.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nviews() const;

    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified plane
     *
     * The HasElement() method is overloaded and its meaning depends on the type of ID.
     *
     */
    bool HasPlane(PlaneID const& planeid) const
    {
      TPCGeo const* pTPC = TPCPtr(planeid);
      return pTPC ? pTPC->HasPlane(planeid) : false;
    }
    bool HasElement(PlaneID const& planeid) const { return HasPlane(planeid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param planeid ID of the plane
     * @param p plane number within the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return a constant reference to the specified plane
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @throw cet::exception (`PlaneOutOfRange` category) if no such plane
     *
     * The GetElement() method is overloaded and its return depends on the type of ID.
     */
    PlaneGeo const& Plane(PlaneID const& planeid) const { return TPC(planeid).Plane(planeid); }
    PlaneGeo const& GetElement(PlaneID const& planeid) const { return Plane(planeid); }
    //@}

    //@{
    /**
     * @brief Returns the specified plane
     * @param planeid plane ID
     * @return a constant pointer to the specified plane, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the type of ID.
     */
    PlaneGeo const* PlanePtr(PlaneID const& planeid) const
    {
      TPCGeo const* pTPC = TPCPtr(planeid);
      return pTPC ? pTPC->PlanePtr(planeid) : nullptr;
    }
    PlaneGeo const* GetElementPtr(PlaneID const& planeid) const { return PlanePtr(planeid); }
    //@}

    /// @} Plane access and information

    /// @name Wire access and information
    /// @{

    //
    // group features
    //

    //@{
    /**
     * @brief Returns the total number of wires in the specified plane
     * @param planeid plane ID
     * @return number of wires in specified plane, or 0 if no plane found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their return
     * depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nwires(PlaneID const& planeid) const
    {
      PlaneGeo const* pPlane = GetElementPtr(planeid);
      return pPlane ? pPlane->NElements() : 0;
    }
    unsigned int NElements(PlaneID const& planeid) const { return Nwires(planeid); }
    unsigned int NSiblingElements(WireID const& wireid) const { return Nwires(wireid); }

    /// Returns the largest number of wires among all planes in this detector
    unsigned int MaxWires() const;

    //@}

    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified wire
     *
     * The HasElement() method is overloaded and its meaning depends on the type of ID.
     */
    bool HasWire(WireID const& wireid) const
    {
      PlaneGeo const* pPlane = PlanePtr(wireid);
      return pPlane ? pPlane->HasWire(wireid) : false;
    }
    bool HasElement(WireID const& wireid) const { return HasWire(wireid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param wireid wire ID
     * @return a constant pointer to the specified wire, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the type of ID.
     */
    WireGeo const* WirePtr(WireID const& wireid) const
    {
      PlaneGeo const* pPlane = PlanePtr(wireid);
      return pPlane ? pPlane->WirePtr(wireid) : nullptr;
    } // WirePtr()
    WireGeo const* GetElementPtr(WireID const& wireid) const { return WirePtr(wireid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param wireid ID of the wire
     * @return a constant reference to the specified wire
     * @throw cet::exception if not found
     *
     * The GetElement() method is overloaded and its return depends on the type of ID.
     */
    WireGeo const& Wire(WireID const& wireid) const { return Plane(wireid).Wire(wireid); }
    WireGeo const& GetElement(WireID const& wireid) const { return Wire(wireid); }
    //@}

    //
    // single object features
    //

    //@{
    /**
     * @brief Returns the angle of the wires in the specified view from vertical
     * @param view the view
     * @param TPC the index of the TPC in the specified cryostat
     * @param Cryo the cryostat
     * @param tpcid ID of the TPC
     * @return the angle [radians]
     * @throw cet::exception ("GeometryCore" category) if no such view
     *
     * The angle is defined as in WireGeo::ThetaZ().
     *
     * This method assumes all wires in the view have the same angle (it queries for the
     * first).
     *
     * @deprecated This does not feel APA-ready
     */
    double WireAngleToVertical(View_t view, TPCID const& tpcid) const;
    //@}

    /// @} Wire access and information

    /**
     * @name Wire geometry queries
     *
     * Please note the differences between functions:
     * ChannelsIntersect(), WireIDsIntersect() and IntersectionPoint()
     * all calculate wires intersection using the same equation.
     * ChannelsIntersect() and WireIdsIntersect() will return true
     * if the two wires cross, return false if they don't.
     * IntersectionPoint() does not check if the two wires cross.
     */
    /// @{

    //
    // simple geometry queries
    //

    /**
     * @brief Fills two arrays with the coordinates of the wire end points
     * @param wireid ID of the wire
     * @param xyzStart (output) an array with the start coordinate
     * @param xyzEnd (output) an array with the end coordinate
     * @throws cet::exception wire not present
     *
     * The starting point is the wire end with lower z coordinate.
     *
     * @deprecated use the wire ID interface instead (but note that it does not
     *             sort the ends)
     */
    void WireEndPoints(WireID const& wireid, double* xyzStart, double* xyzEnd) const;

    //@{
    /**
     * @brief Returns a segment whose ends are the wire end points
     * @param wireid ID of the wire
     * @return a segment whose ends are the wire end points
     * @throws cet::exception wire not present
     *
     * The start and end are assigned as returned from the geo::WireGeo object.
     * The rules for this assignment are documented in that class.
     *
     * @deprecated use the wire ID interface instead (but note that it does not
     *             sort the ends)
     */
    Segment<Point_t> WireEndPoints(WireID const& wireID) const;

    //@}

    //
    // wire intersections
    //

    //@{
    /**
     * @brief Computes the intersection between two wires.
     * @param wid1 ID of the first wire
     * @param wid2 ID of the other wire
     * @param[out] intersection the intersection point (global coordinates)
     * @return whether an intersection was found inside the TPC the wires belong
     * @see `geo::WiresIntersection()`, `geo::LineClosestPoint()`
     *
     * The wires identified by `wid1` and `wid2` are intersected, and the 3D intersection
     * point is written into the `intersection` parameter.  The "intersection" point is
     * actually the point belonging to the first wire (`wid2`) which is the closest (in
     * Euclidean 3D metric) to the second wire.
     *
     * The intersection is computed only if the wires belong to different planes of the
     * same TPC. If that is not the case (i.e. they belong to different TPC or cryostat,
     * or if they belong to the same plane), `false` is returned and `intersection` is set
     * with all components to infinity (`std::numeric_limits<>::infinity()`).
     *
     * When the intersection is computed, it is always stored in the `intersection` output
     * parameter. Return value is `true` if this intersection lies within the physical
     * boundaries first wire, while it is instead `false` if it lies on the extrapolation
     * of the wire direction, but not within the wire physical extension.
     *
     * To test that the result is not infinity (nor NaN), use
     * `geo::vect::isfinite(intersection)` etc.
     *
     * @note If `geo::WireGeo` objects are already available, using instead the free
     *       function `geo::WiresIntersection()` or the method
     *       `geo::WireGeo::IntersectionWith()` is faster (and _recommended_).  For purely
     *       geometric intersection, `geo::LineClosestPoint()` is also available.
     */
    bool WireIDsIntersect(WireID const& wid1, WireID const& wid2, Point_t& intersection) const;
    //@}

    //@{
    /**
     * @brief Computes the intersection between two wires.
     * @param wid1 ID of the first wire
     * @param wid2 ID of the other wire
     * @param widIntersect (output) the coordinate of the intersection point
     * @return whether an intersection was found within the TPC
     *
     * The "intersection" refers to the projection of the wires into the same
     * @f$ x = 0 @f$ plane.
     * Wires are assumed to have at most one intersection.
     * If wires are parallel, `widIntersect` will have the two components set to
     * infinity (`std::numeric_limits<>::infinity()`) and the TPC number set to
     * invalid (`geo::TPCID::InvalidID`). Also, `false` is returned.
     * If the intersection is outside the TPC, `false` is also returned, but the
     * `widIntersect` will contain the coordinates of that intersection. The TPC
     * number is still set to invalid, although the intersection _might_ belong
     * to a valid TPC somewhere else.
     *
     *
     * @deprecated This method uses arbitrary assumptions and should not be
     *             used. Use the interface returning a full vector instead.
     */
    std::optional<WireIDIntersection> WireIDsIntersect(WireID const& wid1,
                                                       WireID const& wid2) const;
    //@}

    /**
     * @brief Returns the plane that is not in the specified arguments
     * @param pid1 a plane
     * @param pid2 another plane
     * @return the ID to the third plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     * @throws cet::exception (category: "GeometryCore") if pid1 and pid2 match
     *
     * This function requires a geometry with exactly three planes.  If the two input
     * planes are not on the same TPC, the result is undefined.
     */
    PlaneID ThirdPlane(PlaneID const& pid1, PlaneID const& pid2) const;

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param pid2 ID of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @param output_plane ID of the plane on which to calculate the slope
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if input planes match
     *
     * Given a slope as projected in two planes, returns the slope as projected in the
     * specified output plane.  The slopes are defined in uniform units; they should be
     * computed as distance ratios (or tangent of a geometrical angle; the formula is
     * still valid using dt/dw directly in case of equal wire pitch in all planes and
     * uniform drift velocity.
     */
    double ThirdPlaneSlope(PlaneID const& pid1,
                           double slope1,
                           PlaneID const& pid2,
                           double slope2,
                           PlaneID const& output_plane) const;

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param pid2 ID of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a slope as projected in two planes, returns the slope as projected in the
     * third plane.  This function is a shortcut assuming exactly three wire planes in the
     * TPC, in which case the output plane is chosen as the one that is neither of the
     * input planes.
     */
    double ThirdPlaneSlope(PlaneID const& pid1,
                           double slope1,
                           PlaneID const& pid2,
                           double slope2) const;

    //@{
    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param plane1 index of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param plane2 index of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @param tpcid TPC where the two planes belong
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a slope as projected in two planes, returns the slope as projected in the
     * third plane.
     */
    double ThirdPlaneSlope(PlaneID::PlaneID_t plane1,
                           double slope1,
                           PlaneID::PlaneID_t plane2,
                           double slope2,
                           TPCID const& tpcid) const
    {
      return ThirdPlaneSlope(PlaneID(tpcid, plane1), slope1, PlaneID(tpcid, plane2), slope2);
    }
    //@}

    /**
     * @brief Returns dT/dW on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first dT/dW
     * @param dTdW1 dT/dW as seen on the first plane
     * @param pid2 ID of the plane of the second dT/dW
     * @param dTdW2 dT/dW  as seen on the second plane
     * @param output_plane ID of the plane on which to calculate the slope
     * @return dT/dW on the third plane, or -999. if dT/dW would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a dT/dW as projected in two planes, returns the dT/dW as projected in the
     * third plane.  The dT/dW are defined in time ticks/wide number units.
     */
    double ThirdPlane_dTdW(PlaneID const& pid1,
                           double slope1,
                           PlaneID const& pid2,
                           double slope2,
                           PlaneID const& output_plane) const;

    /**
     * @brief Returns dT/dW on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first dT/dW
     * @param dTdW1 dT/dW as seen on the first plane
     * @param pid2 ID of the plane of the second dT/dW
     * @param dTdW2 dT/dW  as seen on the second plane
     * @return dT/dW on the third plane, or -999. if dT/dW would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a dT/dW as projected in two planes, returns the dT/dW as projected in the
     * third plane.  This function is a shortcut assuming exactly three wire planes in the
     * TPC, in which case the output plane is chosen as the one that is neither of the
     * input planes.
     */
    double ThirdPlane_dTdW(PlaneID const& pid1,
                           double slope1,
                           PlaneID const& pid2,
                           double slope2) const;

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param angle1 angle or the wires on the first plane
     * @param slope1 slope as observed on the first plane
     * @param angle2 angle or the wires on the second plane
     * @param slope2 slope as observed on the second plane
     * @param angle_target angle or the wires on the target plane
     * @return the slope as measure on the third plane, or 999 if infinity
     *
     * This function will return a small slope if both input slopes are small.
     */
    static double ComputeThirdPlaneSlope(double angle1,
                                         double slope1,
                                         double angle2,
                                         double slope2,
                                         double angle_target);

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param angle1 angle or the wires on the first plane
     * @param pitch1 wire pitch on the first plane
     * @param dTdW1 slope in dt/dw units as observed on the first plane
     * @param angle2 angle or the wires on the second plane
     * @param pitch2 wire pitch on the second plane
     * @param dTdW2 slope in dt/dw units as observed on the second plane
     * @param angle_target angle or the wires on the target plane
     * @param pitch_target wire pitch on the target plane
     * @return dt/dw slope as measured on the third plane, or 999 if infinity
     *
     * The input slope must be specified in dt/dw non-homogeneous coordinates.
     *
     * This function will return a small slope if both input slopes are small.
     */
    static double ComputeThirdPlane_dTdW(double angle1,
                                         double pitch1,
                                         double dTdW1,
                                         double angle2,
                                         double pitch2,
                                         double dTdW2,
                                         double angle_target,
                                         double pitch_target);

    /// @} Wire geometry queries

    /**
     * @name Optical detector geometry access and information
     * @anchor GeometryCoreOpDetGeometry
     * @see @ref GeometryCoreOpDetChannel "optical detector channel information"
     *
     * There are a number of ways to identify an optical detector or channel:
     *
     * * geometric:
     *     * cryostat (e.g. `geo::CryostatID`) and relative optical detector
     *       number within it
     *     * unique optical detector number
     * * readout:
     *     * optical detector channel
     *     * "hardware" channel
     *
     * And they all should be better documented!
     */
    /// @{

    //
    // group features
    //

    /// Number of OpDets in the whole detector
    unsigned int NOpDets() const;

    //
    // access
    //
    /**
     * @brief Returns the `geo::OpDetGeo` object for the given detector number.
     * @param OpDet optical detector unique number
     * @see GeometryCoreOpDetGeometry "optical detector identification"
     */
    OpDetGeo const& OpDetGeoFromOpDet(unsigned int OpDet) const;

    //@{
    /**
     * @brief Find the nearest OpChannel to some point
     * @param xyz point to be queried, in world coordinates
     * @return the nearest OpChannel to the point,
     *         or `std::numeric_limits<unsigned int>::max()` if invalid point
     *
     * @deprecated This method does not tell in which cryostat the detector is;
     *             use `geo::CryostatGeo::GetClosestOpDet()` instead
     *             (find the cryostat with `PositionToCryostatPtr()`).
     */
    unsigned int GetClosestOpDet(Point_t const& point) const;
    //@}

    //
    // object description
    //

    /**
     * @brief Returns gdml string which gives sensitive opdet name
     * @param c ID of the cryostat the detector is in
     *
     * This name is defined in the geometry (GDML) description.
     *
     * @todo Change to use CryostatID
     */
    std::string OpDetGeoName(CryostatID const& cid = cryostat_zero) const;

    /// @} Optical detector access and information

    /// @name Auxiliary detectors access and information
    /// @{

    /// @todo use a AutDetID_t instead of unsigned int?

    //
    // group features
    //

    /**
     * @brief Returns the number of auxiliary detectors
     *
     * This method returns the total number of scintillator paddles (Auxiliary Detectors
     * aka AuxDet) outside of the cryostat
     *
     * @todo Change return type to size_t
     */
    unsigned int NAuxDets() const { return AuxDets().size(); }

    /**
     * @brief Returns the number of sensitive components of auxiliary detector
     * @param aid ID of the auxiliary detector
     * @return number of sensitive components in the auxiliary detector aid
     * @thrws cet::exception (category "Geometry") if aid does not exist
     */
    unsigned int NAuxDetSensitive(size_t const& aid) const;

    //
    // access
    //

    /**
     * @brief Returns the specified auxiliary detector
     * @param ad the auxiliary detector index
     * @return a constant reference to the specified auxiliary detector
     *
     * @todo what happens if it does not exist?
     * @todo remove the default parameter?
     */
    AuxDetGeo const& AuxDet(unsigned int const ad = 0) const;

    /// @} Auxiliary detectors access and information

    /// @name TPC readout channels and views
    /// @{

    //
    // group features
    //

    //
    /**
     * @brief Returns a list of possible views in the detector.
     * @return the set of views
     */
    std::set<View_t> const& Views() const { return allViews; }

    //
    // geometry queries
    //

    /// @} TPC readout channels

    /**
     * @name Optical readout channels
     * @anchor GeometryCoreOpDetChannel
     * @see @ref GeometryCoreOpDetGeometry "optical detector geometry information"
     */
    /// @todo add explanation of the different IDs
    /// @{

    //
    // group features
    //

    /// Get unique opdet number from cryo and internal count
    unsigned int OpDetFromCryo(unsigned int o, unsigned int c) const;

    /// @} Optical readout channels

    //
    // unsorted methods
    //

    /// @name Geometry initialization
    /// @{

    /**
     * @brief Loads the geometry information from the specified files
     * @param gdmlfile path to file to be used for Geant4 simulation
     * @param rootfile path to file for internal geometry representation
     * @param builder algorithm to be used for the interpretation of geometry
     * @param bForceReload reload even if there is already a valid geometry
     *
     * Both paths must directly resolve to an available file, as no search is performed
     * for them.
     *
     * The gdmlfile parameter does not have to necessarily be in GDML format, as long as
     * it's something supported by Geant4. This file is not used by the geometry, but its
     * path is provided on request by the simulation modules (see LArSoft `LArG4` module).
     * The rootfile also does not need to be a ROOT file, but just anything that
     * TGeoManager::Import() supports. This file is parsed immediately and the internal
     * geometry representation is built out of it.
     *
     * @note After calling this method, the detector geometry information can be
     * considered complete, but the geometry service provider is not fully initialized
     * yet, since it's still necessary to provide or update the channel mapping.
     */
    void LoadGeometryFile(std::string gdmlfile,
                          std::string rootfile,
                          GeometryBuilder& builder,
                          bool bForceReload = false);

    /**
     * @brief Loads the geometry information from the specified files
     * @param gdmlfile path to file to be used for Geant4 simulation
     * @param rootfile path to file for internal geometry representation
     * @param bForceReload reload even if there is already a valid geometry
     *
     * This legacy version of `LoadGeometryFile()` uses a standard `geo::GeometryBuilder`
     * implementation.  Do not rely on it if you can avoid it.
     */
    void LoadGeometryFile(std::string gdmlfile, std::string rootfile, bool bForceReload = false);

    CryostatList_t const& Cryostats() const noexcept { return fCryostats; }
    AuxDetList_t const& AuxDets() const noexcept { return fAuxDets; }

  private:
    void SortGeometry();

    std::unique_ptr<GeoObjectSorter const> fSorter;

    CryostatList_t fCryostats{};
    AuxDetList_t fAuxDets{};

    double fSurfaceY;          ///< The point where air meets earth for this detector.
    std::string fDetectorName; ///< Name of the detector.
    std::string fGDMLfile;     ///< path to geometry file used for Geant4 simulation
    std::string fROOTfile;     ///< path to geometry file for geometry in GeometryCore
    double fMinWireZDist;      ///< Minimum distance in Z from a point in which
                               ///< to look for the closest wire
    double fPositionWiggle;    ///< accounting for rounding errors when testing positions

    /// Configuration for the geometry builder
    /// (needed since builder is created after construction).
    fhicl::ParameterSet fBuilderParameters;

    // cached values
    std::set<View_t> allViews; ///< All views in the detector.

    std::vector<TGeoNode const*> FindDetectorEnclosure(
      std::string const& name = "volDetEnclosure") const;

    bool FindFirstVolume(std::string const& name, std::vector<const TGeoNode*>& path) const;

    /// Parses ROOT geometry nodes and builds LArSoft geometry representation.
    void BuildGeometry(GeometryBuilder& builder);

    /// Wire ID check for WireIDsIntersect methods
    bool WireIDIntersectionCheck(const WireID& wid1, const WireID& wid2) const;

    /// Deletes the detector geometry structures
    void ClearGeometry();

  }; // class GeometryCore

  /// @}
  // END Geometry group --------------------------------------------------------

} // namespace geo

//******************************************************************************
//***  template implementation
//***
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// template member function specializations
namespace geo {} // namespace geo

inline geo::GeometryCore::Segment<geo::Point_t> geo::GeometryCore::WireEndPoints(
  WireID const& wireid) const
{
  WireGeo const& wire = Wire(wireid);
  return {wire.GetStart(), wire.GetEnd()};
}

//------------------------------------------------------------------------------
template <typename Stream>
void geo::GeometryCore::Print(Stream&& out, std::string indent /* = "  " */) const
{

  out << "Detector " << DetectorName() << " has " << Ncryostats() << " cryostats and " << NAuxDets()
      << " auxiliary detectors:";

  auto const& detEnclosureBox = DetectorEnclosureBox();
  out << "\n"
      << indent << "Detector enclosure: " << detEnclosureBox.Min() << " -- "
      << detEnclosureBox.Max() << " cm => ( " << detEnclosureBox.SizeX() << " x "
      << detEnclosureBox.SizeY() << " x " << detEnclosureBox.SizeZ() << " ) cm^3";

  for (auto const& cryostat : Iterate<CryostatGeo>()) {
    out << "\n" << indent;
    cryostat.PrintCryostatInfo(std::forward<Stream>(out), indent + "  ", cryostat.MaxVerbosity);

    const unsigned int nTPCs = cryostat.NTPC();
    for (unsigned int t = 0; t < nTPCs; ++t) {
      const TPCGeo& tpc = cryostat.TPC(t);

      out << "\n" << indent << "  ";
      tpc.PrintTPCInfo(std::forward<Stream>(out), indent + "    ", tpc.MaxVerbosity);

      const unsigned int nPlanes = tpc.Nplanes();
      for (unsigned int p = 0; p < nPlanes; ++p) {
        const PlaneGeo& plane = tpc.Plane(p);
        const unsigned int nWires = plane.Nwires();

        out << "\n" << indent << "    ";
        plane.PrintPlaneInfo(std::forward<Stream>(out), indent + "      ", plane.MaxVerbosity);

        for (unsigned int w = 0; w < nWires; ++w) {
          const WireGeo& wire = plane.Wire(w);
          WireID wireID(plane.ID(), w);

          // the wire should be aligned on z axis, half on each side of 0,
          // in its local frame
          out << "\n" << indent << "      " << wireID << " ";
          wire.PrintWireInfo(std::forward<Stream>(out), indent + "      ", wire.MaxVerbosity);
        } // for wire
      }   // for plane
    }     // for TPC

    unsigned int nOpDets = cryostat.NOpDet();
    for (unsigned int iOpDet = 0; iOpDet < nOpDets; ++iOpDet) {
      OpDetGeo const& opDet = cryostat.OpDet(iOpDet);
      out << "\n" << indent << "  [OpDet #" << iOpDet << "] ";
      opDet.PrintOpDetInfo(std::forward<Stream>(out), indent + "  ", opDet.MaxVerbosity);
    } // for
  }   // for cryostat

  unsigned int const nAuxDets = NAuxDets();
  for (unsigned int iDet = 0; iDet < nAuxDets; ++iDet) {
    AuxDetGeo const& auxDet = AuxDet(iDet);

    out << "\n" << indent << "[#" << iDet << "] ";
    auxDet.PrintAuxDetInfo(std::forward<Stream>(out), indent + "  ", auxDet.MaxVerbosity);

    unsigned int const nSensitive = auxDet.NSensitiveVolume();
    switch (nSensitive) {
    case 0: break;
    case 1: {
      AuxDetSensitiveGeo const& auxDetS = auxDet.SensitiveVolume(0U);
      out << "\n" << indent << "  ";
      auxDetS.PrintAuxDetInfo(std::forward<Stream>(out), indent + "    ", auxDetS.MaxVerbosity);
      break;
    }
    default:
      for (unsigned int iSens = 0; iSens < nSensitive; ++iSens) {
        out << "\n" << indent << "[#" << iSens << "] ";
        AuxDetSensitiveGeo const& auxDetS = auxDet.SensitiveVolume(iSens);
        auxDetS.PrintAuxDetInfo(std::forward<Stream>(out), indent + "    ", auxDetS.MaxVerbosity);
      } // for
      break;
    } // if sensitive detectors

  } // for auxiliary detector

  out << '\n';

} // geo::GeometryCore::Print()

#endif // LARCOREALG_GEOMETRY_GEOMETRYCORE_H
