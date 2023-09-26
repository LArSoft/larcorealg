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
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/Iterable.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/details/GeometryIterationPolicy.h"
#include "larcorealg/Geometry/details/ToGeometryElement.h"
#include "larcorealg/Geometry/details/ZeroIDs.h"
#include "larcorealg/Geometry/fwd.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"       // geo::vect namespace
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// Framework and infrastructure libraries
#include "fhiclcpp/fwd.h"

// External libraries
#include "range/v3/view.hpp"

// C/C++ standard libraries
#include <cstddef> // size_t
#include <memory>  // std::shared_ptr<>
#include <set>
#include <string>
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

  class GeometryCore : Iterable<details::GeometryIterationPolicy, details::ToGeometryElement> {
    using Iteration = Iterable<details::GeometryIterationPolicy, details::ToGeometryElement>;

  public:
    /// Type of list of cryostats
    using CryostatList_t = std::vector<CryostatGeo>;
    /// Type of list of auxiliary detectors
    using AuxDetList_t = std::vector<AuxDetGeo>;

    /**
     * @brief Initialize geometry from a given configuration
     * @param pset configuration parameters
     *
     * This constructor does not load any geometry description.  The next step is to do
     * exactly that, by GeometryCore::LoadGeometryFile().
     */
    GeometryCore(fhicl::ParameterSet const& pset,
                 std::unique_ptr<GeometryBuilder> builder,
                 std::unique_ptr<GeoObjectSorter> sorter);
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
    unsigned int Ncryostats() const { return fCryostats.size(); }
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
    CryostatGeo const& Cryostat(CryostatID const& cryoid = details::cryostat_zero) const;
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
      return HasCryostat(cryoid) ? &fCryostats[cryoid.Cryostat] : nullptr;
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
    unsigned int NTPC(CryostatID const& cryoid = details::cryostat_zero) const
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
    TPCGeo const& TPC(TPCID const& tpcid = details::tpc_zero) const
    {
      return Cryostat(tpcid).TPC(tpcid);
    }
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
    std::string OpDetGeoName(CryostatID const& cid = details::cryostat_zero) const;

    /// @} Optical detector access and information

    /// @name Auxiliary detectors access and information
    /// @{

    /// @todo use a AutDetID_t instead of unsigned int?

    //
    // group features
    //

    // //
    // // access
    // //

    // /**
    //  * @brief Returns the specified auxiliary detector
    //  * @param ad the auxiliary detector index
    //  * @return a constant reference to the specified auxiliary detector
    //  *
    //  * @todo what happens if it does not exist?
    //  * @todo remove the default parameter?
    //  */
    // AuxDetGeo const& AuxDet(unsigned int const ad = 0) const;

    /// @} Auxiliary detectors access and information

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

    CryostatList_t const& Cryostats() const noexcept { return fCryostats; }

  private:
    /**
     * @brief Loads the geometry information from the specified files
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
    void LoadGeometryFile();

    void SortGeometry();

    CryostatList_t fCryostats{};

    std::unique_ptr<GeometryBuilder> fBuilder;
    std::unique_ptr<GeoObjectSorter> fSorter;

    TGeoManager* fManager{nullptr};
    std::string fGDMLfile;     ///< path to geometry file used for Geant4 simulation
    std::string fDetectorName; ///< Name of the detector.
    double fSurfaceY;          ///< The point where air meets earth for this detector.
    double fPositionWiggle;    ///< accounting for rounding errors when testing positions

    std::vector<GeoNodePathEntry> FindDetectorEnclosure(
      std::string const& name = "volDetEnclosure") const;

    bool FindFirstVolume(std::string const& name, std::vector<GeoNodePathEntry>& path) const;

    /// Parses ROOT geometry nodes and builds LArSoft geometry representation.
    void BuildGeometry();

  }; // class GeometryCore

  /// @}
  // END Geometry group --------------------------------------------------------

} // namespace geo

//------------------------------------------------------------------------------
template <typename Stream>
void geo::GeometryCore::Print(Stream&& out, std::string indent /* = "  " */) const
{

  out << "Detector " << DetectorName() << " has " << Ncryostats() << " cryostats:";

  auto const& detEnclosureBox = DetectorEnclosureBox();
  out << "\n"
      << indent << "Detector enclosure: " << detEnclosureBox.Min() << " -- "
      << detEnclosureBox.Max() << " cm => ( " << detEnclosureBox.SizeX() << " x "
      << detEnclosureBox.SizeY() << " x " << detEnclosureBox.SizeZ() << " ) cm^3";

  for (auto const& cryostat : Iterate<CryostatGeo>()) {
    out << "\n" << indent;
    cryostat.PrintCryostatInfo(std::forward<Stream>(out), indent + "  ", cryostat.MaxVerbosity);

    unsigned const int nTPCs = cryostat.NTPC();
    for (unsigned int t = 0; t < nTPCs; ++t) {
      TPCGeo const& tpc = cryostat.TPC(t);

      out << "\n" << indent << "  ";
      tpc.PrintTPCInfo(std::forward<Stream>(out), indent + "    ", tpc.MaxVerbosity);
    } // for TPC

    unsigned int nOpDets = cryostat.NOpDet();
    for (unsigned int iOpDet = 0; iOpDet < nOpDets; ++iOpDet) {
      OpDetGeo const& opDet = cryostat.OpDet(iOpDet);
      out << "\n" << indent << "  [OpDet #" << iOpDet << "] ";
      opDet.PrintOpDetInfo(std::forward<Stream>(out), indent + "  ", opDet.MaxVerbosity);
    } // for
  }   // for cryostat

  out << '\n';

} // geo::GeometryCore::Print()

#endif // LARCOREALG_GEOMETRY_GEOMETRYCORE_H
