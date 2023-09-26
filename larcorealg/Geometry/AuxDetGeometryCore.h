/**
 * @file   AuxDetGeometryCore.h
 * @brief  Access the description of auxiliary detector geometry
 * @author brebel@fnal.gov
 * @see    AuxDetGeometryCore.cxx
 * @ingroup Geometry
 */
#ifndef GEO_AUXDETGEOMETRYCORE_H
#define GEO_AUXDETGEOMETRYCORE_H

// LArSoft libraries
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetGeoObjectSorter.h"
#include "larcorealg/Geometry/AuxDetReadoutGeom.h"
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// Framework and infrastructure libraries
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <cstddef> // size_t
#include <cstdint> // uint32_t
#include <memory>  // std::shared_ptr<>
#include <string>
#include <vector>

namespace geo {

  /** **************************************************************************
   * @brief Description of geometry of one set of auxiliary detectors
   * @ingroup Geometry
   *
   * @note All lengths are specified in centimetres
   *
   *
   * How to correctly instantiate a GeometryCore object
   * ---------------------------------------------------
   *
   * Instantiation is a multi-step procedure:

   * 1. construct a GeometryCore object (the "service provider"), with the full
   *    configuration; at this step, configuration is just stored

   * 2. load a geometry with GeometryCore::LoadGeometryFile(); this loads the detector
   *    geometry information

   * 3. prepare a channel map algorithm object (might use for example
   *    GeometryCore::DetectorName() or the detector geometry from the newly created
   *    object, but any use of channel mapping related functions is forbidden and it would
   *    yield undefined behaviour (expected to be catastrophic)

   * 4. acquire the channel mapping algorithm with GeometryCore::ApplyChannelMap(); at
   *    this point, the WireReadoutGeom object is asked to initialize itself and to
   *    perform whatever modifications to the geometry provider is needed.
   *
   * Step 3 (creation of the channel mapping algorithm object) can be performed
   * at any time before step 4, provided that no GeometryCore instance is needed
   * for it.
   *
   *
   * Configuration parameters
   * -------------------------
   *
   * - *Name* (string; mandatory): string identifying the detector; it can be different
   *   from the base name of the file used to initialize the geometry; standard names are
   *   recommended by each experiment.  This name can be used, for example, to select
   *   which channel mapping algorithm to use.
   * - *SurfaceY* (real; mandatory): depth of the detector, in centimetrs; see SurfaceY()
   *   for details
   * - *MinWireZDist* (real; default: 3)
   * - *PositionEpsilon* (real; default: 0.01%) set the default tolerance
   *   (see DefaultWiggle())
   */
  class AuxDetGeometryCore {
  public:
    /**
     * @brief Initialize geometry from a given configuration
     * @param pset configuration parameters
     *
     * This constructor does not load any geometry description.
     * The next step is to do exactly that, by GeometryCore::LoadGeometryFile().
     */
    AuxDetGeometryCore(fhicl::ParameterSet const& pset,
                       std::unique_ptr<AuxDetGeoObjectSorter> sorter = nullptr,
                       std::unique_ptr<AuxDetInitializer> initializer = nullptr);

    // You shall not copy or move or assign me!
    AuxDetGeometryCore(AuxDetGeometryCore const&) = delete;
    AuxDetGeometryCore(AuxDetGeometryCore&&) = delete;
    AuxDetGeometryCore& operator=(AuxDetGeometryCore const&) = delete;
    AuxDetGeometryCore& operator=(AuxDetGeometryCore&&) = delete;

    /**
     * @brief Returns the full directory path to the GDML file source
     * @return the full directory path to the GDML file source
     *
     * This is the full path of the source of the detector geometry handed to the detector
     * simulation (GEANT).
     */
    std::string const& GDMLFile() const { return fGDMLfile; }

    /// Returns a string with the name of the detector, as configured
    std::string const& DetectorName() const { return fDetectorName; }

    //
    // object description and information
    //

    /// @todo use a AutDetID_t instead of unsigned int?

    //
    // group features
    //

    /**
     * @brief Returns the number of auxiliary detectors
     *
     * This method returns the total number of scintillator paddles
     * (Auxiliary Detectors aka AuxDet) outside of the cryostat
     */
    std::size_t NAuxDets() const { return fAuxDets.size(); }

    /**
     * @brief Returns the number of sensitive components of auxiliary detector
     * @param aid ID of the auxiliary detector
     * @return number of sensitive components in the auxiliary detector aid
     * @throws cet::exception (category "AuxDetGeometry") if aid does not exist
     */
    std::size_t NAuxDetSensitive(size_t aid) const;

    //
    // access
    //

    /// Returns the full list of pointer to the auxiliary detectors
    std::vector<AuxDetGeo> const& AuxDetGeoVec() const { return fAuxDets; }

    /**
     * @brief Returns the specified auxiliary detector
     * @param ad the auxiliary detector index
     * @return a constant reference to the specified auxiliary detector
     *
     * @todo what happens if it does not exist?
     * @todo remove the default parameter?
     */
    AuxDetGeo const& AuxDet(std::size_t const ad = 0) const;

    Point_t AuxDetChannelToPosition(std::string const& auxDetName, uint32_t channel) const;

    AuxDetSensitiveGeo const& ChannelToAuxDetSensitive(
      std::string const& auxDetName,
      uint32_t channel) const; // return the AuxDetSensitiveGeo for the given

    /**
     * @brief Returns the index of the auxiliary detector at specified location.
     * @param point location to be tested
     * @param tolerance tolerance (cm) for matches. Default 0
     * @return the index of the detector, or `std::numeric_limits<unsigned int>::max()` if
     *        no detector is there
     *
     * @bug Actually, an exception is thrown.
     */
    std::size_t FindAuxDetAtPosition(Point_t const& point, double tolerance = 0) const;

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param point location to be tested
     * @param ad _(output)_ the auxiliary detector index
     * @param tolerance tolerance (cm) for matches. Default 0.
     * @return constant reference to AuxDetGeo object of the auxiliary detector
     *
     * @todo what happens if it does not exist?
     */
    [[nodiscard]] AuxDetGeo const& PositionToAuxDet(Point_t const& point,
                                                    double tolerance = 0.) const;

    /**
     * @brief Fills the indices of the sensitive auxiliary detector at location
     * @param point location to be tested
     * @param adg _(output)_ auxiliary detector index
     * @param sv _(output)_ sensitive volume index
     * @param tolerance tolerance (cm) for matches. Default 0.
     */
    void FindAuxDetSensitiveAtPosition(Point_t const& point,
                                       std::size_t& adg,
                                       std::size_t& sv,
                                       double tolerance = 0) const;

    /// @name Geometry initialization
    /// @{

    /// Returns whether we have a channel map
    bool hasAuxDetChannelMap() const { return bool(fReadoutGeom); }

    /// @}

  private:
    /**
     * @brief Loads the geometry information from the specified files
     * @see ApplyChannelMap()
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

    /**
     * @brief Initializes the geometry to work with this channel map
     * @see LoadGeometryFile()
     *
     * The specified channel mapping is used with this geometry.  The algorithm object is
     * asked and allowed to make the necessary modifications to the geometry description.
     * These modifications typically involve some resorting of the objects.
     *
     * The ownership of the algorithm object is shared, usually with a calling framework:
     * we maintain it alive as long as we need it (and no other code can delete it), and
     * we delete it only if no other code is sharing the ownership.
     *
     * This method needs to be called after LoadGeometryFile() to complete the geometry
     * initialization.
     */

    void ApplyChannelMap();
    std::vector<AuxDetGeo> fAuxDets;

    std::unique_ptr<AuxDetGeoObjectSorter> fSorter;
    std::unique_ptr<AuxDetInitializer> fInitializer;
    std::string fGDMLfile; ///< path to geometry file used for Geant4 simulation
    std::string fDetectorName;
    fhicl::ParameterSet fBuilderParameters; ///< Configuration of geometry builder.
    std::unique_ptr<AuxDetReadoutGeom const>
      fReadoutGeom; ///< Object containing the channel to wire mapping
    bool fThrowIfAbsent;
  }; // class GeometryCore

} // namespace geo

#endif // GEO_AUXDETGEOMETRYCORE_H
