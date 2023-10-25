/**
 * @file   AuxDetGeometryCore.h
 * @brief  Access the description of auxiliary detector geometry
 * @see    AuxDetGeometryCore.cxx
 * @ingroup Geometry
 */

#ifndef LARCOREALG_GEOMETRY_AUXDETGEOMETRYCORE_H
#define LARCOREALG_GEOMETRY_AUXDETGEOMETRYCORE_H

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
   * @brief Description of physical geometry of one set of auxiliary detectors
   * @ingroup Geometry
   *
   * @note All lengths are specified in centimeters
   *
   *
   * AuxDetGeometryCore construction
   * -------------------------------
   *
   * The constructor of the AuxDetGeometryCore performs two steps:
   *
   * 1. It initializes a AuxDetGeometryCoreobject (the "service provider") with the configuration
   *    parameters, a user-provided geometry builder, and a user-provided geometry sorter.
   *
   * 2. It loads the geometry with AuxDetGeometryCore::LoadGeometryFile(), which:
   *    a. imports the GDML file into ROOT's TGeo* system
   *    b. builds all LArSoft geometry constructs using ROOT's geometry system according
   *       to the user-provided builder
   *    c. sorts all LArSoft geometry constructs according to the provided sorter.
   *
   *
   * Configuration parameters
   * ------------------------
   *
   * - *GDML* (string; mandatory): string identifying GDML file that represents the
        detector goemetry.
   * - *RelativePath* (string; default: (empty)): string identifying path of GDML file
        relative to any paths on the `FW_SEARCH_PATH` environment variable.
   * - *Name* (string; default: stem of GDML file): string identifying the detector; it
   *    can be different from the base name of the file used to initialize the geometry;
   *    standard names are recommended by each experiment.
   * - *Builder* (parameter set; default: {}): parameters used in constructing
        GeometryBuilderStandard
   * - *ThrowIfAbsent* (boolean; default: true) if true throws when the requested AuxDet
        element is not available; otherwise logs a message and returns a nonsense value
        when possible.
   */

  class AuxDetGeometryCore {
  public:
    /**
     * @brief Initialize geometry from a given configuration
     * @param pset configuration parameters
     * @param sorter sorter used for sorting AuxDet elements
     * @param initializer initializer used to construct AuxDetReadoutGeom object
     * @see   AuxDetReadoutGeom.cxx
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
     * This method returns the total number of scintillator paddles (Auxiliary Detectors
     * aka AuxDet) outside of the cryostat
     */
    std::size_t NAuxDets() const { return fAuxDets.size(); }

    /**
     * @brief Returns the number of sensitive components of auxiliary detector
     * @param ad ID of the auxiliary detector
     * @return number of sensitive components in the auxiliary detector aid
     * @throws cet::exception (category "AuxDetGeometryCore") if ad is out of range
     */
    std::size_t NAuxDetSensitive(size_t ad) const;

    //
    // access
    //

    /// Returns the full list of pointer to the auxiliary detectors
    std::vector<AuxDetGeo> const& AuxDetGeoVec() const { return fAuxDets; }

    /**
     * @brief Returns the specified auxiliary detector
     * @param ad the auxiliary detector index
     * @return a constant reference to the specified auxiliary detector
     * @throws cet::exception (category "AuxDetGeometryCore") if ad is out of range
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
     * @param tolerance tolerance (cm) for matches. Default 0.
     * @return constant reference to AuxDetGeo object of the auxiliary detector
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
    void LoadGeometryFile();
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

#endif // LARCOREALG_GEOMETRY_AUXDETGEOMETRYCORE_H
