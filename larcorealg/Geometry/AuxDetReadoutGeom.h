////////////////////////////////////////////////////////////////////////
/// @file  AuxDetReadoutGeom.h
/// @brief Interface to auxiliary-detector geometry for wire readouts.
/// @ingroup Geometry
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_AUXDETREADOUTGEOM_H
#define LARCOREALG_GEOMETRY_AUXDETREADOUTGEOM_H

#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// C/C++ standard libraries
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace geo {

  using chanAndSV = std::pair<std::uint32_t, std::size_t>;

  struct AuxDetReadoutInitializers {
    std::map<std::size_t, std::string> ADGeoToName; ///< map the AuxDetGeo index to the name
    std::map<std::string, std::size_t> NameToADGeo; ///< map the names to the AuxDetGeo index
    std::map<std::size_t, std::vector<chanAndSV>>
      ADGeoToChannelAndSV; ///< map the AuxDetGeo index to a vector of pairs corresponding
                           ///< to the channel and AuxDetSensitiveGeo index
  };

  class AuxDetInitializer {
  public:
    AuxDetReadoutInitializers init(std::vector<AuxDetGeo> const& ads) const;
    virtual ~AuxDetInitializer();

  private:
    virtual AuxDetReadoutInitializers initialize(std::vector<AuxDetGeo> const& ads) const = 0;
  };

  /// @ingroup Geometry
  class AuxDetReadoutGeom {
  public:
    explicit AuxDetReadoutGeom(AuxDetReadoutInitializers initializers);

    // method returns the entry in the sorted AuxDetGeo vector so that the Geometry in
    // turn can return that object
    std::size_t NearestAuxDet(Point_t const& point,
                              std::vector<AuxDetGeo> const& auxDets,
                              double tolerance = 0,
                              bool throwIfAbsent = true) const;
    std::size_t NearestSensitiveAuxDet(Point_t const& point,
                                       std::vector<AuxDetGeo> const& auxDets,
                                       double tolerance = 0,
                                       bool throwIfAbsent = true) const;
    std::pair<std::size_t, std::size_t> ChannelToSensitiveAuxDet(std::string const& detName,
                                                                 std::uint32_t channel) const;

    Point_t AuxDetChannelToPosition(std::uint32_t channel,
                                    std::string const& auxDetName,
                                    std::vector<AuxDetGeo> const& auxDets) const;

  private:
    std::size_t DetNameToAuxDet(std::string const& detName) const;

    std::map<std::size_t, std::string> fADGeoToName; ///< map the AuxDetGeo index to the name
    std::map<std::string, std::size_t> fNameToADGeo; ///< map the names to the AuxDetGeo index
    std::map<std::size_t, std::vector<chanAndSV>>
      fADGeoToChannelAndSV; ///< map the AuxDetGeo index to a vector of pairs
                            ///< corresponding to the channel and AuxDetSensitiveGeo index
  };
}
#endif // LARCOREALG_GEOMETRY_AUXDETREADOUTGEOM_H
