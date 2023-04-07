////////////////////////////////////////////////////////////////////////
/// \file  WireReadoutStandardGeom.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARCOREALG_GEOMETRY_CHANNELSTANDARDMAPALG_H
#define LARCOREALG_GEOMETRY_CHANNELSTANDARDMAPALG_H

#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcorealg/Geometry/WireReadoutSorter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h" // readout::TPCsetID, ...

#include "fhiclcpp/fwd.h"

#include <set>
#include <vector>

namespace geo {

  class WireReadoutStandardGeom : public WireReadoutGeom {
  public:
    WireReadoutStandardGeom(fhicl::ParameterSet const& pset,
                            GeometryCore const* geom,
                            std::unique_ptr<WireReadoutSorter> sorter);

    std::vector<WireID> ChannelToWire(raw::ChannelID_t channel) const override;
    unsigned int Nchannels() const override;

    /// @brief Returns the number of channels in the specified ROP
    /// @return number of channels in the specified ROP, 0 if non-existent
    unsigned int Nchannels(readout::ROPID const& ropid) const override;

    //@{
    double WireCoordinate(double YPos, double ZPos, PlaneID const& planeID) const override;
    //@}

    //@{
    WireID NearestWireID(Point_t const& worldPos, PlaneID const& planeID) const override;
    //@}

    //@{
    raw::ChannelID_t PlaneWireToChannel(WireID const& wireID) const override;
    //@}

    std::set<PlaneID> const& PlaneIDs() const override;

    //
    // TPC set interface
    //
    /// @name TPC set mapping
    /// @{
    /**
     * @brief Returns the total number of TPC sets in the specified cryostat
     * @param cryoid cryostat ID
     * @return number of TPC sets in the cryostat, or 0 if no cryostat found
     *
     * In this mapping, TPCs have independent readout and there is one TPC in each TPC set
     * and one TPC set for each TPC.
     */
    unsigned int NTPCsets(readout::CryostatID const& cryoid) const override;

    /// Returns the largest number of TPC sets any cryostat in the detector has
    unsigned int MaxTPCsets() const override;

    /// Returns whether we have the specified TPC set
    /// @return whether the TPC set is valid and exists
    bool HasTPCset(readout::TPCsetID const& tpcsetid) const override;

    /**
     * @brief Returns the ID of the TPC set the specified TPC belongs to
     * @param tpcid ID of the TPC
     * @return the ID of the corresponding TPC set, or invalid ID when tpcid is
     *
     * In this mapping, TPC sets and TPCs are mapped one-to-one.  The returned value
     * mirrors the TPC ID in the readout space.  If the TPC ID is not valid, an invalid
     * TPC set ID is returned.  Note that this check is performed on the validity of the
     * TPC ID, that does not necessarily imply that the TPC specified by the ID actually
     * exists.
     */
    readout::TPCsetID TPCtoTPCset(TPCID const& tpcid) const override;

    /**
     * @brief Returns a list of ID of TPCs belonging to the specified TPC set
     * @param tpcsetid ID of the TPC set to convert into TPC IDs
     * @return the list of TPCs, empty if TPC set is invalid
     *
     * In this mapping, TPC sets and TPCs are mapped one-to-one.  The returned list
     * contains always one entry, unless the specified TPC set ID is invalid, in which
     * case the list is empty.  Note that the check is performed on the validity of the
     * TPC set ID, that does not necessarily imply that the TPC set specified by the ID
     * actually exists. Check the existence of the TPC set first (HasTPCset()).  Behaviour
     * on valid, non-existent TPC set IDs is undefined.
     */
    std::vector<TPCID> TPCsetToTPCs(readout::TPCsetID const& tpcsetid) const override;

    /// Returns the ID of the first TPC belonging to the specified TPC set
    TPCID FirstTPCinTPCset(readout::TPCsetID const& tpcsetid) const override;

    /// @} TPC set mapping

    //
    // Readout plane interface
    //
    /// @name Readout plane mapping
    /// @{

    /**
     * @brief Returns the total number of ROPs in the specified TPC set
     * @param tpcsetid TPC set ID
     * @return number of readout planes in the TPC sets, or 0 if ID is invalid
     *
     * Note that this methods explicitly check the existence of the TPC set.
     *
     * In this mapping, planes have independent readout and there is one wire plane in
     * each readout plane and one readout plane for each wire plane.
     */
    unsigned int NROPs(readout::TPCsetID const& tpcsetid) const override;

    /// Returns the largest number of ROPs a TPC set in the detector has
    unsigned int MaxROPs() const override;

    /// Returns whether we have the specified ROP
    /// @return whether the readout plane is valid and exists
    bool HasROP(readout::ROPID const& ropid) const override;

    /**
     * @brief Returns the ID of the ROP planeid belongs to, or invalid if none
     * @param planeid ID of the plane
     * @return the ID of the corresponding ROP, or invalid ID when planeid is
     *
     * In this mapping, readout planes and wire planes are mapped one-to-one.  The
     * returned value mirrors the plane ID in the readout space.  If the plane ID is not
     * valid, an invalid readout plane ID is returned.  Note that this check is performed
     * on the validity of the plane ID, that does not necessarily imply that the plane
     * specified by the ID actually exists.
     */
    readout::ROPID WirePlaneToROP(PlaneID const& planeid) const override;

    /**
     * @brief Returns a list of ID of wire planes belonging to the specified ROP
     * @param ropid ID of the readout plane to convert into wire planes
     * @return the list of wire plane IDs, empty if readout plane ID is invalid
     *
     * In this mapping, readout planes and wire planes are mapped one-to-one.  The
     * returned list contains always one entry, unless the specified readout plane ID is
     * invalid, in which case the list is empty.  Note that this check is performed on the
     * validity of the readout plane ID, that does not necessarily imply that the readout
     * plane specified by the ID actually exists.
     */
    std::vector<PlaneID> ROPtoWirePlanes(readout::ROPID const& ropid) const override;

    /**
     * @brief Returns a list of ID of TPCs the specified ROP spans
     * @param ropid ID of the readout plane
     * @return the list of TPC IDs, empty if readout plane ID is invalid
     *
     * In this mapping, readout planes and wire planes are mapped one-to-one.  The
     * returned list contains always one entry, unless the specified readout plane ID is
     * invalid, in which case the list is empty.  Note that this check is performed on the
     * validity of the readout plane ID, that does not necessarily imply that the readout
     * plane specified by the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour on non-existing readout planes is undefined.
     */
    std::vector<TPCID> ROPtoTPCs(readout::ROPID const& ropid) const override;

    /// Returns the ID of the ROP the channel belongs to (invalid if none)
    readout::ROPID ChannelToROP(raw::ChannelID_t channel) const override;

    /**
     * @brief Returns the ID of the first channel in the specified readout plane
     * @param ropid ID of the readout plane
     * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane ID, that
     * does not necessarily imply that the readout plane specified by the ID actually
     * exists. Check if the ROP exists with HasROP().  The behaviour for non-existing
     * readout planes is undefined.
     */
    raw::ChannelID_t FirstChannelInROP(readout::ROPID const& ropid) const override;

    /// Returns the ID of the first plane belonging to the specified ROP
    PlaneID FirstWirePlaneInROP(readout::ROPID const& ropid) const override;

    /// @} readout plane mapping

  private:
    unsigned int fNcryostat;              ///< number of cryostats in the detector
    unsigned int fNchannels;              ///< number of channels in the detector
    raw::ChannelID_t fTopChannel;         ///< book keeping highest channel #
    std::vector<unsigned int> fNTPC;      ///< number of TPCs in each cryostat
    std::set<View_t> fViews;              ///< vector of the views present in the detector
    std::set<PlaneID> fPlaneIDs;          ///< vector of the PlaneIDs present in the detector
    PlaneInfoMap_t<float> fFirstWireProj; ///< Distance (0,0,0) to first wire
                                          ///< along orth vector per plane per TPC
    PlaneInfoMap_t<float> fOrthVectorsY;  ///< Unit vectors orthogonal to wires in
    PlaneInfoMap_t<float> fOrthVectorsZ;  ///< each plane - stored as 2 components
                                          ///< to avoid having to invoke any bulky
                                          ///< TObjects / CLHEP vectors etc
    PlaneInfoMap_t<float> fWireCounts;    ///< Number of wires in each plane - for
                                          ///< range checking after calculation
    TPCInfoMap_t<unsigned int> fNPlanes;  ///< Number of planes in each TPC - for
                                          ///< range checking after calculation
    PlaneInfoMap_t<unsigned int> fPlaneBaselines; ///< The number of wires in all the
                                                  ///< tpcs and planes up to this one
                                                  ///< in the heirachy
    PlaneInfoMap_t<unsigned int> fWiresPerPlane;  ///< The number of wires in this plane
                                                  ///< in the heirachy

    SigType_t SignalTypeForChannelImpl(raw::ChannelID_t const channel) const override;

    /// Retrieved the wire cound for the specified plane ID
    unsigned int WireCount(PlaneID const& id) const { return AccessElement(fWireCounts, id); }

    /// Returns the largest number of TPCs in a single cryostat
    unsigned int MaxTPCs() const;
  };

}
#endif // LARCOREALG_GEOMETRY_CHANNELSTANDARDMAPALG_H
