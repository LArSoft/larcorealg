#ifndef GEO_CHANNELMAPALG_H
#define GEO_CHANNELMAPALG_H
// vim: sw=2 expandtab :
////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/ChannelMapAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// LArSoft  libraries
#include "larcorealg/Geometry/Iterable.h"
#include "larcorealg/Geometry/details/ReadoutIterationPolicy.h"
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

// C/C++ standard libraries
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace geo {

  /**
   * @brief Interface for a class providing readout channel mapping to geometry
   * @ingroup Geometry
   *
   * @note A number of methods react specifically when provided with invalid IDs
   * as arguments. An invalid ID is an ID with the `isValid` flag unset, or, in
   * case of channel IDs, an ID with value `raw::InvalidChannelID`.
   * An ID that does not present this feature is by definition "valid"; this
   * does not imply that the represented entity (channel, geometry entity or
   * readout group) actually exists. *The behaviour of the methods to valid,
   * non-existent IDs should be considered undefined*, and it is recommended
   * that the existence of the entity is checked beforehand.
   * Unless the documentation explicitly defines a behaviour, an undefined
   * behaviour should be assumed; nevertheless, the documentation of some of the
   * methods still reminds of this.
   */
  class ChannelMapAlg : Iterable<details::ReadoutIterationPolicy> {
    using Iteration = Iterable<details::ReadoutIterationPolicy>;

  public:
    explicit ChannelMapAlg(GeometryCore const* geom);
    virtual ~ChannelMapAlg();

    //--------------------------------------------------------------------------
    // Iteration facilities
    using Iteration::begin;
    using Iteration::end;
    using Iteration::Iterate;

    //--------------------------------------------------------------------------
    /// @{
    /// @name TPC channel mapping

    /// Returns the total number of channels present (not necessarily contiguous)
    virtual unsigned int Nchannels() const = 0;

    /// @brief Returns the number of channels in the specified ROP
    /// @return number of channels in the specified ROP, 0 if non-existent
    virtual unsigned int Nchannels(readout::ROPID const& ropid) const = 0;

    /// @brief Returns an std::vector<ChannelID_t> in all TPCs in a TPCSet
    std::vector<raw::ChannelID_t> ChannelsInTPCs() const;

    /// @brief Returns whether the specified channel is valid
    /// This default implementation assumes all channels up to Nchannels() valid.
    virtual bool HasChannel(raw::ChannelID_t channel) const
    {
      return raw::isValidChannelID(channel) ? channel < Nchannels() : false;
    }

    /// Returns a list of TPC wires connected to the specified readout channel ID
    /// @throws cet::exception (category: "Geometry") if non-existent channel
    virtual std::vector<WireID> ChannelToWire(raw::ChannelID_t channel) const = 0;

    /**
     * @brief Returns the view (wire orientation) on the specified TPC channel
     * @param channel TPC channel ID
     * @return the type of signal on the specified channel, or geo::kUnknown
     *
     * The view of the readout plane `channel` belongs to is returned, as in
     * `View(readout::ROPID const&) const`.
     */
    View_t View(raw::ChannelID_t const channel) const;

    /**
     * @brief Returns the view of the channels in the specified readout plane
     * @param ropid readout plane ID
     * @return the type of signal on the specified ROP
     *
     * Returns the view (wire orientation) on the channels of specified readout
     * plane.
     * If ropid is an invalid ID, geo::kUnknown is returned.
     * If ropid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent readout plane, the result is undefined.
     * Use HasROP() to check if the readout plane actually exists.
     */
    View_t View(readout::ROPID const& ropid) const;

    /**
     * @brief Returns the type of signal on the channels of specified TPC plane
     * @param plane TPC plane ID
     * @return the type of signal on the specified plane, or geo::kMysteryType
     *
     * Assumes that all the channels on the plane have the same signal type.
     *
     * @todo verify that kMysteryType is returned on invalid plane
     */
    SigType_t SignalType(PlaneID const& pid) const;

    /**
     * @brief Return the signal type of the specified channel
     * @param channel ID of the channel
     * @return signal type of the channel, or geo::kMysteryType if not known
     *
     * On any type of error (e.g., invalid or unknown channel ID),
     * geo::kMysteryType is returned.
     */
    SigType_t SignalType(raw::ChannelID_t const channel) const;

    /**
     * @brief Return the signal type on the specified readout plane
     * @param ropid ID of the readout plane
     * @return signal type on the plane, or geo::kMysteryType if not known
     *
     * If the readout plane ID is marked invalid, geo::kMysteryType is returned.
     * If the readout plane is not marked invalid, but it does not match an
     * existing readout plane, the result is undefined.
     *
     * The default implementation uses readout plane to channel mapping.
     * Other implementation may decide to do the opposite.
     */
    SigType_t SignalType(readout::ROPID const& ropid) const;

    /**
     * @brief Returns the ID of the channel nearest to the specified position
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param planeid ID of the wire plane the channel must belong to
     * @return the ID of the channel, or `raw::InvalidChannelID` if invalid wire
     * @bug on invalid wire, a `geo::InvalidWireError` exception is thrown
     *
     */
    raw::ChannelID_t NearestChannel(Point_t const& worldLoc, PlaneID const& planeid) const;

    /**
     * @brief Returns an intersection point of two channels
     * @param c1 one channel ID
     * @param c2 the other channel ID
     * @param y (output) y coordinate of the intersection
     * @param z (output) z coordinate of the intersection
     * @return whether a intersection point was found
     *
     * @todo what happens for channels from different TPCs?
     * @todo what happens for channels with multiple intersection points?
     *
     * @deprecated This is clearly not APA-aware
     */
    std::optional<WireIDIntersection> ChannelsIntersect(raw::ChannelID_t c1,
                                                        raw::ChannelID_t c2) const;

    /// Returns a list of the plane IDs in the whole detector
    virtual std::set<PlaneID> const& PlaneIDs() const = 0;

    /**
     * @brief Returns the channel ID a wire is connected to
     * @param wireID ID of the wire
     * @return the ID of the channel
     * @see PlaneWireToChannel(geo::WireID const&)
     *
     * Behaviour on an invalid or not present wires is undefined.
     */
    virtual raw::ChannelID_t PlaneWireToChannel(WireID const& wireID) const = 0;

    /// @}

    //--------------------------------------------------------------------------
    /// @{
    /// @name Optical detector channel mapping

    /**
     * @brief Returns the number of optical channels contained in some detectors
     * @param NOpDets number of optical detectors
     * @return optical channels contained in NOpDets detectors
     *
     * This function returns how many channels can be expected to be present
     * in a detector with NOpDets optical detectors. This is an upper limit, as
     * not all channels have necessarily to be present.
     *
     * For example: if a detector has four channels per optical detector, the
     * returned value will be four times the argument NOpDets. If there is a
     * single channel on each optical detector, the return value will be the
     * value NOpDets (this is in fact the fallback implementation). If each
     * optical detector can have anywhere between 2 and 12 channels, the
     * returned value is 12 times NOpDets, and it will be an overestimation of
     * the number of channels.
     */
    virtual unsigned int NOpChannels(unsigned int NOpDets) const;

    /// Number of electronics channels for all the optical detectors
    unsigned int NOpChannels() const;

    /// Largest optical channel number
    unsigned int MaxOpChannel() const;

    /// Is this a valid OpChannel number?
    bool IsValidOpChannel(int opChannel) const;

    /**
     * @brief Returns the number of optical channels contained in some detectors
     * @param NOpDets number of optical detectors
     * @return optical channels contained in NOpDets detectors
     *
     * This function returns the first optical channel ID larger than the last
     * channel ID in a detector with NOpDets optical detectors (with the same
     * logic as `NOpChannels()`).
     * For example, in a detector with 32 channels with contiguous
     * IDs starting at 0, this function would return 32. If the channels started
     * with ID 1, this function would instead return 33 and if there were a 16
     * channel gap, so that valid channels are from 0 to 15 and from 32 to 47,
     * this function would return 48.
     */
    virtual unsigned int MaxOpChannel(unsigned int NOpDets) const;

    /**
     * @brief Returns the number of channels in the specified optical detectors
     * @param opDet ID of the chosen optical detector
     * @return optical channels contained in optical detector with ID opDet
     *
     * This function returns how many channels are actually present in the
     * optical detector with the specified ID.
     *
     * For example: if a detector has four channels per optical detector, the
     * returned value will be four, regardless opDet, and . If there is a
     * single channel on each optical detector, the return value will be 1,
     * again ignoring opDet (this is in fact the fallback implementation). If
     * each optical detector can have anywhere between 2 and 12 channels, the
     * returned value will be 2, 12, etc., that is the exact number of channels
     * in opDet.
     *
     * Although implementations are encouraged to return 0 on invalid optical
     * detectors, the actual return value in that case is undefined.
     */
    virtual unsigned int NOpHardwareChannels(unsigned int opDet) const;

    /**
     * @brief Returns whether the ID identifies a valid optical detector channel
     * @param opChannel channel number
     * @param NOpDets number of optical detectors in the detector
     * @return whether opChannel would be a valid channel
     *
     * The specification of the number of optical channels reflects the logic
     * described in `NOpChannel()`.
     */
    virtual bool IsValidOpChannel(unsigned int opChannel, unsigned int NOpDets) const;

    /**
     * @brief Returns the channel ID of the specified hardware channel
     * @param detNum optical detector ID
     * @param hwchannel hardware channel within the specified optical detector
     * @return ID of the channel identified by detector and hardware channel IDs
     *
     * If the input IDs identify a non-existing channel, the result is
     * undefined.
     */
    virtual unsigned int OpChannel(unsigned int detNum, unsigned int hwchannel = 0) const;

    /**
     * @brief Returns the optical detector the specified optical channel belongs
     * @param opChannel the optical detector channel being queried
     * @return the optical detector the specified optical channel belongs to
     *
     * If the specified optical channel is invalid, behaviour is undefined.
     */
    virtual unsigned int OpDetFromOpChannel(unsigned int opChannel) const;

    /**
     * @brief Returns the optical detector the specified optical channel belongs
     * @param opChannel the optical detector channel being queried
     * @return the optical detector geometry object the specified optical channel belongs to
     *
     * If the specified optical channel is invalid, behaviour is undefined.
     */
    OpDetGeo const& OpDetGeoFromOpChannel(unsigned int opChannel) const;

    /**
     * @brief Returns the hardware channel number of specified optical channel
     * @param opChannel the optical detector channel being queried
     * @return the optical detector the specified optical channel belongs to
     *
     * If the specified optical channel is invalid, behaviour is undefined.
     */
    virtual unsigned int HardwareChannelFromOpChannel(unsigned int opChannel) const;

    /// @}

    //--------------------------------------------------------------------------
    /// @{
    /// @name Mapping of position to wires

    /**
     * @brief Returns the index of the wire nearest to the specified position
     * @param YPos y coordinate on the wire plane
     * @param ZPos z coordinate on the wire plane
     * @param planeID ID of the plane
     * @return an index interpolation between the two nearest wires
     * @see NearestWireID()
     *
     * Respect to NearestWireID(), this method returns a real number,
     * representing a continuous coordinate in the wire axis, with the round
     * values corresponding to the actual wires.
     *
     * The plane is required to be valid and exist in the detector. Otherwise,
     * the behaviour is undefined.
     */
    virtual double WireCoordinate(double YPos, double ZPos, PlaneID const& planeID) const = 0;

    /**
     * @brief Returns the ID of the wire nearest to the specified position
     * @param worldPos position to be tested
     * @param planeID plane containing the wire
     * @return the ID of the wire closest to worldPos in the specified plane
     * @throw InvalidWireIDError the ID found is not present in the detector
     * @see WireCoordinate(double, double, geo::PlaneID const&)
     *
     * The plane is required to be valid and exist in the detector. Otherwise,
     * the behaviour is undefined.
     * An exception is thrown if the wire that would be the closest is actually
     * not present; but no check is performed whether the specified position is
     * outside the wire plane: wires are extrapolated to be infinitely long.
     * In other words, the result can be trusted only as long as the position
     * is within the specified wire plane.
     */
    virtual WireID NearestWireID(Point_t const& worldPos, PlaneID const& planeID) const = 0;

    /// @}

    //--------------------------------------------------------------------------
    /// @{
    /// @name Auxiliary detectors

    /**
     * @brief Returns the auxiliary detector closest to the specified point
     * @param point coordinates of the position to be investigated (x, y, z)
     * @param auxDets list of the sought auxiliary detectors
     * @param tolerance tolerance for comparison. Default 0.
     * @return index of auxiliary detector within auxDets
     */
    virtual size_t NearestAuxDet(Point_t const& point,
                                 std::vector<AuxDetGeo> const& auxDets,
                                 double tolerance = 0) const;

    /**
     * @brief Returns sensitive auxiliary detector closest to specified point
     * @param point coordinates of the position to be investigated (x, y, z)
     * @param auxDets list of the auxiliary detectors
     * @param tolerance tolerance for comparison. Default 0.
     * @return index of sought sensitive auxiliary detector within auxDets
     */
    virtual size_t NearestSensitiveAuxDet(Point_t const& point,
                                          std::vector<AuxDetGeo> const& auxDets,
                                          double tolerance = 0) const;

    /**
     * @brief Returns the index of the detector containing the specified channel
     * @param auxDets list of the auxiliary detectors
     * @param detName name of the auxiliary detector being investigated
     * @param channel number of the channel within that auxiliary detector
     * @return index of the sought auxiliary detector within auxDets
     *
     * @bug This function is somehow broken in that it ignores the `auxDets`
     *      in the arguments and instead relies on a cache that is never filled
     *      by this class (derived classes can fill it though).
     */
    virtual size_t ChannelToAuxDet(std::vector<AuxDetGeo> const& auxDets,
                                   std::string const& detName,
                                   uint32_t const& channel) const;

    /**
     * @brief Returns the index of the sensitive detector containing the channel
     * @param auxDets list of the sensitive auxiliary detectors
     * @param detName name of the auxiliary detector being investigated
     * @param channel number of the channel within that auxiliary detector
     * @return index of the sought sensitive auxiliary detector within auxDets
     */
    virtual std::pair<size_t, size_t> ChannelToSensitiveAuxDet(
      std::vector<AuxDetGeo> const& auxDets,
      std::string const& detName,
      uint32_t const& channel) const;

    /**
     * @brief Returns the index of the auxiliary detector at specified location.
     * @param point location to be tested
     * @param tolerance tolerance (cm) for matches. Default 0
     * @return the index of the detector, or
     *        `std::numeric_limits<unsigned int>::max()` if no detector is there
     *
     * @bug Actually, an exception is thrown.
     */
    unsigned int FindAuxDetAtPosition(Point_t const& point, double tolerance = 0) const;

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param point location to be tested
     * @param ad _(output)_ the auxiliary detector index
     * @param tolerance tolerance (cm) for matches. Default 0.
     * @return constant reference to AuxDetGeo object of the auxiliary detector
     *
     * @todo what happens if it does not exist?
     */
    AuxDetGeo const& PositionToAuxDet(Point_t const& point,
                                      unsigned int& ad,
                                      double tolerance = 0) const;

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

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param point location to be tested
     * @param ad _(output)_ the auxiliary detector index
     * @param sv _(output)_ the auxiliary detector sensitive volume index
     * @param tolerance tolerance (cm) for matches. Default 0.
     * @return reference to AuxDetSensitiveGeo object of the auxiliary detector
     *
     * @todo what happens if it does not exist?
     */
    AuxDetSensitiveGeo const& PositionToAuxDetSensitive(Point_t const& point,
                                                        size_t& ad,
                                                        size_t& sv,
                                                        double tolerance = 0) const;

    AuxDetGeo const& ChannelToAuxDet(std::string const& auxDetName, uint32_t channel) const;

    AuxDetSensitiveGeo const& ChannelToAuxDetSensitive(std::string const& auxDetName,
                                                       uint32_t channel) const;

    /// @}

    //--------------------------------------------------------------------------
    /// @name TPC set mapping
    /// @{
    /**
     * @brief Returns the total number of TPC sets in the specified cryostat
     * @param cryoid cryostat ID
     * @return number of TPC sets in the cryostat, or 0 if no cryostat found
     */
    virtual unsigned int NTPCsets(readout::CryostatID const& cryoid) const = 0;

    /// Returns the largest number of TPC sets any cryostat in the detector has
    virtual unsigned int MaxTPCsets() const = 0;

    /// Returns whether we have the specified TPC set
    /// @return whether the TPC set is valid and exists
    virtual bool HasTPCset(readout::TPCsetID const& tpcsetid) const = 0;

    readout::TPCsetID FindTPCsetAtPosition(Point_t const& worldLoc) const;

    /// Returns the ID of the TPC set tpcid belongs to
    virtual readout::TPCsetID TPCtoTPCset(TPCID const& tpcid) const = 0;

    /**
     * @brief Returns a list of ID of TPCs belonging to the specified TPC set
     * @param tpcsetid ID of the TPC set to convert into TPC IDs
     * @return the list of TPCs, empty if TPC set is invalid
     *
     * Note that the check is performed on the validity of the TPC set ID, that
     * does not necessarily imply that the TPC set specified by the ID actually
     * exists. Check the existence of the TPC set first (HasTPCset()).
     * Behaviour on valid, non-existent TPC set IDs is undefined.
     */
    virtual std::vector<TPCID> TPCsetToTPCs(readout::TPCsetID const& tpcsetid) const = 0;

    /// Returns the ID of the first TPC belonging to the specified TPC set
    virtual TPCID FirstTPCinTPCset(readout::TPCsetID const& tpcsetid) const = 0;

    /// @} TPC set mapping

    //--------------------------------------------------------------------------
    /// @name Readout plane mapping
    /// @{
    /**
     * @brief Returns the total number of ROP in the specified TPC set
     * @param tpcsetid TPC set ID
     * @return number of readout planes in the TPC set, or 0 if no TPC set found
     *
     * Note that this methods explicitly check the existence of the TPC set.
     */
    virtual unsigned int NROPs(readout::TPCsetID const& tpcsetid) const = 0;

    /// Returns the largest number of ROPs a TPC set in the detector has
    virtual unsigned int MaxROPs() const = 0;

    /// Returns whether we have the specified ROP
    /// @return whether the readout plane is valid and exists
    virtual bool HasROP(readout::ROPID const& ropid) const = 0;

    /// Returns the ID of the ROP planeid belongs to
    virtual readout::ROPID WirePlaneToROP(PlaneID const& planeid) const = 0;

    /// Returns a list of ID of planes belonging to the specified ROP
    virtual std::vector<PlaneID> ROPtoWirePlanes(readout::ROPID const& ropid) const = 0;

    /// Returns the ID of the first plane belonging to the specified ROP
    virtual PlaneID FirstWirePlaneInROP(readout::ROPID const& ropid) const = 0;

    /**
     * @brief Returns a list of ID of TPCs the specified ROP spans
     * @param ropid ID of the readout plane
     * @return the list of TPC IDs, empty if readout plane ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour on non-existing readout planes is undefined.
     */
    virtual std::vector<TPCID> ROPtoTPCs(readout::ROPID const& ropid) const = 0;

    /**
     * @brief Returns the ID of the ROP the channel belongs to
     * @return the ID of the ROP the channel belongs to (invalid if channel is)
     * @see HasChannel()
     *
     * The channel must exist, or be the invalid channel value.
     * With a channel that is not present in the mapping and that is not the
     * invalid channel (`raw::InvalidChannelID`), the result is undefined.
     */
    virtual readout::ROPID ChannelToROP(raw::ChannelID_t channel) const = 0;

    /**
     * @brief Returns the ID of the first channel in the specified readout plane
     * @param ropid ID of the readout plane
     * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour for non-existing readout planes is undefined.
     */
    virtual raw::ChannelID_t FirstChannelInROP(readout::ROPID const& ropid) const = 0;

    /// @}

    //--------------------------------------------------------------------------
    /// @{
    /// @name Testing (not in the interface)

    /// Retrieve the private fFirstChannelInNextPlane vector for testing.
    std::vector<std::vector<std::vector<raw::ChannelID_t>>> const& FirstChannelInNextPlane() const
    {
      return fFirstChannelInThisPlane;
    }

    /// Retrieve the private fFirstChannelInThisPlane vector for testing.
    std::vector<std::vector<std::vector<raw::ChannelID_t>>> const& FirstChannelInThisPlane() const
    {
      return fFirstChannelInNextPlane;
    }

    /// @}

    //--------------------------------------------------------------------------

  protected:
    /**
     * @brief Return the signal type of the specified channel
     * @param channel ID of the channel
     * @return signal type of the channel, or geo::kMysteryType if not known
     *
     * On any type of error (e.g., invalid or unknown channel ID),
     * geo::kMysteryType is returned.
     */
    virtual SigType_t SignalTypeForChannelImpl(raw::ChannelID_t const channel) const = 0;

    /**
     * @brief Return the signal type on the specified readout plane
     * @param ropid ID of the readout plane
     * @return signal type on the plane, or geo::kMysteryType if not known
     *
     * If the readout plane ID is marked invalid, geo::kMysteryType is returned.
     * If the readout plane is not marked invalid, but it does not match an
     * existing readout plane, the result is undefined.
     *
     * The default implementation uses readout plane to channel mapping.
     * Other implementation may decide to do the opposite.
     */
    virtual SigType_t SignalTypeForROPIDImpl(readout::ROPID const& ropid) const;

    /// Data type for per-TPC information
    template <typename T>
    using TPCInfoMap_t = std::vector<std::vector<T>>;

    /// Data type for per-plane information
    template <typename T>
    using PlaneInfoMap_t = TPCInfoMap_t<std::vector<T>>;

    /**
    * @name Internal structure data access
    *
    * These functions allow access to the XxxInfoMap_t types based on geometry
    * element IDs.
    * They are strictly internal.
    */
    /// @{

    /// Returns the specified element of the TPC map
    template <typename T>
    T const& AccessElement(TPCInfoMap_t<T> const& map, TPCID const& id) const
    {
      return map[id.Cryostat][id.TPC];
    }

    /// Returns the number of elements in the specified cryostat of the TPC map
    template <typename T>
    size_t AccessElementSize(TPCInfoMap_t<T> const& map, CryostatID const& id) const
    {
      return map[id.Cryostat].size();
    }

    //@{
    /// Returns whether the ID specifies a valid entry
    template <typename T>
    bool isValidElement(TPCInfoMap_t<T> const& map, CryostatID const& id) const
    {
      return id.Cryostat < map.size();
    }
    template <typename T>
    bool isValidElement(TPCInfoMap_t<T> const& map, TPCID const& id) const
    {
      return isValidElement(map, id.asCryostatID()) && (id.TPC < map[id.Cryostat].size());
    }
    //@}

    /// Returns the specified element of the plane map
    template <typename T>
    T const& AccessElement(PlaneInfoMap_t<T> const& map, PlaneID const& id) const
    {
      return map[id.Cryostat][id.TPC][id.Plane];
    }

    /// Returns the number of elements in the specified TPC of the plane map
    template <typename T>
    size_t AccessElementSize(PlaneInfoMap_t<T> const& map, TPCID const& id) const
    {
      return map[id.Cryostat][id.TPC].size();
    }

    //@{
    /// Returns whether the ID specifies a valid entry
    template <typename T>
    bool isValidElement(PlaneInfoMap_t<T> const& map, CryostatID const& id) const
    {
      return id.Cryostat < map.size();
    }
    template <typename T>
    bool isValidElement(PlaneInfoMap_t<T> const& map, TPCID const& id) const
    {
      return isValidElement(map, id.asCryostatID()) && (id.TPC < map[id.Cryostat].size());
    }
    template <typename T>
    bool isValidElement(PlaneInfoMap_t<T> const& map, PlaneID const& id) const
    {
      return isValidElement(map, id.asTPCID()) && (id.Plane < AccessSize(map, id.asTPCID()));
    }
    //@}

    /// Returns a pointer to the specified element, or nullptr if invalid
    template <typename T>
    T const* GetElementPtr(PlaneInfoMap_t<T> const& map, PlaneID const& id) const
    {
      if (id.Cryostat >= map.size()) return nullptr;
      auto const& cryo_map = map[id.Cryostat];
      if (id.TPC >= cryo_map.size()) return nullptr;
      auto const& TPC_map = cryo_map[id.TPC];
      if (id.Plane >= TPC_map.size()) return nullptr;
      auto const& plane_map = TPC_map[id.Plane];
      return &plane_map;
    } // GetElementPtr()

    ///@} Internal structure data access

    // These 3D vectors are used in initializing the Channel map.
    // Only a 1D vector is really needed so far, but these are more general.
    PlaneInfoMap_t<raw::ChannelID_t> fFirstChannelInThisPlane;
    PlaneInfoMap_t<raw::ChannelID_t> fFirstChannelInNextPlane;

    std::map<std::string, size_t>
      fADNameToGeo; ///< map the names of the dets to the AuxDetGeo objects
    std::map<size_t, std::vector<size_t>>
      fADChannelToSensitiveGeo; ///< map the AuxDetGeo index to a vector of
                                ///< indices corresponding to the AuxDetSensitiveGeo index
  private:
    GeometryCore const* fGeom;
  };

}
#endif // GEO_CHANNELMAPALG_H
