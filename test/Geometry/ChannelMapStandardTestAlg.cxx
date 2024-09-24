/**
 * @file   ChannelMapStandardTestAlg.cxx
 * @brief  Tests the standard channel mapping algorithm.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 26th, 2015
 *
 * The methods require a Boost test enviroment.
 */

// LArSoft libraries
#include "ChannelMapStandardTestAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

// framework
#include "cetlib_except/exception.h"

// Boost libraries
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <string>
#include <vector>

namespace {

  /// Checks that there is a 1-to-1 match between ID data members
  void CheckMatchingTPClevelIDs(readout::TPCsetID const& tpcsetID, geo::TPCID const& tpcID)
  {
    BOOST_TEST(tpcID.isValid == tpcsetID.isValid);
    BOOST_TEST(tpcID.Cryostat == tpcsetID.Cryostat);
    BOOST_TEST(tpcID.TPC == tpcsetID.TPCset);
  } // CheckMatchingTPClevelIDs()

  /// Checks that there is a 1-to-1 match between ID data members
  void CheckMatchingPlaneLevelIDs(readout::ROPID const& ropID, geo::PlaneID const& planeID)
  {
    CheckMatchingTPClevelIDs(ropID, planeID);
    BOOST_TEST(planeID.Plane == ropID.ROP);
  } // CheckMatchingPlaneLevelIDs()

} // local namespace

//-----------------------------------------------------------------------------
unsigned int geo::ChannelMapStandardTestAlg::Run()
{
  // All the tests

  TPCsetMappingTest();
  ROPMappingTest();
  ChannelMappingTest();

  return 0;
} // ChannelMapStandardTestAlg::Run()

//-----------------------------------------------------------------------------
void geo::ChannelMapStandardTestAlg::TPCsetMappingTest() const
{

  /*
   * WireReadoutStandardGeom interface to be tested (via GeometryCore):
   *
   *     unsigned int NTPCsets(readout::CryostatID const& cryoid) const
   *
   *     unsigned int MaxTPCsets() const
   *
   *     bool HasTPCset(readout::TPCsetID const& tpcsetid) const
   *
   *     readout::TPCsetID TPCtoTPCset(geo::TPCID const& tpcid) const
   *
   *     std::vector<geo::TPCID> TPCsetToTPCs
   *       (readout::TPCsetID const& tpcsetid) const
   *
   *     geo::TPCID FirstTPCinTPCset(readout::TPCsetID const& tpcsetid) const
   *       can't do (not exposed!)
   *
   *     unsigned int NROPs(readout::TPCsetID const& tpcsetid) const
   *
   *     unsigned int MaxROPs() const
   *
   *     bool HasROP(readout::ROPID const& ropid) const
   *       (invalid input is checked in ROPMappingTest())
   *
   *     std::vector<geo::TPCID> ROPtoTPCs(readout::ROPID const& ropid) const
   *
   */

  // check for invalid input
  BOOST_TEST(wireReadoutGeom->NTPCsets({}) == 0U);
  BOOST_TEST(!wireReadoutGeom->HasTPCset({}));
  BOOST_TEST(!wireReadoutGeom->TPCtoTPCset({}).isValid);
  BOOST_TEST(wireReadoutGeom->TPCsetToTPCs({}).empty());
  BOOST_TEST(wireReadoutGeom->NROPs({}) == 0U);
  BOOST_TEST(!wireReadoutGeom->HasROP({}));
  BOOST_TEST(wireReadoutGeom->ROPtoTPCs({}).empty());

  //
  // detector-wide checks
  //
  // check that the maximum size of TPC and TPC sets in the detector match
  BOOST_TEST(wireReadoutGeom->MaxTPCsets() == geom->MaxTPCs());

  // check that the maximum size of planes and ROPs in the detector match
  BOOST_TEST(wireReadoutGeom->MaxROPs() == wireReadoutGeom->MaxPlanes());

  // check that we have no TPC set in cryostats after the last one
  readout::CryostatID const NonexistingCryostatID(geom->Ncryostats());
  BOOST_TEST(wireReadoutGeom->NTPCsets(NonexistingCryostatID) == 0U);

  //
  // cryostat-wide checks
  //
  for (geo::CryostatID const& cryostatID : geom->Iterate<CryostatID>()) {
    BOOST_TEST_CHECKPOINT("cryostat: " << std::string(cryostatID));

    readout::CryostatID const ROcryostatID = static_cast<readout::CryostatID const&>(cryostatID);

    // check that the number of TPC and TPC sets in each cryostat match
    unsigned int const NTPCsets = wireReadoutGeom->NTPCsets(ROcryostatID);
    BOOST_TEST(NTPCsets == geom->NTPC(cryostatID));

    // check that we have no TPC set after the last one
    readout::TPCsetID const NonexistingTPCsetID(ROcryostatID,
                                                (readout::TPCsetID::TPCsetID_t)NTPCsets);
    BOOST_TEST(!wireReadoutGeom->HasTPCset(NonexistingTPCsetID));
    BOOST_TEST(wireReadoutGeom->NROPs(NonexistingTPCsetID) == 0U);
  } // for cryostats

  //
  // TPC-wide checks
  //
  for (geo::TPCID const& tpcID : geom->Iterate<TPCID>()) {
    BOOST_TEST_CHECKPOINT("TPC: " << std::string(tpcID));

    // check that the IDs of this TPC and of the TPC set including it match
    readout::TPCsetID const tpcsetID = wireReadoutGeom->TPCtoTPCset(tpcID);
    CheckMatchingTPClevelIDs(tpcsetID, tpcID);

    // check that this TPC set maps to only one (and the right one) TPC
    std::vector<geo::TPCID> TPCs = wireReadoutGeom->TPCsetToTPCs(tpcsetID);
    BOOST_TEST(TPCs.size() == 1U);
    BOOST_TEST(TPCs.front() == tpcID);

    // check that the number of ROP in the TPC set matches the planes in the TPC
    unsigned int const NROPs = wireReadoutGeom->NROPs(tpcsetID);
    BOOST_TEST(NROPs == wireReadoutGeom->Nplanes(tpcID));

    // check all the ROPs:
    for (unsigned int iROPinTPCset = 0; iROPinTPCset < NROPs; ++iROPinTPCset) {
      readout::ROPID const ropID(tpcsetID, (readout::ROPID::ROPID_t)iROPinTPCset);

      BOOST_TEST_CHECKPOINT("ROP: " << std::string(ropID));

      // do we have it?
      BOOST_TEST(wireReadoutGeom->HasROP(ropID));

      // is it in the right TPC? and only one?
      std::vector<geo::TPCID> TPCs = wireReadoutGeom->ROPtoTPCs(ropID);
      BOOST_TEST(TPCs.size() == 1U);
      BOOST_TEST(TPCs.front() == tpcID);

    } // for channels

  } // for TPCs

} // ChannelMapStandardTestAlg::TPCsetMappingTest()

//-----------------------------------------------------------------------------
void geo::ChannelMapStandardTestAlg::ROPMappingTest() const
{

  /*
   * WireReadoutStandardGeom interface to be tested (via GeometryCore):
   *
   *     unsigned int Nchannels(readout::ROPID const& ropid) const
   *
   *     readout::ROPID WirePlaneToROP(geo::PlaneID const& planeid) const
   *
   *     std::vector<geo::PlaneID> ROPtoWirePlanes
   *       (readout::ROPID const& ropid) const
   *
   *     std::vector<geo::TPCID> ROPtoTPCs(readout::ROPID const& ropid) const
   *
   *     readout::ROPID ChannelToROP(raw::ChannelID_t channel) const
   *
   *     raw::ChannelID_t FirstChannelInROP(readout::ROPID const& ropid) const
   *
   *     geo::PlaneID FirstWirePlaneInROP(readout::ROPID const& ropid) const
   *       can't do (not exposed!)
   *
   */

  // check for invalid input
  BOOST_TEST(wireReadoutGeom->Nchannels({}) == 0U);
  BOOST_TEST(!wireReadoutGeom->WirePlaneToROP({}).isValid);
  BOOST_TEST(wireReadoutGeom->ROPtoWirePlanes({}).empty());
  BOOST_TEST(wireReadoutGeom->ROPtoTPCs({}).empty());
  BOOST_TEST(!wireReadoutGeom->ChannelToROP(raw::InvalidChannelID).isValid);
  BOOST_TEST(!raw::isValidChannelID(wireReadoutGeom->FirstChannelInROP({})));

  //
  // TPC-wide checks
  //
  for (geo::TPCID const& tpcID : geom->Iterate<TPCID>()) {
    BOOST_TEST_CHECKPOINT("TPC: " << std::string(tpcID));

    // build a non-existent ROP ID (but we pretend it valid)
    readout::TPCsetID const tpcsetID = wireReadoutGeom->TPCtoTPCset(tpcID);
    unsigned int const NROPs = wireReadoutGeom->NROPs(tpcsetID);
    readout::ROPID NonexistingROPID(tpcsetID, (readout::ROPID::ROPID_t)NROPs);

    // check that we don't have ROPs beyond the last one
    BOOST_TEST(!wireReadoutGeom->HasROP(NonexistingROPID));
  } // for TPCs

  //
  // plane-wide checks
  //
  for (geo::PlaneID const& planeID : wireReadoutGeom->Iterate<PlaneID>()) {
    BOOST_TEST_MESSAGE("plane: " << std::string(planeID));

    // check that the ROP of this plane matches the plane itself
    readout::ROPID const ropID = wireReadoutGeom->WirePlaneToROP(planeID);
    CheckMatchingPlaneLevelIDs(ropID, planeID);

    // check that there is only one plane in this ROP
    std::vector<geo::PlaneID> const PlanesInROP = wireReadoutGeom->ROPtoWirePlanes(ropID);
    BOOST_TEST(PlanesInROP.size() == 1U);
    BOOST_TEST(PlanesInROP.front() == planeID);

    // check that the number of channels in the ROP matches the number of wires
    unsigned int const NChannels = wireReadoutGeom->Nchannels(ropID);
    BOOST_TEST(NChannels == wireReadoutGeom->Nwires(planeID));

    // check that the TPC is one, and the right one
    std::vector<geo::TPCID> const TPCs = wireReadoutGeom->ROPtoTPCs(ropID);
    BOOST_TEST(TPCs.size() == 1U);
    BOOST_TEST(TPCs.front() == planeID.asTPCID());

    // check that the first channel is valid
    raw::ChannelID_t const FirstChannelID = wireReadoutGeom->FirstChannelInROP(ropID);
    BOOST_TEST(raw::isValidChannelID(FirstChannelID) == ropID.isValid);

    // check all the channels:
    for (unsigned int iChannelInROP = 0; iChannelInROP < NChannels; ++iChannelInROP) {
      raw::ChannelID_t const channelID = FirstChannelID + iChannelInROP;

      BOOST_TEST_CHECKPOINT("channel: " << channelID);

      // is it in the right plane? and only one?
      std::vector<geo::WireID> const ChannelWires = wireReadoutGeom->ChannelToWire(channelID);
      BOOST_TEST(ChannelWires.size() == 1U);
      BOOST_TEST(ChannelWires.front() == planeID);

      // does the channel map back to the right ROP?
      readout::ROPID const ChannelROPID = wireReadoutGeom->ChannelToROP(channelID);
      BOOST_TEST(ChannelROPID == ropID);

    } // for channels

  } // for planes

} // ChannelMapStandardTestAlg::ROPMappingTest()

//-----------------------------------------------------------------------------
void geo::ChannelMapStandardTestAlg::ChannelMappingTest() const
{

  /*
   * WireReadoutStandardGeom interface being tested (via GeometryCore):
   *
   *     unsigned int Nchannels() const
   *       (only called, not checked)
   *
   *     unsigned int HasChannel(raw::ChannelID_t) const
   *
   * The rest is not tested here.
   */

  // check for invalid input
  BOOST_TEST(!wireReadoutGeom->HasChannel(raw::InvalidChannelID));

  //
  // channel-wide checks
  //
  unsigned int const NChannels = wireReadoutGeom->Nchannels();
  for (unsigned int iChannel = 0; iChannel < NChannels; ++iChannel) {
    raw::ChannelID_t channel = (raw::ChannelID_t)iChannel;

    BOOST_TEST_MESSAGE("channel: " << channel);

    BOOST_TEST(wireReadoutGeom->HasChannel(iChannel));

  } // for channels
  BOOST_TEST(!wireReadoutGeom->HasChannel((raw::ChannelID_t)NChannels));

} // ChannelMapStandardTestAlg::ChannelMappingTest()

//-----------------------------------------------------------------------------
