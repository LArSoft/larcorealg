/**
 * @file   GeometryIteratorTestAlg.cxx
 * @brief  Unit test for geometry iterators
 * @date   May 7th, 2015
 * @author petrillo@fnal.gov
 *
 * The methods require a Boost test enviroment.
 */

// LArSoft libraries
#include "GeometryIteratorTestAlg.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "test/Geometry/IteratorTypes.h"

// Boost libraries
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <string>
#include <type_traits>

using namespace geo::details;

namespace {

  /**
   * @brief Checks that the geometry iterator behaves like the corresponding ID
   *        iterator
   *
   * The following methods are tested:
   *
   *     geometry_element_iterator(geometry_element_iterator const&)
   *     bool operator== (iterator const& as) const // result always true
   *     bool operator== (id_iterator_t const& as) const // result always true
   *     bool operator!= (iterator const& as) const [result always false
   *     bool operator!= (id_iterator_t const& as) const // result always false
   *     Element_t const& operator* () const
   *     Element_t const* operator-> () const
   *     iterator& operator++ ()
   *     iterator operator++ (int)
   *     operator bool() const
   *     ElementPtr_t get() const
   *     LocalID_t const& ID() const
   *
   * Checks are implemented by Boost tests.
   */
  template <typename ITER, typename ITERID>
  void CompareIteratorAndIteratorID(ITER iter, ITERID id_iter)
  {
    static_assert(std::is_same<typename ITER::id_iterator_t, ITERID>{},
                  "CompareIteratorAndIteratorID() requires compatible iterator types");

    // ID comparison
    BOOST_TEST(iter.ID() == *id_iter);

    // check copy assignment
    ITER iter_copy(iter);
    ITERID id_iter_copy(id_iter);

    // check comparisons too
    BOOST_TEST(iter == iter_copy);
    BOOST_TEST(iter_copy == iter);
    BOOST_TEST(!(iter != iter_copy));
    BOOST_TEST(!(iter_copy != iter));

    // check increment operator
    BOOST_TEST((iter++).ID() == *(id_iter++));
    BOOST_TEST((++iter_copy).ID() == *(++id_iter_copy));

    BOOST_TEST(iter == iter_copy);

  } // CompareIteratorAndIteratorID()

} // local namespace

//-----------------------------------------------------------------------------
unsigned int geo::GeometryIteratorTestAlg::Run() const
{
  // All the tests

  // - geometry ID iterators
  CryostatIDIteratorsTest();
  TPCIDIteratorsTest();
  PlaneIDIteratorsTest();
  WireIDIteratorsTest();
  TPCsetIDIteratorsTest();
  ROPIDIteratorsTest();

  // - geometry element iterators
  CryostatIteratorsTest();
  TPCIteratorsTest();
  PlaneIteratorsTest();
  WireIteratorsTest();

  return 0;
} // GeometryIteratorTestAlg::Run()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::CryostatIDIteratorsTest() const
{
  GeometryIterationPolicy const policy{geom};
  /*
   * public interface (cryostat_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     cryostat_id_iterator_base();
   *
   *     /// Constructor: points to the specified cryostat
   *     cryostat_id_iterator_base
   *       (geo::GeometryCore const* geom, GEOID const& start_from);
   *
   *     /// Returns true if the two iterators point to the same cryostat
   *     template <typename OTHERID>
   *     bool operator== (cryostat_id_iterator_base<OTHERID> const& as) const;
   *
   *     /// Returns true if the two iterators point to different cryostats
   *     template <typename OTHERID>
   *     bool operator!= (cryostat_id_iterator_base<OTHERID> const& as) const;
   *
   *     /// Returns the ID the iterator points to
   *     LocalID_t const& operator* () const;
   *
   *     /// Returns a pointer to the ID the iterator points to
   *     LocalID_t const* operator-> () const;
   *
   *
   *     /// Returns whether the iterator is pointing to a valid cryostat
   *     operator bool() const;
   *
   *     /// Returns a pointer to cryostat, or nullptr if invalid
   *     ElementPtr_t get() const;
   *
   */

  //
  // default constructed
  //
  {
    geo::cryostat_id_iterator iCryo;
    BOOST_TEST_CHECKPOINT("Default created cryostat ID iterator: " << std::string(*iCryo));

    BOOST_TEST(iCryo->Cryostat == geo::CryostatID::getInvalidID());
  }

  //
  // begin-constructed
  //
  {
    // initialize to the beginning directly; this has probably ID's isValid true
    auto BeginID = geo::CryostatID::first();
    geo::cryostat_id_iterator iCryo(BeginID, policy);
    BOOST_TEST(iCryo->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iCryo == iCryo);

    // construct from explicit begin position
    geo::cryostat_id_iterator iCryoBC{BeginID, policy};
    BOOST_TEST(iCryoBC == iCryo);

    // construct at begin position by geometry
    geo::cryostat_id_iterator iCryoGB = geom->begin<CryostatID>();
    BOOST_TEST(iCryoGB == iCryo);

    // check access to ID
    BOOST_TEST(*iCryo == BeginID);
    BOOST_TEST(iCryo->Cryostat == BeginID.Cryostat);

    // check access to geometry element
    geo::CryostatGeo const* pCryo = geom->CryostatPtr(BeginID);
    geo::cryostat_iterator const iCryoElem{geom, iCryo};
    BOOST_TEST(iCryoElem.get() == pCryo);

    // test copy and postfix increment
    geo::cryostat_id_iterator iCryoI(iCryo++);

    BOOST_TEST(iCryo->Cryostat == geo::CryostatID::CryostatID_t(1));
    BOOST_TEST(iCryoI->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iCryoI != iCryo);

    // test copy and prefix increment
    ++iCryoI;
    BOOST_TEST(iCryoI->Cryostat == geo::CryostatID::CryostatID_t(1));
    BOOST_TEST(iCryoI == iCryo);

    if (geom->Ncryostats() > 1) {
      ++iCryoI;
      BOOST_TEST(iCryoI->Cryostat == geo::CryostatID::CryostatID_t(2));
      BOOST_TEST(iCryoI != iCryo);
    }
  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last TPC
    geo::CryostatID LastID(geom->Ncryostats() - 1); // last cryostat

    geo::cryostat_id_iterator iLastCryo(LastID, policy);
    BOOST_TEST_CHECKPOINT(
      "Position-created iterator to last cryostat ID: " << std::string(*iLastCryo));

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastCryo == LastID);
    BOOST_TEST(iLastCryo->Cryostat == LastID.Cryostat);
    geo::cryostat_iterator iLastCryoElem{geom, iLastCryo};
    BOOST_TEST(iLastCryoElem);
    BOOST_TEST(iLastCryoElem.get() == geom->CryostatPtr(LastID));

    // test increment to past-the-end
    geo::cryostat_id_iterator iEndCryo = iLastCryo;
    ++iEndCryo;
    ++iLastCryoElem;

    // check that the iterator tests false
    BOOST_TEST(!iLastCryoElem);

    BOOST_TEST(iEndCryo->Cryostat == geom->Ncryostats());
    BOOST_TEST(iEndCryo == geom->end<CryostatID>());
    BOOST_TEST(!iLastCryoElem.get());
  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::cryostat_id_iterator iCryo{policy.GetEndID<geo::CryostatID>(), policy};
    BOOST_TEST_CHECKPOINT("End-created cryostat ID iterator: " << std::string(*iCryo));

    BOOST_TEST(iCryo->Cryostat == geom->Ncryostats());

    // check access to geometry element (result of operator* is not defined)
    geo::cryostat_iterator iCryoElem{geom, iCryo};
    BOOST_TEST(!iCryoElem);
    BOOST_TEST(!(iCryoElem.get())); // should get nullptr

    // construct at end position by geometry
    geo::cryostat_id_sentinel iCryoGE = geom->end<CryostatID>();
    BOOST_TEST(iCryoGE == iCryo);

    // initialize to the end directly; this has probably ID's isValid true
    geo::cryostat_id_iterator iCryo2(geo::CryostatID(geom->Ncryostats()), policy);
    BOOST_TEST(iCryo2->Cryostat == geom->Ncryostats());
    BOOST_TEST(iCryo2 == iCryo);
  }

} // GeometryIteratorTestAlg::CryostatIDIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::CryostatIteratorsTest() const
{
  GeometryIterationPolicy const policy{geom};
  /*
   * This test is extensively based on the assumption that the iterators should
   * behave like the corresponding ID iterators, including "corner cases".
   *
   * public interface (geometry_element_iterator<>)
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     geometry_element_iterator()
   *
   *     /// Constructor: points to the same element as the specified ID iterator
   *     geometry_element_iterator(id_iterator_t const& iter)
   *
   *     /// Constructor: points to the same element as the specified ID iterator
   *     geometry_element_iterator(id_iterator_t&& iter)
   *
   *     /// Constructor: points to the specified geometry element
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, GeoID_t const& start_from)
   *
   *     /// Returns true if the two iterators point to the same object
   *     bool operator== (iterator const& as) const
   *
   *     /// Returns true if the two iterators point to the same object
   *     bool operator== (id_iterator_t const& as) const
   *
   *     /// Returns true if the two iterators point to different objects
   *     bool operator!= (iterator const& as) const
   *
   *     /// Returns true if the two iterators point to different objects
   *     bool operator!= (id_iterator_t const& as) const
   *
   *     /// Returns the geometry element the iterator points to
   *     Element_t const& operator* () const
   *
   *     /// Returns a pointer to the element the iterator points to (or nullptr)
   *     Element_t const* operator-> () const
   *
   *     /// Prefix increment: returns this iterator pointing to the next element
   *     iterator& operator++ ()
   *
   *     /// Postfix increment: returns the current iterator, then increments it
   *     iterator operator++ (int)
   *
   *     /// Returns whether the iterator is pointing to a valid geometry element
   *     operator bool() const
   *
   *     /// Returns a pointer to the geometry element, or nullptr if invalid
   *     ElementPtr_t get() const
   *
   *     /// Returns the ID of the pointed geometry element
   *     LocalID_t const& ID() const
   *
   */

  //
  // default constructed
  //
  {
    geo::cryostat_id_iterator iCryoID;
    BOOST_TEST_CHECKPOINT("Default created cryostat iterator: " << std::string(*iCryoID));

    geo::cryostat_iterator iCryo;

    // ID comparison
    BOOST_TEST(iCryo.ID() == *iCryoID);

    // check copy assignment
    geo::cryostat_iterator iCryo_copy(iCryo);
    geo::cryostat_id_iterator iCryoID_copy(iCryoID);

    // check comparisons too
    BOOST_TEST(iCryo == iCryo_copy);
    BOOST_TEST(iCryo_copy == iCryo);
    BOOST_TEST(!(iCryo != iCryo_copy));
    BOOST_TEST(!(iCryo_copy != iCryo));

    BOOST_TEST(iCryo == iCryo_copy);
    BOOST_TEST(iCryoID_copy == iCryoID);
  }

  //
  // begin-constructed
  //
  {
    auto BeginID = geo::CryostatID::first();

    geo::cryostat_id_iterator iCryoID(BeginID, policy);

    BOOST_TEST_CHECKPOINT("Begin-created cryostat iterator (" << std::string(BeginID) << ")");

    // initialize to the beginning directly
    geo::cryostat_iterator iCryoD(geom, iCryoID);
    CompareIteratorAndIteratorID(iCryoD, iCryoID);

    // construct from explicit begin position
    geo::cryostat_iterator iCryoBC{geom, geom->begin<CryostatID>()};
    CompareIteratorAndIteratorID(iCryoBC, iCryoID);

    // construct at begin position by geometry
    geo::cryostat_iterator iCryoGB = geom->begin<CryostatGeo>();
    CompareIteratorAndIteratorID(iCryoGB, iCryoID);
  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last TPC
    geo::CryostatID LastID(geom->Ncryostats() - 1); // last cryostat
    geo::cryostat_id_iterator iLastCryoID(LastID, policy);

    BOOST_TEST_CHECKPOINT("Position-created iterator to last cryostat: " << std::string(LastID));
    geo::cryostat_iterator iLastCryo(geom, iLastCryoID);

    // test increment to past-the-end
    geo::cryostat_id_iterator iEndCryoID = iLastCryoID;
    ++iEndCryoID;

    geo::cryostat_iterator iEndCryo = iLastCryo;
    ++iEndCryo;

    CompareIteratorAndIteratorID(iEndCryo, iEndCryoID);
  }

  //
  // end-constructed
  //
  {
    geo::cryostat_id_sentinel iCryoID = geom->end<CryostatID>();

    // construct from end position
    BOOST_TEST_CHECKPOINT("End-created cryostat iterator");
    geo::cryostat_sentinel iCryo{geom->end<CryostatID>()};
    BOOST_TEST(iCryo == iCryoID);

    // construct at end position by geometry
    geo::cryostat_sentinel iCryoGE = geom->end<CryostatGeo>();
    BOOST_TEST(iCryoGE == iCryoID);
  }

} // GeometryIteratorTestAlg::CryostatIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::TPCIDIteratorsTest() const
{
  GeometryIterationPolicy const policy{geom};
  /*
   * public interface (TPC_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     TPC_id_iterator_base()
   *
   *     /// Constructor: points to begin
   *     TPC_id_iterator_base(geo::GeometryCore const* geom)
   *
   *     /// Constructor: points to the specified cryostat
   *     TPC_id_iterator_base
   *       (geo::GeometryCore const* geom, GEOID const& start_from)
   *
   *     /// Constructor: points to begin
   *     TPC_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t)
   *
   *     /// Constructor: points to end
   *     TPC_id_iterator_base(geo::GeometryCore const* geom, EndPos_t)
   *
   *     /// Returns true if the two iterators point to the same TPC
   *     template <typename OTHERID>
   *     bool operator== (TPC_id_iterator_base<OTHERID> const& as) const
   *
   *     /// Returns true if the two iterators point to different TPCs
   *     template <typename OTHERID>
   *     bool operator!= (TPC_id_iterator_base<OTHERID> const& as) const
   *
   *     /// Returns the TPCID the iterator points to
   *     LocalID_t const& operator* () const
   *
   *     /// Returns the TPCID the iterator points to
   *     LocalID_t const* operator-> () const
   *
   *     /// Prefix increment: returns this iterator pointing to the next TPC
   *     iterator& operator++ ()
   *
   *     /// Postfix increment: returns the current iterator, then increments it
   *     iterator operator++ (int)
   *
   *     /// Returns whether the iterator is pointing to a valid TPC
   *     operator bool() const
   *
   *     /// Returns a pointer to TPC, or nullptr if invalid
   *     ElementPtr_t get() const
   *
   */

  //
  // default constructed
  //
  {
    geo::TPC_id_iterator iTPC;
    BOOST_TEST_CHECKPOINT("Default created TPC ID iterator: " << std::string(*iTPC));

    BOOST_TEST(iTPC->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iTPC->TPC == geo::TPCID::getInvalidID());
  }

  //
  // begin-constructed
  //
  {
    // initialize to the beginning directly; this has probably ID's isValid true
    auto BeginID = geo::TPCID::first();
    geo::TPC_id_iterator iTPC(BeginID, policy);
    BOOST_TEST(iTPC->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iTPC->TPC == geo::TPCID::TPCID_t(0));

    // construct from explicit begin position
    geo::TPC_id_iterator iTPCBC{BeginID, policy};
    BOOST_TEST(iTPCBC == iTPC);

    // construct at begin position by geometry
    geo::TPC_id_iterator iTPCGB = geom->begin<TPCID>();
    BOOST_TEST(iTPCGB == iTPC);

    // check access to ID
    BOOST_TEST(*iTPC == BeginID);
    BOOST_TEST(iTPC->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iTPC->TPC == BeginID.TPC);

    // check access to geometry element
    geo::TPCGeo const* pTPC = geom->TPCPtr(BeginID);
    geo::TPC_iterator const iTPCElem{geom, iTPC};
    BOOST_TEST(iTPCElem);
    BOOST_TEST(iTPCElem.get() == pTPC);

    // test copy and postfix increment
    geo::TPC_id_iterator iTPCI(iTPC++);

    const unsigned int nTPCsInC0 = geom->NTPC(geo::CryostatID(0));
    if (nTPCsInC0 > 1) {
      BOOST_TEST(iTPCI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iTPCI->TPC == geo::TPCID::TPCID_t(0));
      BOOST_TEST(iTPC->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iTPC->TPC == geo::TPCID::TPCID_t(1));
    }
    BOOST_TEST(iTPCI != iTPC);

    // test copy and prefix increment
    ++iTPCI;
    BOOST_TEST(iTPCI == iTPC); // arguable if both are end-iterators by now
  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping cryostat
    geo::TPCID ID(0, 0);
    ID.TPC = geom->NTPC(ID) - 1; // last TPC of first cryostat

    geo::TPC_id_iterator iTPC(ID, policy);
    geo::TPC_iterator const iTPCElem{geom, iTPC};
    BOOST_TEST(iTPCElem);

    // check that the pointed ID is as expected
    BOOST_TEST(*iTPC == ID);
    BOOST_TEST(iTPC->Cryostat == ID.Cryostat);
    BOOST_TEST(iTPC->TPC == ID.TPC);
    BOOST_TEST(iTPCElem.get() == geom->TPCPtr(ID));

    ++iTPC;
    // check that the pointed ID is as expected
    BOOST_TEST(iTPC->Cryostat == geo::CryostatID::CryostatID_t(ID.Cryostat + 1));
    BOOST_TEST(iTPC->TPC == geo::TPCID::TPCID_t(0));

    // test iterator to last TPC
    geo::TPCID LastID(geom->Ncryostats() - 1, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    geo::TPC_id_iterator iLastTPC(LastID, policy);
    geo::TPC_iterator iLastTPCElem{geom, iLastTPC};
    BOOST_TEST_CHECKPOINT("Position-created iterator to last TPC ID: " << std::string(*iLastTPC));

    // check that the iterator tests true
    BOOST_TEST(iLastTPCElem);

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastTPC == LastID);
    BOOST_TEST(iLastTPC->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastTPC->TPC == LastID.TPC);
    BOOST_TEST(iLastTPCElem.get() == geom->TPCPtr(LastID));

    // test increment to past-the-end
    geo::TPC_id_iterator iEndTPC = iLastTPC;
    ++iEndTPC;
    ++iLastTPCElem;

    // check that the iterator tests false
    BOOST_TEST(!iLastTPCElem);

    BOOST_TEST(iEndTPC->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndTPC->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iEndTPC == geom->end<TPCID>());
    BOOST_TEST(!iLastTPCElem.get());
  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::TPC_id_iterator iTPC{policy.GetEndID<geo::TPCID>(), policy};
    geo::TPC_iterator const iTPCElem{geom, iTPC};
    BOOST_TEST_CHECKPOINT("End-created TPC ID iterator: " << std::string(*iTPC));

    BOOST_TEST(iTPC->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iTPC->TPC == geo::TPCID::TPCID_t(0));

    // check that the iterator tests false
    BOOST_TEST(!iTPCElem);

    // check access to geometry element (result of operator* is not defined)
    BOOST_TEST(!(iTPCElem.get())); // should get nullptr

    // construct at end position by geometry
    geo::TPC_id_sentinel iTPCGE = geom->end<TPCID>();
    BOOST_TEST(iTPCGE == iTPC);

    // initialize to the end directly; this has probably ID's isValid true
    geo::TPC_id_iterator iTPC2(geo::TPCID(geom->Ncryostats(), 0), policy);
    BOOST_TEST(iTPC2->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iTPC2->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iTPC2 == iTPC);
  }

} // GeometryIteratorTestAlg::TPCIDIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::TPCIteratorsTest() const
{
  GeometryIterationPolicy const policy{geom};
  /*
   * This test is extensively based on the assumption that the iterators should
   * behave like the corresponding ID iterators, including "corner cases".
   *
   * public interface (geometry_element_iterator<>)
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     geometry_element_iterator()
   *
   *     /// Constructor: points to the same element as the specified ID iterator
   *     geometry_element_iterator(id_iterator_t const& iter)
   *
   *     /// Constructor: points to the same element as the specified ID iterator
   *     geometry_element_iterator(id_iterator_t&& iter)
   *
   *     /// Constructor: points to the specified geometry element
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, GeoID_t const& start_from)
   *
   *     /// Returns true if the two iterators point to the same object
   *     bool operator== (iterator const& as) const
   *
   *     /// Returns true if the two iterators point to the same object
   *     bool operator== (id_iterator_t const& as) const
   *
   *     /// Returns true if the two iterators point to different objects
   *     bool operator!= (iterator const& as) const
   *
   *     /// Returns true if the two iterators point to different objects
   *     bool operator!= (id_iterator_t const& as) const
   *
   *     /// Returns the geometry element the iterator points to
   *     Element_t const& operator* () const
   *
   *     /// Returns a pointer to the element the iterator points to (or nullptr)
   *     Element_t const* operator-> () const
   *
   *     /// Prefix increment: returns this iterator pointing to the next element
   *     iterator& operator++ ()
   *
   *     /// Postfix increment: returns the current iterator, then increments it
   *     iterator operator++ (int)
   *
   *     /// Returns whether the iterator is pointing to a valid geometry element
   *     operator bool() const
   *
   *     /// Returns a pointer to the geometry element, or nullptr if invalid
   *     ElementPtr_t get() const
   *
   *     /// Returns the ID of the pointed geometry element
   *     LocalID_t const& ID() const
   *
   */

  //
  // default constructed
  //
  {
    geo::TPC_id_iterator iTPCID;
    BOOST_TEST_CHECKPOINT("Default created TPC iterator: " << std::string(*iTPCID));

    geo::TPC_iterator iTPC;

    // ID comparison
    BOOST_TEST(iTPC.ID() == *iTPCID);

    // check copy assignment
    geo::TPC_iterator iTPC_copy(iTPC);
    geo::TPC_id_iterator iTPCID_copy(iTPCID);

    // check comparisons too
    BOOST_TEST(iTPC == iTPC_copy);
    BOOST_TEST(iTPC_copy == iTPC);
    BOOST_TEST(!(iTPC != iTPC_copy));
    BOOST_TEST(!(iTPC_copy != iTPC));

    BOOST_TEST(iTPC == iTPC_copy);
    BOOST_TEST(iTPCID_copy == iTPCID);
  }

  //
  // begin-constructed
  //
  {
    auto BeginID = geo::TPCID::first();

    geo::TPC_id_iterator iTPCID(BeginID, policy);

    BOOST_TEST_CHECKPOINT("Begin-created TPC iterator (" << std::string(BeginID) << ")");

    // initialize to the beginning directly
    geo::TPC_iterator iTPCD(geom, iTPCID);

    // construct from explicit begin position
    geo::TPC_iterator iTPCBC{geom, geom->begin<TPCID>()};
    CompareIteratorAndIteratorID(iTPCBC, iTPCID);

    // construct at begin position by geometry
    geo::TPC_iterator iTPCGB = geom->begin<TPCGeo>();
    CompareIteratorAndIteratorID(iTPCGB, iTPCID);
  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last TPC
    geo::TPCID LastID(geom->Ncryostats() - 1, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    geo::TPC_id_iterator iLastTPCID(LastID, policy);

    BOOST_TEST_CHECKPOINT("Position-created iterator to last TPC: " << std::string(LastID));
    geo::TPC_iterator iLastTPC(geom, iLastTPCID);

    // test increment to past-the-end
    geo::TPC_id_iterator iEndTPCID = iLastTPCID;
    ++iEndTPCID;

    geo::TPC_iterator iEndTPC = iLastTPC;
    ++iEndTPC;

    CompareIteratorAndIteratorID(iEndTPC, iEndTPCID);
  }

  //
  // end-constructed
  //
  {
    geo::TPC_id_sentinel iTPCID = geom->end<TPCID>();

    // construct from end position
    BOOST_TEST_CHECKPOINT("End-created TPC iterator");
    geo::TPC_sentinel iTPC{geom->end<TPCID>()};
    BOOST_TEST(iTPC == iTPCID);

    // construct at end position by geometry
    geo::TPC_sentinel iTPCGE = geom->end<TPCGeo>();
    BOOST_TEST(iTPCGE == iTPCID);
  }

} // GeometryIteratorTestAlg::TPCIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::PlaneIDIteratorsTest() const
{
  GeometryIterationPolicy const policy{geom};
  /*
   * public interface (plane_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     plane_id_iterator_base()
   *
   *     /// Constructor: points to the specified cryostat
   *     plane_id_iterator_base
   *       (geo::GeometryCore const* geom, GEOID const& start_from)
   *
   *     // TODO reconsider if the additional template is indeed needed
   *     /// Returns true if the two iterators point to the same plane
   *     template <typename OTHERID>
   *     bool operator== (plane_id_iterator_base<OTHERID> const& as) const
   *
   *     /// Returns true if the two iterators point to different planes
   *     template <typename OTHERID>
   *     bool operator!= (plane_id_iterator_base<OTHERID> const& as) const
   *
   *     /// Returns the PlaneID the iterator points to
   *     LocalID_t const& operator* () const
   *
   *     /// Returns the PlaneID the iterator points to
   *     LocalID_t const* operator-> () const
   *
   *     /// Prefix increment: returns this iterator pointing to the next plane
   *     iterator& operator++ ()
   *
   *     /// Postfix increment: returns the current iterator, then increments it
   *     iterator operator++ (int)
   *
   *     /// Returns whether the iterator is pointing to a valid plane
   *     operator bool() const
   *
   *     /// Returns a pointer to plane, or nullptr if invalid
   *     ElementPtr_t get() const
   *
   */

  //
  // default constructed
  //
  {
    geo::plane_id_iterator iPlane;
    BOOST_TEST_CHECKPOINT("Default created plane ID iterator: " << std::string(*iPlane));

    BOOST_TEST(iPlane->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iPlane->TPC == geo::TPCID::getInvalidID());
    BOOST_TEST(iPlane->Plane == geo::PlaneID::getInvalidID());
  }

  //
  // begin-constructed
  //
  {
    // initialize to the beginning directly; this has probably ID's isValid true
    auto BeginID = geo::PlaneID::first();
    geo::plane_id_iterator iPlane(BeginID, policy);
    BOOST_TEST(iPlane->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iPlane->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iPlane->Plane == geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iPlane == iPlane);

    // construct from explicit begin position
    geo::plane_id_iterator iPlaneBC{BeginID, policy};
    BOOST_TEST(iPlaneBC == iPlane);

    // construct at begin position by geometry
    geo::plane_id_iterator iPlaneGB = geom->begin<PlaneID>();
    BOOST_TEST(iPlaneGB == iPlane);

    // check access to ID
    BOOST_TEST(*iPlane == BeginID);
    BOOST_TEST(iPlane->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iPlane->TPC == BeginID.TPC);
    BOOST_TEST(iPlane->Plane == BeginID.Plane);

    // check access to geometry element
    geo::PlaneGeo const* pPlane = geom->PlanePtr(BeginID);
    geo::plane_iterator const iPlaneElem{geom, iPlane};
    BOOST_TEST(iPlaneElem);
    BOOST_TEST(iPlaneElem.get() == pPlane);

    // test copy and postfix increment
    geo::plane_id_iterator iPlaneI(iPlane++);

    const unsigned int nPlanesInC0T0 = geom->Nplanes(geo::TPCID(0, 0));
    if (nPlanesInC0T0 > 1) {
      BOOST_TEST(iPlaneI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iPlaneI->TPC == geo::TPCID::TPCID_t(0));
      BOOST_TEST(iPlaneI->Plane == geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iPlane->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iPlane->TPC == geo::TPCID::TPCID_t(0));
      BOOST_TEST(iPlane->Plane == geo::PlaneID::PlaneID_t(1));
    }
    BOOST_TEST(iPlaneI != iPlane);

    // test copy and prefix increment
    ++iPlaneI;
    BOOST_TEST(iPlaneI == iPlane); // arguable if both are end-iterators by now
  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping TPC
    geo::PlaneID ID(0, 0, 0);
    ID.Plane = geom->Nplanes(ID) - 1; // last plane of first TPC

    geo::plane_id_iterator iPlane(ID, policy);
    geo::plane_iterator iPlaneElem{geom, iPlane};

    // check that the iterator tests true
    BOOST_TEST(iPlaneElem);

    // check that the pointed ID is as expected
    BOOST_TEST(*iPlane == ID);
    BOOST_TEST(iPlane->Cryostat == ID.Cryostat);
    BOOST_TEST(iPlane->TPC == ID.TPC);
    BOOST_TEST(iPlane->Plane == ID.Plane);
    BOOST_TEST(iPlaneElem.get() == geom->PlanePtr(ID));

    // check that the pointed ID is as expected
    ++iPlane;
    if (ID.TPC + 1 < geom->NTPC(ID)) {
      BOOST_TEST(iPlane->Cryostat == ID.Cryostat);
      BOOST_TEST(iPlane->TPC == ID.TPC + 1);
      BOOST_TEST(iPlane->Plane == geo::PlaneID::PlaneID_t(0));
    }
    else {
      BOOST_TEST(iPlane->Cryostat == ID.Cryostat + 1);
      BOOST_TEST(iPlane->TPC == geo::TPCID::TPCID_t(0));
      BOOST_TEST(iPlane->Plane == geo::PlaneID::PlaneID_t(0));
    }

    // test iterator to last plane
    geo::PlaneID LastID(geom->Ncryostats() - 1, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1;      // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    geo::plane_id_iterator iLastPlane(LastID, policy);
    geo::plane_iterator iLastPlaneElem{geom, iLastPlane};
    BOOST_TEST_CHECKPOINT(
      "Position-created iterator to last plane ID: " << std::string(*iLastPlane));

    // check that the iterator tests true
    BOOST_TEST(iLastPlaneElem);

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastPlane == LastID);
    BOOST_TEST(iLastPlane->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastPlane->TPC == LastID.TPC);
    BOOST_TEST(iLastPlane->Plane == LastID.Plane);
    BOOST_TEST(iLastPlaneElem.get() == geom->PlanePtr(LastID));

    // test increment to past-the-end
    geo::plane_id_iterator iEndPlane = iLastPlane;
    ++iEndPlane;
    ++iLastPlaneElem;

    // check that the iterator tests false
    BOOST_TEST(!iLastPlaneElem);

    BOOST_TEST(iEndPlane->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndPlane->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iEndPlane->Plane == geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iEndPlane == geom->end<PlaneID>());
    BOOST_TEST(!iLastPlaneElem.get());
  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::plane_id_iterator iPlane{policy.GetEndID<geo::PlaneID>(), policy};
    geo::plane_iterator iPlaneElem{geom, iPlane};
    BOOST_TEST_CHECKPOINT("End-created plane ID iterator: " << std::string(*iPlane));

    BOOST_TEST(iPlane->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iPlane->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iPlane->Plane == geo::PlaneID::PlaneID_t(0));

    // check that the iterator tests false
    BOOST_TEST(!iPlaneElem);

    // check access to geometry element (result of operator* is not defined)
    BOOST_TEST(!(iPlaneElem.get())); // should get nullptr

    // construct at end position by geometry
    geo::plane_id_sentinel iPlaneGE = geom->end<PlaneID>();
    BOOST_TEST(iPlaneGE == iPlane);

    // initialize to the end directly; this has probably ID's isValid true
    geo::plane_id_iterator iPlane2(geo::PlaneID(geom->Ncryostats(), 0, 0), policy);
    BOOST_TEST(iPlane2->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iPlane2->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iPlane2->Plane == geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iPlane2 == iPlane);
  }
} // GeometryIteratorTestAlg::PlaneIDIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::PlaneIteratorsTest() const
{
  GeometryIterationPolicy const policy{geom};
  /*
   * This test is extensively based on the assumption that the iterators should
   * behave like the corresponding ID iterators, including "corner cases".
   *
   * public interface (geometry_element_iterator<>)
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     geometry_element_iterator()
   *
   *     /// Constructor: points to the same element as the specified ID iterator
   *     geometry_element_iterator(id_iterator_t const& iter)
   *
   *     /// Constructor: points to the same element as the specified ID iterator
   *     geometry_element_iterator(id_iterator_t&& iter)
   *
   *     /// Constructor: points to the specified geometry element
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, GeoID_t const& start_from)
   *
   *     /// Returns true if the two iterators point to the same object
   *     bool operator== (iterator const& as) const
   *
   *     /// Returns true if the two iterators point to the same object
   *     bool operator== (id_iterator_t const& as) const
   *
   *     /// Returns true if the two iterators point to different objects
   *     bool operator!= (iterator const& as) const
   *
   *     /// Returns true if the two iterators point to different objects
   *     bool operator!= (id_iterator_t const& as) const
   *
   *     /// Returns the geometry element the iterator points to
   *     Element_t const& operator* () const
   *
   *     /// Returns a pointer to the element the iterator points to (or nullptr)
   *     Element_t const* operator-> () const
   *
   *     /// Prefix increment: returns this iterator pointing to the next element
   *     iterator& operator++ ()
   *
   *     /// Postfix increment: returns the current iterator, then increments it
   *     iterator operator++ (int)
   *
   *     /// Returns whether the iterator is pointing to a valid geometry element
   *     operator bool() const
   *
   *     /// Returns a pointer to the geometry element, or nullptr if invalid
   *     ElementPtr_t get() const
   *
   *     /// Returns the ID of the pointed geometry element
   *     LocalID_t const& ID() const
   *
   */

  //
  // default constructed
  //
  {
    geo::plane_id_iterator iPlaneID;
    BOOST_TEST_CHECKPOINT("Default created plane iterator: " << std::string(*iPlaneID));

    geo::plane_iterator iPlane;

    // ID comparison
    BOOST_TEST(iPlane.ID() == *iPlaneID);

    // check copy assignment
    geo::plane_iterator iPlane_copy(iPlane);
    geo::plane_id_iterator iPlaneID_copy(iPlaneID);

    // check comparisons too
    BOOST_TEST(iPlane == iPlane_copy);
    BOOST_TEST(iPlane_copy == iPlane);
    BOOST_TEST(!(iPlane != iPlane_copy));
    BOOST_TEST(!(iPlane_copy != iPlane));

    BOOST_TEST(iPlane == iPlane_copy);
    BOOST_TEST(iPlaneID_copy == iPlaneID);
  }

  //
  // begin-constructed
  //
  {
    auto BeginID = geo::PlaneID::first();

    geo::plane_id_iterator iPlaneID(BeginID, policy);

    BOOST_TEST_CHECKPOINT("Begin-created plane iterator (" << std::string(BeginID) << ")");

    // initialize to the beginning directly
    geo::plane_iterator iPlaneD(geom, iPlaneID);

    // construct from explicit begin position
    geo::plane_iterator iPlaneBC(geom, geom->begin<PlaneID>());
    CompareIteratorAndIteratorID(iPlaneBC, iPlaneID);

    // construct at begin position by geometry
    geo::plane_iterator iPlaneGB = geom->begin<PlaneGeo>();
    CompareIteratorAndIteratorID(iPlaneGB, iPlaneID);
  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last plane
    geo::PlaneID LastID(geom->Ncryostats() - 1, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1;      // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    geo::plane_id_iterator iLastPlaneID(LastID, policy);

    BOOST_TEST_CHECKPOINT("Position-created iterator to last plane: " << std::string(LastID));
    geo::plane_iterator iLastPlane(geom, iLastPlaneID);

    // test increment to past-the-end
    geo::plane_id_iterator iEndPlaneID = iLastPlaneID;
    ++iEndPlaneID;

    geo::plane_iterator iEndPlane = iLastPlane;
    ++iEndPlane;

    CompareIteratorAndIteratorID(iEndPlane, iEndPlaneID);
  }

  //
  // end-constructed
  //
  {
    geo::plane_id_sentinel iPlaneID = geom->end<PlaneID>();

    // construct from end position
    BOOST_TEST_CHECKPOINT("End-created plane iterator");
    geo::plane_sentinel iPlane(geom->end<PlaneID>());
    BOOST_TEST(iPlane == iPlaneID);

    // construct at end position by geometry
    geo::plane_sentinel iPlaneGE = geom->end<PlaneGeo>();
    BOOST_TEST(iPlaneGE == iPlaneID);
  }

} // GeometryIteratorTestAlg::PlaneIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::WireIDIteratorsTest() const
{
  GeometryIterationPolicy const policy{geom};
  /*
   * public interface (wire_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     wire_id_iterator_base()
   *
   *     /// Constructor: points to the specified cryostat
   *     wire_id_iterator_base
   *       (geo::GeometryCore const* geom, GEOID const& start_from)
   *
   *     // TODO reconsider if the additional template is indeed needed
   *     /// Returns true if the two iterators point to the same wire
   *     template <typename OTHERID>
   *     bool operator== (wire_id_iterator_base<OTHERID> const& as) const
   *
   *     /// Returns true if the two iterators point to different wires
   *     template <typename OTHERID>
   *     bool operator!= (wire_id_iterator_base<OTHERID> const& as) const
   *
   *     /// Returns the WireID the iterator points to
   *     LocalID_t const& operator* () const
   *
   *     /// Returns the WireID the iterator points to
   *     LocalID_t const* operator-> () const
   *
   *     /// Prefix increment: returns this iterator pointing to the next wire
   *     iterator& operator++ ()
   *
   *     /// Postfix increment: returns the current iterator, then increments it
   *     iterator operator++ (int)
   *
   *     /// Returns whether the iterator is pointing to a valid wire
   *     operator bool() const
   *
   *     /// Returns a pointer to wire, or nullptr if invalid
   *     ElementPtr_t get() const
   *
   */

  //
  // default constructed
  //
  {
    geo::wire_id_iterator iWire;
    BOOST_TEST_CHECKPOINT("Default created wire ID iterator: " << std::string(*iWire));

    BOOST_TEST(iWire->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iWire->TPC == geo::TPCID::getInvalidID());
    BOOST_TEST(iWire->Plane == geo::PlaneID::getInvalidID());
    BOOST_TEST(iWire->Wire == geo::WireID::getInvalidID());
  }

  //
  // begin-constructed
  //
  {
    // initialize to the beginning directly; this has probably ID's isValid true
    auto BeginID = geo::WireID::first();
    geo::wire_id_iterator iWire(BeginID, policy);
    BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iWire->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iWire->Plane == geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iWire->Wire == geo::WireID::WireID_t(0));

    // construct from explicit begin position
    auto iWireBC = geom->begin<WireID>();
    BOOST_TEST(iWireBC == iWire);

    // construct at begin position by geometry
    geo::wire_id_iterator iWireGB = geom->begin<WireID>();
    BOOST_TEST(iWireGB == iWire);

    // check access to ID
    BOOST_TEST(*iWire == BeginID);
    BOOST_TEST(iWire->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iWire->TPC == BeginID.TPC);
    BOOST_TEST(iWire->Plane == BeginID.Plane);
    BOOST_TEST(iWire->Wire == BeginID.Wire);

    // check access to geometry element
    geo::WireGeo const* pWire = geom->WirePtr(BeginID);
    geo::wire_iterator const iWireElem{geom, iWire};
    BOOST_TEST(iWireElem);
    BOOST_TEST(iWireElem.get() == pWire);

    // test copy and postfix increment
    geo::wire_id_iterator iWireI(iWire++);

    const unsigned int nWiresInC0T0P0 = geom->Nwires(geo::PlaneID(0, 0, 0));
    if (nWiresInC0T0P0 > 1) {
      BOOST_TEST(iWireI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iWireI->TPC == geo::TPCID::TPCID_t(0));
      BOOST_TEST(iWireI->Plane == geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iWireI->Wire == geo::WireID::WireID_t(0));
      BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iWire->TPC == geo::TPCID::TPCID_t(0));
      BOOST_TEST(iWire->Plane == geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iWire->Wire == geo::WireID::WireID_t(1));
    }
    BOOST_TEST(iWireI != iWire);

    // test copy and prefix increment
    ++iWireI;
    BOOST_TEST(iWireI == iWire); // arguable if both are end-iterators by now
  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping plane
    geo::WireID ID(0, 0, 0, 0);
    ID.Wire = geom->Nwires(ID) - 1; // last wire of first plane

    geo::wire_id_iterator iWire(ID, policy);
    geo::wire_iterator iWireElem{geom, iWire};

    // check that the iterator tests true
    BOOST_TEST(iWireElem);

    // check that the pointed ID is as expected
    BOOST_TEST(*iWire == ID);
    BOOST_TEST(iWire->Cryostat == ID.Cryostat);
    BOOST_TEST(iWire->TPC == ID.TPC);
    BOOST_TEST(iWire->Plane == ID.Plane);
    BOOST_TEST(iWire->Wire == ID.Wire);
    BOOST_TEST(iWireElem.get() == geom->WirePtr(ID));

    ++iWire;
    // check that the pointed ID is as expected
    if (ID.Plane + 1 < geom->Nplanes(ID)) {
      BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iWire->TPC == geo::TPCID::TPCID_t(0));
      BOOST_TEST(iWire->Plane == ID.Plane + 1);
      BOOST_TEST(iWire->Wire == geo::WireID::WireID_t(0));
    }
    else if (ID.TPC + 1 < geom->NTPC(ID)) {
      BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iWire->TPC == ID.TPC + 1);
      BOOST_TEST(iWire->Plane == geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iWire->Wire == geo::WireID::WireID_t(0));
    }
    else {
      BOOST_TEST(iWire->Cryostat == ID.Cryostat + 1);
      BOOST_TEST(iWire->TPC == geo::TPCID::TPCID_t(0));
      BOOST_TEST(iWire->Plane == geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iWire->Wire == geo::WireID::WireID_t(0));
    }

    // test iterator to last wire
    geo::WireID LastID(geom->Ncryostats() - 1, 0, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1;      // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    LastID.Wire = geom->Nwires(LastID) - 1;   // last wire of last plane
    geo::wire_id_iterator iLastWire(LastID, policy);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last wire ID: " << std::string(*iLastWire));

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastWire == LastID);
    BOOST_TEST(iLastWire->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastWire->TPC == LastID.TPC);
    BOOST_TEST(iLastWire->Plane == LastID.Plane);
    BOOST_TEST(iLastWire->Wire == LastID.Wire);
    geo::wire_iterator iLastWireElem{geom, iLastWire};
    BOOST_TEST(iLastWireElem);
    BOOST_TEST(iLastWireElem.get() == geom->WirePtr(LastID));

    // test increment to past-the-end
    geo::wire_id_iterator iEndWire = iLastWire;
    ++iEndWire;
    ++iLastWireElem;

    // check that the iterator tests false
    BOOST_TEST(!iLastWireElem);

    BOOST_TEST(iEndWire->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndWire->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iEndWire->Plane == geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iEndWire->Wire == geo::WireID::WireID_t(0));
    BOOST_TEST(iEndWire == geom->end<WireID>());
    BOOST_TEST(!iLastWireElem.get());
  }

  //
  // end-constructed
  //
  {
    // construct from end position
    auto iWire = geom->end<WireID>();
    geo::wire_sentinel const iWireElem [[maybe_unused]]{iWire};
    BOOST_TEST_CHECKPOINT("End-created end ID iterator");

    // construct at end position by geometry
    geo::wire_id_sentinel iWireGE = geom->end<WireID>();
    BOOST_TEST(iWireGE == iWire);

    // initialize to the end directly; this has probably ID's isValid true
    geo::wire_id_iterator iWire2(geo::WireID(geom->Ncryostats(), 0, 0, 0), policy);
    BOOST_TEST(iWire2->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iWire2->TPC == geo::TPCID::TPCID_t(0));
    BOOST_TEST(iWire2->Plane == geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iWire2->Wire == geo::WireID::WireID_t(0));
    BOOST_TEST(iWire2 == iWire);
  }

} // GeometryIteratorTestAlg::WireIDIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::WireIteratorsTest() const
{
  GeometryIterationPolicy const policy{geom};
  /*
   * This test is extensively based on the assumption that the iterators should
   * behave like the corresponding ID iterators, including "corner cases".
   *
   * public interface (geometry_element_iterator<>)
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     geometry_element_iterator()
   *
   *     /// Constructor: points to the same element as the specified ID iterator
   *     geometry_element_iterator(id_iterator_t const& iter)
   *
   *     /// Constructor: points to the same element as the specified ID iterator
   *     geometry_element_iterator(id_iterator_t&& iter)
   *
   *     /// Constructor: points to the specified geometry element
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, GeoID_t const& start_from)
   *
   *     /// Returns true if the two iterators point to the same object
   *     bool operator== (iterator const& as) const
   *
   *     /// Returns true if the two iterators point to the same object
   *     bool operator== (id_iterator_t const& as) const
   *
   *     /// Returns true if the two iterators point to different objects
   *     bool operator!= (iterator const& as) const
   *
   *     /// Returns true if the two iterators point to different objects
   *     bool operator!= (id_iterator_t const& as) const
   *
   *     /// Returns the geometry element the iterator points to
   *     Element_t const& operator* () const
   *
   *     /// Returns a pointer to the element the iterator points to (or nullptr)
   *     Element_t const* operator-> () const
   *
   *     /// Prefix increment: returns this iterator pointing to the next element
   *     iterator& operator++ ()
   *
   *     /// Postfix increment: returns the current iterator, then increments it
   *     iterator operator++ (int)
   *
   *     /// Returns whether the iterator is pointing to a valid geometry element
   *     operator bool() const
   *
   *     /// Returns a pointer to the geometry element, or nullptr if invalid
   *     ElementPtr_t get() const
   *
   *     /// Returns the ID of the pointed geometry element
   *     LocalID_t const& ID() const
   *
   */

  //
  // default constructed
  //
  {
    geo::wire_id_iterator iWireID;
    BOOST_TEST_CHECKPOINT("Default created wire iterator: " << std::string(*iWireID));

    geo::wire_iterator iWire;

    // ID comparison
    BOOST_TEST(iWire.ID() == *iWireID);

    // check copy assignment
    geo::wire_iterator iWire_copy(iWire);
    geo::wire_id_iterator iWireID_copy(iWireID);

    // check comparisons too
    BOOST_TEST(iWire == iWire_copy);
    BOOST_TEST(iWire_copy == iWire);
    BOOST_TEST(!(iWire != iWire_copy));
    BOOST_TEST(!(iWire_copy != iWire));

    BOOST_TEST(iWire == iWire_copy);
    BOOST_TEST(iWireID_copy == iWireID);
  }

  //
  // begin-constructed
  //
  {
    auto BeginID = geo::WireID::first();

    geo::wire_id_iterator iWireID(BeginID, policy);

    BOOST_TEST_CHECKPOINT("Begin-created wire iterator (" << std::string(BeginID) << ")");

    // initialize to the beginning directly
    geo::wire_iterator iWireD(geom, iWireID);

    // construct from explicit begin position
    geo::wire_iterator iWireBC(geom, geom->begin<WireID>());
    CompareIteratorAndIteratorID(iWireBC, iWireID);

    // construct at begin position by geometry
    geo::wire_iterator iWireGB = geom->begin<WireGeo>();
    CompareIteratorAndIteratorID(iWireGB, iWireID);
  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last wire
    geo::WireID LastID(geom->Ncryostats() - 1, 0, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1;      // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    LastID.Wire = geom->Nwires(LastID) - 1;   // last wire of last plane
    geo::wire_id_iterator iLastWireID(LastID, policy);

    BOOST_TEST_CHECKPOINT("Position-created iterator to last wire: " << std::string(LastID));
    geo::wire_iterator iLastWire(geom, iLastWireID);

    // test increment to past-the-end
    geo::wire_id_iterator iEndWireID = iLastWireID;
    ++iEndWireID;

    geo::wire_iterator iEndWire = iLastWire;
    ++iEndWire;

    CompareIteratorAndIteratorID(iEndWire, iEndWireID);
  }

  //
  // end-constructed
  //
  {
    geo::wire_id_sentinel iWireID = geom->end<WireID>();

    // construct from end position
    BOOST_TEST_CHECKPOINT("End-created wire iterator");
    geo::wire_sentinel iWire{policy.GetEndID<geo::WireID>()};
    BOOST_TEST(iWire == iWireID);

    // construct at end position by geometry
    geo::wire_sentinel iWireGE = geom->end<WireGeo>();
    BOOST_TEST(iWireGE == iWireID);
  }

} // GeometryIteratorTestAlg::WireIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::TPCsetIDIteratorsTest() const
{
  ReadoutIterationPolicy const policy{geom, channelMapAlg};
  /*
   * public interface (TPCset_id_iterator_base):
   *
   *
   *   /// Default constructor; effect not defined: assign to it before using!
   *   TPCset_id_iterator_base()
   *
   *   /// Constructor: points to the specified cryostat.
   *   TPCset_id_iterator_base
   *     (geo::GeometryCore const* geom, GeoID_t const& start_from)
   *
   *   /// Returns true if the two iterators point to the same TPC set.
   *   template <typename OTHERID>
   *   bool operator== (TPCset_id_iterator_base<OTHERID> const& as) const
   *
   *   /// Returns true if the two iterators point to different TPC sets.
   *   template <typename OTHERID>
   *   bool operator!= (TPCset_id_iterator_base<OTHERID> const& as) const
   *
   *   /// Returns the TPCsetID the iterator points to.
   *   LocalID_t const& operator* () const
   *
   *   /// Returns the TPCsetID the iterator points to.
   *   LocalID_t const* operator-> () const
   *
   *   /// Prefix increment: returns this iterator pointing to the next TPC set.
   *   iterator& operator++ ()
   *
   *   /// Postfix increment: returns the current iterator, then increments it.
   *   iterator operator++ (int)
   *
   *   /// Returns whether the iterator is pointing to a valid TPC set.
   *   operator bool() const
   *
   *
   */

  //
  // default constructed
  //
  {
    geo::TPCset_id_iterator iTPCset;
    BOOST_TEST_CHECKPOINT("Default created TPC set ID iterator: " << std::string(*iTPCset));

    BOOST_TEST(iTPCset->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iTPCset->TPCset == readout::TPCsetID::getInvalidID());
  }

  //
  // begin-constructed
  //
  {
    // initialize to the beginning directly; this has probably ID's isValid true
    auto BeginID = readout::TPCsetID::first();
    geo::TPCset_id_iterator iTPCset(BeginID, policy);
    BOOST_TEST(iTPCset->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iTPCset->TPCset == readout::TPCsetID::TPCsetID_t(0));

    // construct from explicit begin position
    geo::TPCset_id_iterator iTPCsetBC{BeginID, policy};
    BOOST_TEST(iTPCsetBC == iTPCset);

    // construct at begin position by geometry
    geo::TPCset_id_iterator iTPCsetGB = channelMapAlg->begin<readout::TPCsetID>();
    BOOST_TEST(iTPCsetGB == iTPCset);

    // check access to ID
    BOOST_TEST(*iTPCset == BeginID);
    BOOST_TEST(iTPCset->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iTPCset->TPCset == BeginID.TPCset);

    // test copy and postfix increment
    geo::TPCset_id_iterator iTPCsetI(iTPCset++);

    const unsigned int nTPCsetsInC0 = channelMapAlg->NTPCsets(geo::CryostatID(0));
    if (nTPCsetsInC0 > 1) {
      BOOST_TEST(iTPCsetI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iTPCsetI->TPCset == readout::TPCsetID::TPCsetID_t(0));
      BOOST_TEST(iTPCset->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iTPCset->TPCset == readout::TPCsetID::TPCsetID_t(1));
    }
    BOOST_TEST(iTPCsetI != iTPCset);

    // test copy and prefix increment
    ++iTPCsetI;
    BOOST_TEST(iTPCsetI == iTPCset); // arguable if both are end-iterators by now
  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping cryostat
    readout::TPCsetID ID(0, 0);
    ID.TPCset = channelMapAlg->NTPCsets(ID) - 1; // last TPC set of first cryostat

    geo::TPCset_id_iterator iTPCset(ID, policy);

    // check that the pointed ID is as expected
    BOOST_TEST(*iTPCset == ID);
    BOOST_TEST(iTPCset->Cryostat == ID.Cryostat);
    BOOST_TEST(iTPCset->TPCset == ID.TPCset);

    ++iTPCset;
    // check that the pointed ID is as expected
    BOOST_TEST(iTPCset->Cryostat == geo::CryostatID::CryostatID_t(ID.Cryostat + 1));
    BOOST_TEST(iTPCset->TPCset == readout::TPCsetID::TPCsetID_t(0));

    // test iterator to last TPC
    readout::TPCsetID LastID(geom->Ncryostats() - 1, 0);
    LastID.TPCset = channelMapAlg->NTPCsets(LastID) - 1; // last TPC set of last cryostat
    geo::TPCset_id_iterator iLastTPCset(LastID, policy);
    BOOST_TEST_CHECKPOINT(
      "Position-created iterator to last TPC set ID: " << std::string(*iLastTPCset));

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastTPCset == LastID);
    BOOST_TEST(iLastTPCset->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastTPCset->TPCset == LastID.TPCset);

    // test increment to past-the-end
    geo::TPCset_id_iterator iEndTPCset = iLastTPCset;
    ++iEndTPCset;

    BOOST_TEST(iEndTPCset->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndTPCset->TPCset == readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iEndTPCset == channelMapAlg->end<readout::TPCsetID>());
  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::TPCset_id_iterator iTPCset{policy.GetEndID<readout::TPCsetID>(), policy};
    BOOST_TEST_CHECKPOINT("End-created TPC set ID iterator: " << std::string(*iTPCset));

    BOOST_TEST(iTPCset->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iTPCset->TPCset == readout::TPCsetID::TPCsetID_t(0));

    // construct at end position by geometry
    BOOST_TEST(iTPCset == channelMapAlg->end<readout::TPCsetID>());

    // initialize to the end directly; this has probably ID's isValid true
    geo::TPCset_id_iterator iTPCset2(readout::TPCsetID(geom->Ncryostats(), 0), policy);
    BOOST_TEST(iTPCset2->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iTPCset2->TPCset == readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iTPCset2 == iTPCset);
  }

} // GeometryIteratorTestAlg::TPCsetIDIteratorsTest()

//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::ROPIDIteratorsTest() const
{
  ReadoutIterationPolicy const policy{geom, channelMapAlg};
  /*
   * public interface (ROP_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     ROP_id_iterator_base()
   *
   *     /// Constructor: points to the specified readout plane.
   *     ROP_id_iterator_base
   *       (geo::GeometryCore const* geom, GeoID_t const& start_from)
   *
   *     /// Returns true if the two iterators point to the same readout plane.
   *     template <typename OTHERID>
   *     bool operator== (ROP_id_iterator_base<OTHERID> const& as) const
   *
   *     /// Returns true if the two iterators point to different readout planes.
   *     template <typename OTHERID>
   *     bool operator!= (ROP_id_iterator_base<OTHERID> const& as) const
   *
   *     /// Returns the PlaneID the iterator points to
   *     LocalID_t const& operator* ()
   *
   *     /// Returns the PlaneID the iterator points to
   *     LocalID_t const* operator-> () const
   *
   *     /// Prefix increment: returns this iterator pointing to the next plane
   *     iterator& operator++ ()
   *
   *     /// Postfix increment: returns the current iterator, then increments it.
   *     iterator operator++ (int)
   *
   *     /// Returns whether the iterator is pointing to a valid plane.
   *     operator bool() const
   *
   */

  //
  // default constructed
  //
  {
    geo::ROP_id_iterator iROP;
    BOOST_TEST_CHECKPOINT("Default created readout plane ID iterator: " << std::string(*iROP));

    BOOST_TEST(iROP->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iROP->TPCset == readout::TPCsetID::getInvalidID());
    BOOST_TEST(iROP->ROP == readout::ROPID::getInvalidID());
  }

  //
  // begin-constructed
  //
  {
    // initialize to the beginning directly; this has probably ID's isValid true
    auto BeginID = readout::ROPID::first();
    geo::ROP_id_iterator iROP(BeginID, policy);
    BOOST_TEST(iROP->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iROP->TPCset == readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iROP->ROP == readout::ROPID::ROPID_t(0));

    // construct at begin position by geometry
    geo::ROP_id_iterator iROPGB = channelMapAlg->begin<readout::ROPID>();
    BOOST_TEST(iROPGB == iROP);

    // check access to ID
    BOOST_TEST(*iROP == BeginID);
    BOOST_TEST(iROP->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iROP->TPCset == BeginID.TPCset);
    BOOST_TEST(iROP->ROP == BeginID.ROP);

    // test copy and postfix increment
    geo::ROP_id_iterator iROPI(iROP++);

    const unsigned int nReadoutPlanesInC0S0 = channelMapAlg->NROPs(readout::TPCsetID(0, 0));
    if (nReadoutPlanesInC0S0 > 1) {
      BOOST_TEST(iROPI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iROPI->TPCset == readout::TPCsetID::TPCsetID_t(0));
      BOOST_TEST(iROPI->ROP == readout::ROPID::ROPID_t(0));
      BOOST_TEST(iROP->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iROP->TPCset == readout::TPCsetID::TPCsetID_t(0));
      BOOST_TEST(iROP->ROP == readout::ROPID::ROPID_t(1));
    }
    BOOST_TEST(iROPI != iROP);

    // test copy and prefix increment
    ++iROPI;
    BOOST_TEST(iROPI == iROP); // arguable if both are end-iterators by now
  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping TPC
    readout::ROPID ID(0, 0, 0);
    ID.ROP = channelMapAlg->NROPs(ID) - 1; // last plane of first TPC set

    geo::ROP_id_iterator iROP(ID, policy);

    // check that the pointed ID is as expected
    BOOST_TEST(*iROP == ID);
    BOOST_TEST(iROP->Cryostat == ID.Cryostat);
    BOOST_TEST(iROP->TPCset == ID.TPCset);
    BOOST_TEST(iROP->ROP == ID.ROP);

    // check that the pointed ID is as expected
    ++iROP;
    if (ID.TPCset + 1 < (int)channelMapAlg->NTPCsets(ID)) {
      BOOST_TEST(iROP->Cryostat == ID.Cryostat);
      BOOST_TEST(iROP->TPCset == ID.TPCset + 1);
      BOOST_TEST(iROP->ROP == readout::ROPID::ROPID_t(0));
    }
    else {
      BOOST_TEST(iROP->Cryostat == ID.Cryostat + 1);
      BOOST_TEST(iROP->TPCset == readout::TPCsetID::TPCsetID_t(0));
      BOOST_TEST(iROP->ROP == readout::ROPID::ROPID_t(0));
    }

    // test iterator to last plane
    readout::ROPID LastID(geom->Ncryostats() - 1, 0, 0);
    LastID.TPCset = channelMapAlg->NTPCsets(LastID) - 1; // last TPC set of last cryostat
    LastID.ROP = channelMapAlg->NROPs(LastID) - 1;       // last readout plane of last TPC set
    geo::ROP_id_iterator iLastROP(LastID, policy);
    BOOST_TEST_CHECKPOINT(
      "Position-created iterator to last readout plane ID: " << std::string(*iLastROP));

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastROP == LastID);
    BOOST_TEST(iLastROP->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastROP->TPCset == LastID.TPCset);
    BOOST_TEST(iLastROP->ROP == LastID.ROP);

    // test increment to past-the-end
    geo::ROP_id_iterator iEndROP = iLastROP;
    ++iEndROP;

    BOOST_TEST(iEndROP->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndROP->TPCset == readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iEndROP->ROP == readout::ROPID::ROPID_t(0));
    BOOST_TEST(iEndROP == channelMapAlg->end<readout::ROPID>());
  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::ROP_id_iterator iROP{policy.GetEndID<readout::ROPID>(), policy};
    BOOST_TEST_CHECKPOINT("End-created readout plane ID iterator: " << std::string(*iROP));

    // construct at end position by geometry
    BOOST_TEST(iROP == channelMapAlg->end<readout::ROPID>());

    // initialize to the end directly; this has probably ID's isValid true
    geo::ROP_id_iterator iROP2(readout::ROPID(geom->Ncryostats(), 0, 0), policy);
    BOOST_TEST(iROP->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iROP->TPCset == readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iROP->ROP == readout::ROPID::ROPID_t(0));
    BOOST_TEST(iROP2 == iROP);
  }
} // GeometryIteratorTestAlg::ROPIDIteratorsTest()
