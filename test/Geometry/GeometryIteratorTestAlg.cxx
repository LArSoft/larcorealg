/**
 * @file   GeometryIteratorTestAlg.cxx
 * @brief  Unit test for geometry iterators
 * @date   May 7th, 2015
 * @author petrillo@fnal.gov
 *
 * The methods require a Boost test enviroment.
 */

// LArSoft libraries
#include "test/Geometry/GeometryIteratorTestAlg.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/GeometryCore.h"

// Boost libraries
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <string>
#include <type_traits>


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
  void CompareIteratorAndIteratorID(ITER iter, ITERID id_iter) {

      static_assert(
        std::is_same<typename ITER::id_iterator_t, ITERID>::value,
        "CompareIteratorAndIteratorID() requires compatible iterator types"
        );

      // direct comparison
      BOOST_CHECK_EQUAL(iter, id_iter);
      BOOST_CHECK(! (iter != id_iter));

      // ID comparison
      BOOST_CHECK_EQUAL(iter.ID(), *id_iter);

      auto pGeoElement = id_iter.get();

      // dereference
      BOOST_CHECK_EQUAL(iter.get(), pGeoElement);
      BOOST_CHECK_EQUAL(iter.operator->(), pGeoElement);

      if (pGeoElement) BOOST_CHECK_EQUAL(&*iter, pGeoElement);
      else             BOOST_CHECK_THROW(*iter, cet::exception);

      // boolean conversions
      BOOST_CHECK_EQUAL(bool(iter), bool(id_iter));

      // check copy assignment
      ITER iter_copy(iter);
      ITERID id_iter_copy(id_iter);

      // check comparisons too
      BOOST_CHECK_EQUAL(iter, iter_copy);
      BOOST_CHECK_EQUAL(iter_copy, iter);
      BOOST_CHECK(!(iter != iter_copy));
      BOOST_CHECK(!(iter_copy != iter));

      // check increment operator
      BOOST_CHECK_EQUAL(iter++, id_iter++);
      BOOST_CHECK_EQUAL(++iter_copy, ++id_iter_copy);

      BOOST_CHECK_EQUAL(iter, iter_copy);

  } // CompareIteratorAndIteratorID()

} // local namespace


//-----------------------------------------------------------------------------
unsigned int geo::GeometryIteratorTestAlg::Run() const {
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
void geo::GeometryIteratorTestAlg::CryostatIDIteratorsTest() const {

  /*
   * public interface (cryostat_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     cryostat_id_iterator_base();
   *
   *     /// Constructor: points to begin
   *     cryostat_id_iterator_base(geo::GeometryCore const* geom);
   *
   *     /// Constructor: points to the specified cryostat
   *     cryostat_id_iterator_base
   *       (geo::GeometryCore const* geom, GEOID const& start_from);
   *
   *     /// Constructor: points to begin
   *     cryostat_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t);
   *
   *     /// Constructor: points to end
   *     cryostat_id_iterator_base(geo::GeometryCore const* geom, EndPos_t);
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
    geo::GeometryCore::cryostat_id_iterator iCryo;
    BOOST_TEST_CHECKPOINT
      ("Default created cryostat ID iterator: " << std::string(*iCryo));

    BOOST_CHECK_EQUAL(iCryo->Cryostat, geo::CryostatID::getInvalidID());

    // check that the iterator tests false
    BOOST_CHECK(!iCryo);
    BOOST_CHECK(!(bool(iCryo)));

  /*
    // this is more of a curiosity, since this behaviour is undefined
    ++iCryo;
    BOOST_TEST_CHECKPOINT("  (incremented: " << std::string(*iCryo) << ")");
  */
  }

  //
  // begin-constructed
  //
  {
    geo::GeometryCore::cryostat_id_iterator iCryo(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created cryostat ID iterator: " << std::string(*iCryo));

    BOOST_CHECK_EQUAL(iCryo->Cryostat, geo::CryostatID::CryostatID_t(0));

    // check that the iterator tests true
    BOOST_CHECK(bool(iCryo));
    BOOST_CHECK(!!iCryo);

    // initialize to the beginning directly; this has probably ID's isValid true
    geo::CryostatID BeginID;
    geom->GetBeginID(BeginID);
    geo::GeometryCore::cryostat_id_iterator iCryoD(geom, BeginID);
    BOOST_CHECK_EQUAL(iCryoD->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iCryoD, iCryo);

    // construct from explicit begin position
    geo::GeometryCore::cryostat_id_iterator iCryoBC
      (geom, geo::GeometryCore::cryostat_id_iterator::begin_pos);
    BOOST_CHECK_EQUAL(iCryoBC, iCryo);

    // construct at begin position by geometry
    geo::GeometryCore::cryostat_id_iterator iCryoGB = geom->begin_cryostat_id();
    BOOST_CHECK_EQUAL(iCryoGB, iCryo);

    // check access to ID
    BOOST_CHECK_EQUAL(*iCryo, BeginID);
    BOOST_CHECK_EQUAL(iCryo->Cryostat, BeginID.Cryostat);

    // check access to geometry element
    geo::CryostatGeo const* pCryo = geom->CryostatPtr(BeginID);
    BOOST_CHECK_EQUAL(iCryo.get(), pCryo);

    // test copy and postfix increment
    geo::GeometryCore::cryostat_id_iterator iCryoI(iCryo++);

    BOOST_CHECK_EQUAL(iCryo->Cryostat, geo::CryostatID::CryostatID_t(1));
    BOOST_CHECK_EQUAL(iCryoI->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_NE(iCryoI, iCryo);

    // test copy and prefix increment
    ++iCryoI;
    BOOST_CHECK_EQUAL(iCryoI->Cryostat, geo::CryostatID::CryostatID_t(1));
    BOOST_CHECK_EQUAL(iCryoI, iCryo);

    if (geom->Ncryostats() > 1) {
      ++iCryoI;
      BOOST_CHECK_EQUAL(iCryoI->Cryostat, geo::CryostatID::CryostatID_t(2));
      BOOST_CHECK_NE(iCryoI, iCryo);
    }

  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last TPC
    geo::CryostatID LastID(geom->Ncryostats() - 1); // last cryostat

    geo::GeometryCore::cryostat_id_iterator iLastCryo(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last cryostat ID: "
      << std::string(*iLastCryo));

    // check that the iterator tests true
    BOOST_CHECK(bool(iLastCryo));
    BOOST_CHECK(!!iLastCryo);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iLastCryo, LastID);
    BOOST_CHECK_EQUAL(iLastCryo->Cryostat, LastID.Cryostat);
    BOOST_CHECK_EQUAL(iLastCryo.get(), geom->CryostatPtr(LastID));

    // test increment to past-the-end
    geo::GeometryCore::cryostat_id_iterator iEndCryo = iLastCryo;
    ++iEndCryo;

    // check that the iterator tests false
    BOOST_CHECK(!bool(iEndCryo));
    BOOST_CHECK(!iEndCryo);

    BOOST_CHECK_EQUAL(iEndCryo->Cryostat, geom->Ncryostats());
    BOOST_CHECK_EQUAL(iEndCryo, geom->end_cryostat_id());
    BOOST_CHECK(!iEndCryo.get());

  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::GeometryCore::cryostat_id_iterator iCryo
      (geom, geo::GeometryCore::cryostat_id_iterator::end_pos);
    BOOST_TEST_CHECKPOINT
      ("End-created cryostat ID iterator: " << std::string(*iCryo));

    BOOST_CHECK_EQUAL(iCryo->Cryostat, geom->Ncryostats());

    // check that the iterator tests false
    BOOST_CHECK(!bool(iCryo));
    BOOST_CHECK(!iCryo);

    // check access to geometry element (result of operator* is not defined)
    BOOST_CHECK(!(iCryo.get())); // should get nullptr

    // construct at end position by geometry
    geo::GeometryCore::cryostat_id_iterator iCryoGE = geom->end_cryostat_id();
    BOOST_CHECK_EQUAL(iCryoGE, iCryo);

    // initialize to the end directly; this has probably ID's isValid true
    geo::GeometryCore::cryostat_id_iterator iCryo2
      (geom, geo::CryostatID(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iCryo2->Cryostat, geom->Ncryostats());
    BOOST_CHECK_EQUAL(iCryo2, iCryo);
  /*
    // this is more of a curiosity, since this behaviour is undefined
    ++iCryo;
    BOOST_TEST_CHECKPOINT("  (incremented: " << std::string(*iCryo) << ")");
  */

  }

} // GeometryIteratorTestAlg::CryostatIDIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::CryostatIteratorsTest() const {

  /*
   * This test is extensively based on the assumption that the iterators should
   * behave like the corresponding ID iterators, including "corner cases".
   *
   * public interface (geometry_element_iterator<>)
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     geometry_element_iterator()
   *
   *     /// Constructor: points to begin
   *     geometry_element_iterator(geo::GeometryCore const* geom)
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
   *     /// Constructor: points to beginning
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, BeginPos_t const pos)
   *
   *     /// Constructor: points to end
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, EndPos_t const pos)
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
    geo::GeometryCore::cryostat_id_iterator iCryoID;
    BOOST_TEST_CHECKPOINT
      ("Default created cryostat iterator: " << std::string(*iCryoID));

    geo::GeometryCore::cryostat_iterator iCryo;

    // direct comparison
    BOOST_CHECK_EQUAL(iCryo, iCryoID);
    BOOST_CHECK(!(iCryo != iCryoID));

    // ID comparison
    BOOST_CHECK_EQUAL(iCryo.ID(), *iCryoID);

    // check copy assignment
    geo::GeometryCore::cryostat_iterator iCryo_copy(iCryo);
    geo::GeometryCore::cryostat_id_iterator iCryoID_copy(iCryoID);

    // check comparisons too
    BOOST_CHECK_EQUAL(iCryo, iCryo_copy);
    BOOST_CHECK_EQUAL(iCryo_copy, iCryo);
    BOOST_CHECK(!(iCryo != iCryo_copy));
    BOOST_CHECK(!(iCryo_copy != iCryo));

    BOOST_CHECK_EQUAL(iCryo, iCryo_copy);

  }

  //
  // begin-constructed
  //
  {
    geo::CryostatID BeginID;
    geom->GetBeginID(BeginID);

    geo::GeometryCore::cryostat_id_iterator iCryoID(geom, BeginID);


    BOOST_TEST_CHECKPOINT
      ("Begin-created cryostat iterator (" << std::string(BeginID) << ")");

    // initialize to the beginning directly
    geo::GeometryCore::cryostat_iterator iCryoD(geom, BeginID);
    CompareIteratorAndIteratorID(iCryoD, iCryoID);

    // initialize to the beginning with implicit initialization
    geo::GeometryCore::cryostat_iterator iCryo(geom);
    CompareIteratorAndIteratorID(iCryo, iCryoID);

    // construct from explicit begin position
    geo::GeometryCore::cryostat_iterator iCryoBC
      (geom, geo::GeometryCore::cryostat_iterator::begin_pos);
    CompareIteratorAndIteratorID(iCryoBC, iCryoID);

    // construct at begin position by geometry
    geo::GeometryCore::cryostat_iterator iCryoGB = geom->begin_cryostat();
    CompareIteratorAndIteratorID(iCryoGB, iCryoID);

  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last TPC
    geo::CryostatID LastID(geom->Ncryostats() - 1); // last cryostat
    geo::GeometryCore::cryostat_id_iterator iLastCryoID(geom, LastID);

    BOOST_TEST_CHECKPOINT("Position-created iterator to last cryostat: "
      << std::string(LastID));
    geo::GeometryCore::cryostat_iterator iLastCryo(geom, LastID);
    CompareIteratorAndIteratorID(iLastCryo, iLastCryoID);

    // test increment to past-the-end
    geo::GeometryCore::cryostat_id_iterator iEndCryoID = iLastCryoID;
    ++iEndCryoID;

    geo::GeometryCore::cryostat_iterator iEndCryo = iLastCryo;
    ++iEndCryo;

    CompareIteratorAndIteratorID(iEndCryo, iEndCryoID);

  }

  //
  // end-constructed
  //
  {
    geo::GeometryCore::cryostat_id_iterator iCryoID = geom->end_cryostat_id();

    // construct from end position
    BOOST_TEST_CHECKPOINT
      ("End-created cryostat iterator: " << std::string(*iCryoID));
    geo::GeometryCore::cryostat_iterator iCryo
      (geom, geo::GeometryCore::cryostat_iterator::end_pos);
    CompareIteratorAndIteratorID(iCryo, iCryoID);

    // construct at end position by geometry
    geo::GeometryCore::cryostat_iterator iCryoGE = geom->end_cryostat();
    CompareIteratorAndIteratorID(iCryoGE, iCryoID);

  }

} // GeometryIteratorTestAlg::CryostatIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::TPCIDIteratorsTest() const {

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
    geo::GeometryCore::TPC_id_iterator iTPC;
    BOOST_TEST_CHECKPOINT
      ("Default created TPC ID iterator: " << std::string(*iTPC));

    BOOST_CHECK_EQUAL(iTPC->Cryostat, geo::CryostatID::getInvalidID());
    BOOST_CHECK_EQUAL(iTPC->TPC, geo::TPCID::getInvalidID());

    // check that the iterator tests false
    BOOST_CHECK(!iTPC);
    BOOST_CHECK(!(bool(iTPC)));

  }

  //
  // begin-constructed
  //
  {
    geo::GeometryCore::TPC_id_iterator iTPC(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created TPC ID iterator: " << std::string(*iTPC));

    BOOST_CHECK_EQUAL(iTPC->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iTPC->TPC,                geo::TPCID::TPCID_t(0));

    // check that the iterator tests true
    BOOST_CHECK(bool(iTPC));
    BOOST_CHECK(!!iTPC);

    // initialize to the beginning directly; this has probably ID's isValid true
    geo::TPCID BeginID;
    geom->GetBeginID(BeginID);
    geo::GeometryCore::TPC_id_iterator iTPCD(geom, BeginID);
    BOOST_CHECK_EQUAL(iTPCD->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iTPCD->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iTPCD, iTPC);

    // construct from explicit begin position
    geo::GeometryCore::TPC_id_iterator iTPCBC
      (geom, geo::GeometryCore::TPC_id_iterator::begin_pos);
    BOOST_CHECK_EQUAL(iTPCBC, iTPC);

    // construct at begin position by geometry
    geo::GeometryCore::TPC_id_iterator iTPCGB = geom->begin_TPC_id();
    BOOST_CHECK_EQUAL(iTPCGB, iTPC);

    // check access to ID
    BOOST_CHECK_EQUAL(*iTPC, BeginID);
    BOOST_CHECK_EQUAL(iTPC->Cryostat, BeginID.Cryostat);
    BOOST_CHECK_EQUAL(iTPC->TPC, BeginID.TPC);

    // check access to geometry element
    geo::TPCGeo const* pTPC = geom->TPCPtr(BeginID);
    BOOST_CHECK_EQUAL(iTPC.get(), pTPC);

    // test copy and postfix increment
    geo::GeometryCore::TPC_id_iterator iTPCI(iTPC++);

    const unsigned int nTPCsInC0 = geom->NTPC(geo::CryostatID(0));
    if (nTPCsInC0 > 1) {
      BOOST_CHECK_EQUAL(iTPCI->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iTPCI->TPC,                geo::TPCID::TPCID_t(0));
      BOOST_CHECK_EQUAL(iTPC->Cryostat,  geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iTPC->TPC,                 geo::TPCID::TPCID_t(1));
    }
    BOOST_CHECK_NE(iTPCI, iTPC);

    // test copy and prefix increment
    ++iTPCI;
    BOOST_CHECK_EQUAL(iTPCI, iTPC); // arguable if both are end-iterators by now

  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping cryostat
    geo::TPCID ID(0, 0);
    ID.TPC = geom->NTPC(ID) - 1; // last TPC of first cryostat

    geo::GeometryCore::TPC_id_iterator iTPC(geom, ID);

    // check that the iterator tests true
    BOOST_CHECK(bool(iTPC));
    BOOST_CHECK(!!iTPC);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iTPC, ID);
    BOOST_CHECK_EQUAL(iTPC->Cryostat, ID.Cryostat);
    BOOST_CHECK_EQUAL(iTPC->TPC, ID.TPC);
    BOOST_CHECK_EQUAL(iTPC.get(), geom->TPCPtr(ID));

    ++iTPC;
    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(iTPC->Cryostat, geo::CryostatID::CryostatID_t(ID.Cryostat + 1));
    BOOST_CHECK_EQUAL(iTPC->TPC,                geo::TPCID::TPCID_t(0));


    // test iterator to last TPC
    geo::TPCID LastID(geom->Ncryostats() - 1, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    geo::GeometryCore::TPC_id_iterator iLastTPC(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last TPC ID: "
      << std::string(*iLastTPC));

    // check that the iterator tests true
    BOOST_CHECK(bool(iLastTPC));
    BOOST_CHECK(!!iLastTPC);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iLastTPC, LastID);
    BOOST_CHECK_EQUAL(iLastTPC->Cryostat, LastID.Cryostat);
    BOOST_CHECK_EQUAL(iLastTPC->TPC, LastID.TPC);
    BOOST_CHECK_EQUAL(iLastTPC.get(), geom->TPCPtr(LastID));

    // test increment to past-the-end
    geo::GeometryCore::TPC_id_iterator iEndTPC = iLastTPC;
    ++iEndTPC;

    // check that the iterator tests false
    BOOST_CHECK(!bool(iEndTPC));
    BOOST_CHECK(!iEndTPC);

    BOOST_CHECK_EQUAL(iEndTPC->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iEndTPC->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iEndTPC, geom->end_TPC_id());
    BOOST_CHECK(!iEndTPC.get());

  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::GeometryCore::TPC_id_iterator iTPC
      (geom, geo::GeometryCore::TPC_id_iterator::end_pos);
    BOOST_TEST_CHECKPOINT
      ("End-created TPC ID iterator: " << std::string(*iTPC));

    BOOST_CHECK_EQUAL(iTPC->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iTPC->TPC,                geo::TPCID::TPCID_t(0));

    // check that the iterator tests false
    BOOST_CHECK(!bool(iTPC));
    BOOST_CHECK(!iTPC);

    // check access to geometry element (result of operator* is not defined)
    BOOST_CHECK(!(iTPC.get())); // should get nullptr

    // construct at end position by geometry
    geo::GeometryCore::TPC_id_iterator iTPCGE = geom->end_TPC_id();
    BOOST_CHECK_EQUAL(iTPCGE, iTPC);

    // initialize to the end directly; this has probably ID's isValid true
    geo::GeometryCore::TPC_id_iterator iTPC2
      (geom, geo::TPCID(geom->Ncryostats(), 0));
    BOOST_CHECK_EQUAL(iTPC2->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iTPC2->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iTPC2, iTPC);

  }

} // GeometryIteratorTestAlg::TPCIDIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::TPCIteratorsTest() const {

  /*
   * This test is extensively based on the assumption that the iterators should
   * behave like the corresponding ID iterators, including "corner cases".
   *
   * public interface (geometry_element_iterator<>)
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     geometry_element_iterator()
   *
   *     /// Constructor: points to begin
   *     geometry_element_iterator(geo::GeometryCore const* geom)
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
   *     /// Constructor: points to beginning
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, BeginPos_t const pos)
   *
   *     /// Constructor: points to end
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, EndPos_t const pos)
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
    geo::GeometryCore::TPC_id_iterator iTPCID;
    BOOST_TEST_CHECKPOINT
      ("Default created TPC iterator: " << std::string(*iTPCID));

    geo::GeometryCore::TPC_iterator iTPC;

    // direct comparison
    BOOST_CHECK_EQUAL(iTPC, iTPCID);
    BOOST_CHECK(!(iTPC != iTPCID));

    // ID comparison
    BOOST_CHECK_EQUAL(iTPC.ID(), *iTPCID);

    // check copy assignment
    geo::GeometryCore::TPC_iterator iTPC_copy(iTPC);
    geo::GeometryCore::TPC_id_iterator iTPCID_copy(iTPCID);

    // check comparisons too
    BOOST_CHECK_EQUAL(iTPC, iTPC_copy);
    BOOST_CHECK_EQUAL(iTPC_copy, iTPC);
    BOOST_CHECK(!(iTPC != iTPC_copy));
    BOOST_CHECK(!(iTPC_copy != iTPC));

    BOOST_CHECK_EQUAL(iTPC, iTPC_copy);

  }

  //
  // begin-constructed
  //
  {
    geo::TPCID BeginID;
    geom->GetBeginID(BeginID);

    geo::GeometryCore::TPC_id_iterator iTPCID(geom, BeginID);


    BOOST_TEST_CHECKPOINT
      ("Begin-created TPC iterator (" << std::string(BeginID) << ")");

    // initialize to the beginning directly
    geo::GeometryCore::TPC_iterator iTPCD(geom, BeginID);
    CompareIteratorAndIteratorID(iTPCD, iTPCID);

    // initialize to the beginning with implicit initialization
    geo::GeometryCore::TPC_iterator iTPC(geom);
    CompareIteratorAndIteratorID(iTPC, iTPCID);

    // construct from explicit begin position
    geo::GeometryCore::TPC_iterator iTPCBC
      (geom, geo::GeometryCore::TPC_iterator::begin_pos);
    CompareIteratorAndIteratorID(iTPCBC, iTPCID);

    // construct at begin position by geometry
    geo::GeometryCore::TPC_iterator iTPCGB = geom->begin_TPC();
    CompareIteratorAndIteratorID(iTPCGB, iTPCID);

  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last TPC
    geo::TPCID LastID(geom->Ncryostats() - 1, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    geo::GeometryCore::TPC_id_iterator iLastTPCID(geom, LastID);

    BOOST_TEST_CHECKPOINT("Position-created iterator to last TPC: "
      << std::string(LastID));
    geo::GeometryCore::TPC_iterator iLastTPC(geom, LastID);
    CompareIteratorAndIteratorID(iLastTPC, iLastTPCID);

    // test increment to past-the-end
    geo::GeometryCore::TPC_id_iterator iEndTPCID = iLastTPCID;
    ++iEndTPCID;

    geo::GeometryCore::TPC_iterator iEndTPC = iLastTPC;
    ++iEndTPC;

    CompareIteratorAndIteratorID(iEndTPC, iEndTPCID);

  }

  //
  // end-constructed
  //
  {
    geo::GeometryCore::TPC_id_iterator iTPCID = geom->end_TPC_id();

    // construct from end position
    BOOST_TEST_CHECKPOINT("End-created TPC iterator: " << std::string(*iTPCID));
    geo::GeometryCore::TPC_iterator iTPC
      (geom, geo::GeometryCore::TPC_iterator::end_pos);
    CompareIteratorAndIteratorID(iTPC, iTPCID);

    // construct at end position by geometry
    geo::GeometryCore::TPC_iterator iTPCGE = geom->end_TPC();
    CompareIteratorAndIteratorID(iTPCGE, iTPCID);

  }

} // GeometryIteratorTestAlg::TPCIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::PlaneIDIteratorsTest() const {

  /*
   * public interface (plane_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     plane_id_iterator_base()
   *
   *     /// Constructor: points to begin
   *     plane_id_iterator_base(geo::GeometryCore const* geom)
   *
   *     /// Constructor: points to the specified cryostat
   *     plane_id_iterator_base
   *       (geo::GeometryCore const* geom, GEOID const& start_from)
   *
   *     /// Constructor: points to begin
   *     plane_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t)
   *
   *     /// Constructor: points to end
   *     plane_id_iterator_base(geo::GeometryCore const* geom, EndPos_t)
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
    geo::GeometryCore::plane_id_iterator iPlane;
    BOOST_TEST_CHECKPOINT
      ("Default created plane ID iterator: " << std::string(*iPlane));

    BOOST_CHECK_EQUAL(iPlane->Cryostat, geo::CryostatID::getInvalidID());
    BOOST_CHECK_EQUAL(iPlane->TPC,           geo::TPCID::getInvalidID());
    BOOST_CHECK_EQUAL(iPlane->Plane,       geo::PlaneID::getInvalidID());

    // check that the iterator tests false
    BOOST_CHECK(!iPlane);
    BOOST_CHECK(!(bool(iPlane)));

  }

  //
  // begin-constructed
  //
  {
    geo::GeometryCore::plane_id_iterator iPlane(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created plane ID iterator: " << std::string(*iPlane));

    BOOST_CHECK_EQUAL(iPlane->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iPlane->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iPlane->Plane,          geo::PlaneID::PlaneID_t(0));

    // check that the iterator tests true
    BOOST_CHECK(bool(iPlane));
    BOOST_CHECK(!!iPlane);

    // initialize to the beginning directly; this has probably ID's isValid true
    geo::PlaneID BeginID;
    geom->GetBeginID(BeginID);
    geo::GeometryCore::plane_id_iterator iPlaneD(geom, BeginID);
    BOOST_CHECK_EQUAL(iPlaneD->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iPlaneD->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iPlaneD->Plane,          geo::PlaneID::PlaneID_t(0));
    BOOST_CHECK_EQUAL(iPlaneD, iPlane);

    // construct from explicit begin position
    geo::GeometryCore::plane_id_iterator iPlaneBC
      (geom, geo::GeometryCore::plane_id_iterator::begin_pos);
    BOOST_CHECK_EQUAL(iPlaneBC, iPlane);

    // construct at begin position by geometry
    geo::GeometryCore::plane_id_iterator iPlaneGB = geom->begin_plane_id();
    BOOST_CHECK_EQUAL(iPlaneGB, iPlane);

    // check access to ID
    BOOST_CHECK_EQUAL(*iPlane, BeginID);
    BOOST_CHECK_EQUAL(iPlane->Cryostat, BeginID.Cryostat);
    BOOST_CHECK_EQUAL(iPlane->TPC, BeginID.TPC);
    BOOST_CHECK_EQUAL(iPlane->Plane, BeginID.Plane);

    // check access to geometry element
    geo::PlaneGeo const* pPlane = geom->PlanePtr(BeginID);
    BOOST_CHECK_EQUAL(iPlane.get(), pPlane);

    // test copy and postfix increment
    geo::GeometryCore::plane_id_iterator iPlaneI(iPlane++);

    const unsigned int nPlanesInC0T0 = geom->Nplanes(geo::TPCID(0, 0));
    if (nPlanesInC0T0 > 1) {
      BOOST_CHECK_EQUAL(iPlaneI->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iPlaneI->TPC,                geo::TPCID::TPCID_t(0));
      BOOST_CHECK_EQUAL(iPlaneI->Plane,          geo::PlaneID::PlaneID_t(0));
      BOOST_CHECK_EQUAL(iPlane->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iPlane->TPC,                geo::TPCID::TPCID_t(0));
      BOOST_CHECK_EQUAL(iPlane->Plane,          geo::PlaneID::PlaneID_t(1));
    }
    BOOST_CHECK_NE(iPlaneI, iPlane);

    // test copy and prefix increment
    ++iPlaneI;
    BOOST_CHECK_EQUAL(iPlaneI, iPlane); // arguable if both are end-iterators by now

  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping TPC
    geo::PlaneID ID(0, 0, 0);
    ID.Plane = geom->Nplanes(ID) - 1; // last plane of first TPC

    geo::GeometryCore::plane_id_iterator iPlane(geom, ID);

    // check that the iterator tests true
    BOOST_CHECK(bool(iPlane));
    BOOST_CHECK(!!iPlane);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iPlane, ID);
    BOOST_CHECK_EQUAL(iPlane->Cryostat, ID.Cryostat);
    BOOST_CHECK_EQUAL(iPlane->TPC, ID.TPC);
    BOOST_CHECK_EQUAL(iPlane->Plane, ID.Plane);
    BOOST_CHECK_EQUAL(iPlane.get(), geom->PlanePtr(ID));

    // check that the pointed ID is as expected
    ++iPlane;
    if (ID.TPC + 1 < geom->NTPC(ID)) {
      BOOST_CHECK_EQUAL(iPlane->Cryostat, ID.Cryostat);
      BOOST_CHECK_EQUAL(iPlane->TPC, ID.TPC + 1);
      BOOST_CHECK_EQUAL(iPlane->Plane, geo::PlaneID::PlaneID_t(0));
    }
    else {
      BOOST_CHECK_EQUAL(iPlane->Cryostat, ID.Cryostat + 1);
      BOOST_CHECK_EQUAL(iPlane->TPC,       geo::TPCID::TPCID_t(0));
      BOOST_CHECK_EQUAL(iPlane->Plane, geo::PlaneID::PlaneID_t(0));
    }

    // test iterator to last plane
    geo::PlaneID LastID(geom->Ncryostats() - 1, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    geo::GeometryCore::plane_id_iterator iLastPlane(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last plane ID: "
      << std::string(*iLastPlane));

    // check that the iterator tests true
    BOOST_CHECK(bool(iLastPlane));
    BOOST_CHECK(!!iLastPlane);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iLastPlane, LastID);
    BOOST_CHECK_EQUAL(iLastPlane->Cryostat, LastID.Cryostat);
    BOOST_CHECK_EQUAL(iLastPlane->TPC, LastID.TPC);
    BOOST_CHECK_EQUAL(iLastPlane->Plane, LastID.Plane);
    BOOST_CHECK_EQUAL(iLastPlane.get(), geom->PlanePtr(LastID));

    // test increment to past-the-end
    geo::GeometryCore::plane_id_iterator iEndPlane = iLastPlane;
    ++iEndPlane;

    // check that the iterator tests false
    BOOST_CHECK(!bool(iEndPlane));
    BOOST_CHECK(!iEndPlane);

    BOOST_CHECK_EQUAL(iEndPlane->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iEndPlane->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iEndPlane->Plane,          geo::PlaneID::PlaneID_t(0));
    BOOST_CHECK_EQUAL(iEndPlane, geom->end_plane_id());
    BOOST_CHECK(!iEndPlane.get());

  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::GeometryCore::plane_id_iterator iPlane
      (geom, geo::GeometryCore::plane_id_iterator::end_pos);
    BOOST_TEST_CHECKPOINT
      ("End-created plane ID iterator: " << std::string(*iPlane));

    BOOST_CHECK_EQUAL(iPlane->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iPlane->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iPlane->Plane,          geo::PlaneID::PlaneID_t(0));

    // check that the iterator tests false
    BOOST_CHECK(!bool(iPlane));
    BOOST_CHECK(!iPlane);

    // check access to geometry element (result of operator* is not defined)
    BOOST_CHECK(!(iPlane.get())); // should get nullptr

    // construct at end position by geometry
    geo::GeometryCore::plane_id_iterator iPlaneGE = geom->end_plane_id();
    BOOST_CHECK_EQUAL(iPlaneGE, iPlane);

    // initialize to the end directly; this has probably ID's isValid true
    geo::GeometryCore::plane_id_iterator iPlane2
      (geom, geo::PlaneID(geom->Ncryostats(), 0, 0));
    BOOST_CHECK_EQUAL(iPlane2->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iPlane2->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iPlane2->Plane,          geo::PlaneID::PlaneID_t(0));
    BOOST_CHECK_EQUAL(iPlane2, iPlane);

  }
} // GeometryIteratorTestAlg::PlaneIDIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::PlaneIteratorsTest() const {

  /*
   * This test is extensively based on the assumption that the iterators should
   * behave like the corresponding ID iterators, including "corner cases".
   *
   * public interface (geometry_element_iterator<>)
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     geometry_element_iterator()
   *
   *     /// Constructor: points to begin
   *     geometry_element_iterator(geo::GeometryCore const* geom)
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
   *     /// Constructor: points to beginning
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, BeginPos_t const pos)
   *
   *     /// Constructor: points to end
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, EndPos_t const pos)
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
    geo::GeometryCore::plane_id_iterator iPlaneID;
    BOOST_TEST_CHECKPOINT
      ("Default created plane iterator: " << std::string(*iPlaneID));

    geo::GeometryCore::plane_iterator iPlane;

    // direct comparison
    BOOST_CHECK_EQUAL(iPlane, iPlaneID);
    BOOST_CHECK(!(iPlane != iPlaneID));

    // ID comparison
    BOOST_CHECK_EQUAL(iPlane.ID(), *iPlaneID);

    // check copy assignment
    geo::GeometryCore::plane_iterator iPlane_copy(iPlane);
    geo::GeometryCore::plane_id_iterator iPlaneID_copy(iPlaneID);

    // check comparisons too
    BOOST_CHECK_EQUAL(iPlane, iPlane_copy);
    BOOST_CHECK_EQUAL(iPlane_copy, iPlane);
    BOOST_CHECK(!(iPlane != iPlane_copy));
    BOOST_CHECK(!(iPlane_copy != iPlane));

    BOOST_CHECK_EQUAL(iPlane, iPlane_copy);

  }

  //
  // begin-constructed
  //
  {
    geo::PlaneID BeginID;
    geom->GetBeginID(BeginID);

    geo::GeometryCore::plane_id_iterator iPlaneID(geom, BeginID);


    BOOST_TEST_CHECKPOINT
      ("Begin-created plane iterator (" << std::string(BeginID) << ")");

    // initialize to the beginning directly
    geo::GeometryCore::plane_iterator iPlaneD(geom, BeginID);
    CompareIteratorAndIteratorID(iPlaneD, iPlaneID);

    // initialize to the beginning with implicit initialization
    geo::GeometryCore::plane_iterator iPlane(geom);
    CompareIteratorAndIteratorID(iPlane, iPlaneID);

    // construct from explicit begin position
    geo::GeometryCore::plane_iterator iPlaneBC
      (geom, geo::GeometryCore::plane_iterator::begin_pos);
    CompareIteratorAndIteratorID(iPlaneBC, iPlaneID);

    // construct at begin position by geometry
    geo::GeometryCore::plane_iterator iPlaneGB = geom->begin_plane();
    CompareIteratorAndIteratorID(iPlaneGB, iPlaneID);

  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last plane
    geo::PlaneID LastID(geom->Ncryostats() - 1, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    geo::GeometryCore::plane_id_iterator iLastPlaneID(geom, LastID);

    BOOST_TEST_CHECKPOINT("Position-created iterator to last plane: "
      << std::string(LastID));
    geo::GeometryCore::plane_iterator iLastPlane(geom, LastID);
    CompareIteratorAndIteratorID(iLastPlane, iLastPlaneID);

    // test increment to past-the-end
    geo::GeometryCore::plane_id_iterator iEndPlaneID = iLastPlaneID;
    ++iEndPlaneID;

    geo::GeometryCore::plane_iterator iEndPlane = iLastPlane;
    ++iEndPlane;

    CompareIteratorAndIteratorID(iEndPlane, iEndPlaneID);

  }

  //
  // end-constructed
  //
  {
    geo::GeometryCore::plane_id_iterator iPlaneID = geom->end_plane_id();

    // construct from end position
    BOOST_TEST_CHECKPOINT
      ("End-created plane iterator: " << std::string(*iPlaneID));
    geo::GeometryCore::plane_iterator iPlane
      (geom, geo::GeometryCore::plane_iterator::end_pos);
    CompareIteratorAndIteratorID(iPlane, iPlaneID);

    // construct at end position by geometry
    geo::GeometryCore::plane_iterator iPlaneGE = geom->end_plane();
    CompareIteratorAndIteratorID(iPlaneGE, iPlaneID);

  }

} // GeometryIteratorTestAlg::PlaneIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::WireIDIteratorsTest() const {

  /*
   * public interface (wire_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     wire_id_iterator_base()
   *
   *     /// Constructor: points to begin
   *     wire_id_iterator_base(geo::GeometryCore const* geom)
   *
   *     /// Constructor: points to the specified cryostat
   *     wire_id_iterator_base
   *       (geo::GeometryCore const* geom, GEOID const& start_from)
   *
   *     /// Constructor: points to begin
   *     wire_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t)
   *
   *     /// Constructor: points to end
   *     wire_id_iterator_base(geo::GeometryCore const* geom, EndPos_t)
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
    geo::GeometryCore::wire_id_iterator iWire;
    BOOST_TEST_CHECKPOINT
      ("Default created wire ID iterator: " << std::string(*iWire));

    BOOST_CHECK_EQUAL(iWire->Cryostat, geo::CryostatID::getInvalidID());
    BOOST_CHECK_EQUAL(iWire->TPC,           geo::TPCID::getInvalidID());
    BOOST_CHECK_EQUAL(iWire->Plane,       geo::PlaneID::getInvalidID());
    BOOST_CHECK_EQUAL(iWire->Wire,         geo::WireID::getInvalidID());

    // check that the iterator tests false
    BOOST_CHECK(!iWire);
    BOOST_CHECK(!(bool(iWire)));

  }

  //
  // begin-constructed
  //
  {
    geo::GeometryCore::wire_id_iterator iWire(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created wire ID iterator: " << std::string(*iWire));

    BOOST_CHECK_EQUAL(iWire->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iWire->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iWire->Plane,          geo::PlaneID::PlaneID_t(0));
    BOOST_CHECK_EQUAL(iWire->Wire,             geo::WireID::WireID_t(0));

    // check that the iterator tests true
    BOOST_CHECK(bool(iWire));
    BOOST_CHECK(!!iWire);

    // initialize to the beginning directly; this has probably ID's isValid true
    geo::WireID BeginID;
    geom->GetBeginID(BeginID);
    geo::GeometryCore::wire_id_iterator iWireD(geom, BeginID);
    BOOST_CHECK_EQUAL(iWireD->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iWireD->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iWireD->Plane,          geo::PlaneID::PlaneID_t(0));
    BOOST_CHECK_EQUAL(iWireD->Wire,             geo::WireID::WireID_t(0));
    BOOST_CHECK_EQUAL(iWireD, iWire);

    // construct from explicit begin position
    geo::GeometryCore::wire_id_iterator iWireBC
      (geom, geo::GeometryCore::wire_id_iterator::begin_pos);
    BOOST_CHECK_EQUAL(iWireBC, iWire);

    // construct at begin position by geometry
    geo::GeometryCore::wire_id_iterator iWireGB = geom->begin_wire_id();
    BOOST_CHECK_EQUAL(iWireGB, iWire);

    // check access to ID
    BOOST_CHECK_EQUAL(*iWire, BeginID);
    BOOST_CHECK_EQUAL(iWire->Cryostat, BeginID.Cryostat);
    BOOST_CHECK_EQUAL(iWire->TPC, BeginID.TPC);
    BOOST_CHECK_EQUAL(iWire->Plane, BeginID.Plane);
    BOOST_CHECK_EQUAL(iWire->Wire, BeginID.Wire);

    // check access to geometry element
    geo::WireGeo const* pWire = geom->WirePtr(BeginID);
    BOOST_CHECK_EQUAL(iWire.get(), pWire);

    // test copy and postfix increment
    geo::GeometryCore::wire_id_iterator iWireI(iWire++);

    const unsigned int nWiresInC0T0P0 = geom->Nwires(geo::PlaneID(0, 0, 0));
    if (nWiresInC0T0P0 > 1) {
      BOOST_CHECK_EQUAL(iWireI->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iWireI->TPC,                geo::TPCID::TPCID_t(0));
      BOOST_CHECK_EQUAL(iWireI->Plane,          geo::PlaneID::PlaneID_t(0));
      BOOST_CHECK_EQUAL(iWireI->Wire,             geo::WireID::WireID_t(0));
      BOOST_CHECK_EQUAL(iWire->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iWire->TPC,                geo::TPCID::TPCID_t(0));
      BOOST_CHECK_EQUAL(iWire->Plane,          geo::PlaneID::PlaneID_t(0));
      BOOST_CHECK_EQUAL(iWire->Wire,             geo::WireID::WireID_t(1));
    }
    BOOST_CHECK_NE(iWireI, iWire);

    // test copy and prefix increment
    ++iWireI;
    BOOST_CHECK_EQUAL(iWireI, iWire); // arguable if both are end-iterators by now

  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping plane
    geo::WireID ID(0, 0, 0, 0);
    ID.Wire = geom->Nwires(ID) - 1; // last wire of first plane

    geo::GeometryCore::wire_id_iterator iWire(geom, ID);

    // check that the iterator tests true
    BOOST_CHECK(bool(iWire));
    BOOST_CHECK(!!iWire);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iWire, ID);
    BOOST_CHECK_EQUAL(iWire->Cryostat, ID.Cryostat);
    BOOST_CHECK_EQUAL(iWire->TPC, ID.TPC);
    BOOST_CHECK_EQUAL(iWire->Plane, ID.Plane);
    BOOST_CHECK_EQUAL(iWire->Wire, ID.Wire);
    BOOST_CHECK_EQUAL(iWire.get(), geom->WirePtr(ID));

    ++iWire;
    // check that the pointed ID is as expected
    if (ID.Plane + 1 < geom->Nplanes(ID)) {
      BOOST_CHECK_EQUAL(iWire->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iWire->TPC,                geo::TPCID::TPCID_t(0));
      BOOST_CHECK_EQUAL(iWire->Plane,                                  ID.Plane + 1);
      BOOST_CHECK_EQUAL(iWire->Wire,             geo::WireID::WireID_t(0));
    } else if (ID.TPC + 1 < geom->NTPC(ID)) {
      BOOST_CHECK_EQUAL(iWire->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iWire->TPC,                                    ID.TPC + 1);
      BOOST_CHECK_EQUAL(iWire->Plane,          geo::PlaneID::PlaneID_t(0));
      BOOST_CHECK_EQUAL(iWire->Wire,             geo::WireID::WireID_t(0));
    } else {
      BOOST_CHECK_EQUAL(iWire->Cryostat,                               ID.Cryostat + 1);
      BOOST_CHECK_EQUAL(iWire->TPC,                geo::TPCID::TPCID_t(0));
      BOOST_CHECK_EQUAL(iWire->Plane,          geo::PlaneID::PlaneID_t(0));
      BOOST_CHECK_EQUAL(iWire->Wire,             geo::WireID::WireID_t(0));
         }

    // test iterator to last wire
    geo::WireID LastID(geom->Ncryostats() - 1, 0, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    LastID.Wire = geom->Nwires(LastID) - 1; // last wire of last plane
    geo::GeometryCore::wire_id_iterator iLastWire(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last wire ID: "
      << std::string(*iLastWire));

    // check that the iterator tests true
    BOOST_CHECK(bool(iLastWire));
    BOOST_CHECK(!!iLastWire);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iLastWire, LastID);
    BOOST_CHECK_EQUAL(iLastWire->Cryostat, LastID.Cryostat);
    BOOST_CHECK_EQUAL(iLastWire->TPC, LastID.TPC);
    BOOST_CHECK_EQUAL(iLastWire->Plane, LastID.Plane);
    BOOST_CHECK_EQUAL(iLastWire->Wire, LastID.Wire);
    BOOST_CHECK_EQUAL(iLastWire.get(), geom->WirePtr(LastID));

    // test increment to past-the-end
    geo::GeometryCore::wire_id_iterator iEndWire = iLastWire;
    ++iEndWire;

    // check that the iterator tests false
    BOOST_CHECK(!bool(iEndWire));
    BOOST_CHECK(!iEndWire);

    BOOST_CHECK_EQUAL(iEndWire->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iEndWire->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iEndWire->Plane,          geo::PlaneID::PlaneID_t(0));
    BOOST_CHECK_EQUAL(iEndWire->Wire,             geo::WireID::WireID_t(0));
    BOOST_CHECK_EQUAL(iEndWire, geom->end_wire_id());
    BOOST_CHECK(!iEndWire.get());

  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::GeometryCore::wire_id_iterator iWire
      (geom, geo::GeometryCore::wire_id_iterator::end_pos);
    BOOST_TEST_CHECKPOINT
      ("End-created end ID iterator: " << std::string(*iWire));

    BOOST_CHECK_EQUAL(iWire->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iWire->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iWire->Plane,          geo::PlaneID::PlaneID_t(0));
    BOOST_CHECK_EQUAL(iWire->Wire,             geo::WireID::WireID_t(0));

    // check that the iterator tests false
    BOOST_CHECK(!bool(iWire));
    BOOST_CHECK(!iWire);

    // check access to geometry element (result of operator* is not defined)
    BOOST_CHECK(!(iWire.get())); // should get nullptr

    // construct at end position by geometry
    geo::GeometryCore::wire_id_iterator iWireGE = geom->end_wire_id();
    BOOST_CHECK_EQUAL(iWireGE, iWire);

    // initialize to the end directly; this has probably ID's isValid true
    geo::GeometryCore::wire_id_iterator iWire2
      (geom, geo::WireID(geom->Ncryostats(), 0, 0, 0));
    BOOST_CHECK_EQUAL(iWire2->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iWire2->TPC,                geo::TPCID::TPCID_t(0));
    BOOST_CHECK_EQUAL(iWire2->Plane,          geo::PlaneID::PlaneID_t(0));
    BOOST_CHECK_EQUAL(iWire2->Wire,             geo::WireID::WireID_t(0));
    BOOST_CHECK_EQUAL(iWire2, iWire);

  }

} // GeometryIteratorTestAlg::WireIDIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::WireIteratorsTest() const {

  /*
   * This test is extensively based on the assumption that the iterators should
   * behave like the corresponding ID iterators, including "corner cases".
   *
   * public interface (geometry_element_iterator<>)
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     geometry_element_iterator()
   *
   *     /// Constructor: points to begin
   *     geometry_element_iterator(geo::GeometryCore const* geom)
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
   *     /// Constructor: points to beginning
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, BeginPos_t const pos)
   *
   *     /// Constructor: points to end
   *     geometry_element_iterator
   *       (geo::GeometryCore const* geom, EndPos_t const pos)
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
    geo::GeometryCore::wire_id_iterator iWireID;
    BOOST_TEST_CHECKPOINT
      ("Default created wire iterator: " << std::string(*iWireID));

    geo::GeometryCore::wire_iterator iWire;

    // direct comparison
    BOOST_CHECK_EQUAL(iWire, iWireID);
    BOOST_CHECK(!(iWire != iWireID));

    // ID comparison
    BOOST_CHECK_EQUAL(iWire.ID(), *iWireID);

    // check copy assignment
    geo::GeometryCore::wire_iterator iWire_copy(iWire);
    geo::GeometryCore::wire_id_iterator iWireID_copy(iWireID);

    // check comparisons too
    BOOST_CHECK_EQUAL(iWire, iWire_copy);
    BOOST_CHECK_EQUAL(iWire_copy, iWire);
    BOOST_CHECK(!(iWire != iWire_copy));
    BOOST_CHECK(!(iWire_copy != iWire));

    BOOST_CHECK_EQUAL(iWire, iWire_copy);

  }

  //
  // begin-constructed
  //
  {
    geo::WireID BeginID;
    geom->GetBeginID(BeginID);

    geo::GeometryCore::wire_id_iterator iWireID(geom, BeginID);


    BOOST_TEST_CHECKPOINT
      ("Begin-created wire iterator (" << std::string(BeginID) << ")");

    // initialize to the beginning directly
    geo::GeometryCore::wire_iterator iWireD(geom, BeginID);
    CompareIteratorAndIteratorID(iWireD, iWireID);

    // initialize to the beginning with implicit initialization
    geo::GeometryCore::wire_iterator iWire(geom);
    CompareIteratorAndIteratorID(iWire, iWireID);

    // construct from explicit begin position
    geo::GeometryCore::wire_iterator iWireBC
      (geom, geo::GeometryCore::wire_iterator::begin_pos);
    CompareIteratorAndIteratorID(iWireBC, iWireID);

    // construct at begin position by geometry
    geo::GeometryCore::wire_iterator iWireGB = geom->begin_wire();
    CompareIteratorAndIteratorID(iWireGB, iWireID);

  }

  //
  // constructed from starting point
  //
  {
    // test iterator to last wire
    geo::WireID LastID(geom->Ncryostats() - 1, 0, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    LastID.Wire = geom->Nwires(LastID) - 1; // last wire of last plane
    geo::GeometryCore::wire_id_iterator iLastWireID(geom, LastID);

    BOOST_TEST_CHECKPOINT
      ("Position-created iterator to last wire: " << std::string(LastID));
    geo::GeometryCore::wire_iterator iLastWire(geom, LastID);
    CompareIteratorAndIteratorID(iLastWire, iLastWireID);

    // test increment to past-the-end
    geo::GeometryCore::wire_id_iterator iEndWireID = iLastWireID;
    ++iEndWireID;

    geo::GeometryCore::wire_iterator iEndWire = iLastWire;
    ++iEndWire;

    CompareIteratorAndIteratorID(iEndWire, iEndWireID);

  }

  //
  // end-constructed
  //
  {
    geo::GeometryCore::wire_id_iterator iWireID = geom->end_wire_id();

    // construct from end position
    BOOST_TEST_CHECKPOINT
      ("End-created wire iterator: " << std::string(*iWireID));
    geo::GeometryCore::wire_iterator iWire
      (geom, geo::GeometryCore::wire_iterator::end_pos);
    CompareIteratorAndIteratorID(iWire, iWireID);

    // construct at end position by geometry
    geo::GeometryCore::wire_iterator iWireGE = geom->end_wire();
    CompareIteratorAndIteratorID(iWireGE, iWireID);

  }

} // GeometryIteratorTestAlg::WireIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::TPCsetIDIteratorsTest() const {

  /*
   * public interface (TPCset_id_iterator_base):
   *
   *
   *   /// Default constructor; effect not defined: assign to it before using!
   *   TPCset_id_iterator_base()
   *
   *   /// Constructor: points to begin.
   *   TPCset_id_iterator_base(geo::GeometryCore const* geom)
   *
   *   /// Constructor: points to the specified cryostat.
   *   TPCset_id_iterator_base
   *     (geo::GeometryCore const* geom, GeoID_t const& start_from)
   *
   *   /// Constructor: points to begin.
   *   TPCset_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t const)
   *
   *   /// Constructor: points to end.
   *   TPCset_id_iterator_base(geo::GeometryCore const* geom, EndPos_t)
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
    BOOST_TEST_CHECKPOINT
      ("Default created TPC set ID iterator: " << std::string(*iTPCset));

    BOOST_CHECK_EQUAL(iTPCset->Cryostat, geo::CryostatID::getInvalidID());
    BOOST_CHECK_EQUAL(iTPCset->TPCset, readout::TPCsetID::getInvalidID());

    // check that the iterator tests false
    BOOST_CHECK(!iTPCset);
    BOOST_CHECK(!(bool(iTPCset)));

  }

  //
  // begin-constructed
  //
  {
    geo::TPCset_id_iterator iTPCset(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created TPC set ID iterator: " << std::string(*iTPCset));

    BOOST_CHECK_EQUAL(iTPCset->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iTPCset->TPCset,   readout::TPCsetID::TPCsetID_t(0));

    // check that the iterator tests true
    BOOST_CHECK(bool(iTPCset));
    BOOST_CHECK(!!iTPCset);

    // initialize to the beginning directly; this has probably ID's isValid true
    readout::TPCsetID BeginID;
    geom->GetBeginID(BeginID);
    geo::TPCset_id_iterator iTPCsetD(geom, BeginID);
    BOOST_CHECK_EQUAL(iTPCsetD->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iTPCsetD->TPCset,   readout::TPCsetID::TPCsetID_t(0));
    BOOST_CHECK_EQUAL(iTPCsetD, iTPCset);

    // construct from explicit begin position
    geo::TPCset_id_iterator iTPCsetBC(geom, geo::iterators::begin_pos);
    BOOST_CHECK_EQUAL(iTPCsetBC, iTPCset);

    // construct at begin position by geometry
    geo::TPCset_id_iterator iTPCsetGB = geom->begin_TPCset_id();
    BOOST_CHECK_EQUAL(iTPCsetGB, iTPCset);

    // check access to ID
    BOOST_CHECK_EQUAL(*iTPCset, BeginID);
    BOOST_CHECK_EQUAL(iTPCset->Cryostat, BeginID.Cryostat);
    BOOST_CHECK_EQUAL(iTPCset->TPCset, BeginID.TPCset);

    // test copy and postfix increment
    geo::TPCset_id_iterator iTPCsetI(iTPCset++);

    const unsigned int nTPCsetsInC0 = geom->NTPCsets(geo::CryostatID(0));
    if (nTPCsetsInC0 > 1) {
      BOOST_CHECK_EQUAL(iTPCsetI->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iTPCsetI->TPCset,   readout::TPCsetID::TPCsetID_t(0));
      BOOST_CHECK_EQUAL(iTPCset->Cryostat,  geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iTPCset->TPCset,    readout::TPCsetID::TPCsetID_t(1));
    }
    BOOST_CHECK_NE(iTPCsetI, iTPCset);

    // test copy and prefix increment
    ++iTPCsetI;
    BOOST_CHECK_EQUAL(iTPCsetI, iTPCset); // arguable if both are end-iterators by now

  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping cryostat
    readout::TPCsetID ID(0, 0);
    ID.TPCset = geom->NTPCsets(ID) - 1; // last TPC set of first cryostat

    geo::TPCset_id_iterator iTPCset(geom, ID);

    // check that the iterator tests true
    BOOST_CHECK(bool(iTPCset));
    BOOST_CHECK(!!iTPCset);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iTPCset, ID);
    BOOST_CHECK_EQUAL(iTPCset->Cryostat, ID.Cryostat);
    BOOST_CHECK_EQUAL(iTPCset->TPCset, ID.TPCset);

    ++iTPCset;
    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL
      (iTPCset->Cryostat, geo::CryostatID::CryostatID_t(ID.Cryostat + 1));
    BOOST_CHECK_EQUAL(iTPCset->TPCset, readout::TPCsetID::TPCsetID_t(0));


    // test iterator to last TPC
    readout::TPCsetID LastID(geom->Ncryostats() - 1, 0);
    LastID.TPCset = geom->NTPCsets(LastID) - 1; // last TPC set of last cryostat
    geo::TPCset_id_iterator iLastTPCset(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last TPC set ID: "
      << std::string(*iLastTPCset));

    // check that the iterator tests true
    BOOST_CHECK(bool(iLastTPCset));
    BOOST_CHECK(!!iLastTPCset);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iLastTPCset, LastID);
    BOOST_CHECK_EQUAL(iLastTPCset->Cryostat, LastID.Cryostat);
    BOOST_CHECK_EQUAL(iLastTPCset->TPCset, LastID.TPCset);

    // test increment to past-the-end
    geo::TPCset_id_iterator iEndTPCset = iLastTPCset;
    ++iEndTPCset;

    // check that the iterator tests false
    BOOST_CHECK(!bool(iEndTPCset));
    BOOST_CHECK(!iEndTPCset);

    BOOST_CHECK_EQUAL
      (iEndTPCset->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iEndTPCset->TPCset, readout::TPCsetID::TPCsetID_t(0));
    BOOST_CHECK_EQUAL(iEndTPCset, geom->end_TPCset_id());

  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::TPCset_id_iterator iTPCset(geom, geo::iterators::end_pos);
    BOOST_TEST_CHECKPOINT
      ("End-created TPC set ID iterator: " << std::string(*iTPCset));

    BOOST_CHECK_EQUAL
      (iTPCset->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iTPCset->TPCset, readout::TPCsetID::TPCsetID_t(0));

    // check that the iterator tests false
    BOOST_CHECK(!bool(iTPCset));
    BOOST_CHECK(!iTPCset);

    // construct at end position by geometry
    geo::TPCset_id_iterator iTPCsetGE = geom->end_TPCset_id();
    BOOST_CHECK_EQUAL(iTPCsetGE, iTPCset);

    // initialize to the end directly; this has probably ID's isValid true
    geo::TPCset_id_iterator iTPCset2
      (geom, readout::TPCsetID(geom->Ncryostats(), 0));
    BOOST_CHECK_EQUAL
      (iTPCset2->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iTPCset2->TPCset, readout::TPCsetID::TPCsetID_t(0));
    BOOST_CHECK_EQUAL(iTPCset2, iTPCset);

  }

} // GeometryIteratorTestAlg::TPCsetIDIteratorsTest()



//-----------------------------------------------------------------------------
void geo::GeometryIteratorTestAlg::ROPIDIteratorsTest() const {

  /*
   * public interface (ROP_id_iterator_base):
   *
   *     /// Default constructor; effect not defined: assign to it before using!
   *     ROP_id_iterator_base()
   *
   *     /// Constructor: points to begin.
   *     ROP_id_iterator_base(geo::GeometryCore const* geom)
   *
   *     /// Constructor: points to the specified readout plane.
   *     ROP_id_iterator_base
   *       (geo::GeometryCore const* geom, GeoID_t const& start_from)
   *
   *     /// Constructor: points to begin.
   *     ROP_id_iterator_base(geo::GeometryCore const* geom, BeginPos_t const)
   *
   *     /// Constructor: points to end.
   *     ROP_id_iterator_base(geo::GeometryCore const* geom, EndPos_t)
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
    BOOST_TEST_CHECKPOINT
      ("Default created readout plane ID iterator: " << std::string(*iROP));

    BOOST_CHECK_EQUAL(iROP->Cryostat, geo::CryostatID::getInvalidID());
    BOOST_CHECK_EQUAL(iROP->TPCset,   readout::TPCsetID::getInvalidID());
    BOOST_CHECK_EQUAL(iROP->ROP,      readout::ROPID::getInvalidID());

    // check that the iterator tests false
    BOOST_CHECK(!iROP);
    BOOST_CHECK(!(bool(iROP)));

  }

  //
  // begin-constructed
  //
  {
    geo::ROP_id_iterator iROP(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created readout plane ID iterator: " << std::string(*iROP));

    BOOST_CHECK_EQUAL(iROP->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iROP->TPCset,   readout::TPCsetID::TPCsetID_t(0));
    BOOST_CHECK_EQUAL(iROP->ROP,      readout::ROPID::ROPID_t(0));

    // check that the iterator tests true
    BOOST_CHECK(bool(iROP));
    BOOST_CHECK(!!iROP);

    // initialize to the beginning directly; this has probably ID's isValid true
    readout::ROPID BeginID;
    geom->GetBeginID(BeginID);
    geo::ROP_id_iterator iROPD(geom, BeginID);
    BOOST_CHECK_EQUAL(iROPD->Cryostat, geo::CryostatID::CryostatID_t(0));
    BOOST_CHECK_EQUAL(iROPD->TPCset,   readout::TPCsetID::TPCsetID_t(0));
    BOOST_CHECK_EQUAL(iROPD->ROP,      readout::ROPID::ROPID_t(0));
    BOOST_CHECK_EQUAL(iROPD, iROP);

    // construct from explicit begin position
    geo::ROP_id_iterator iROPBC(geom, geo::iterators::begin_pos);
    BOOST_CHECK_EQUAL(iROPBC, iROP);

    // construct at begin position by geometry
    geo::ROP_id_iterator iROPGB = geom->begin_ROP_id();
    BOOST_CHECK_EQUAL(iROPGB, iROP);

    // check access to ID
    BOOST_CHECK_EQUAL(*iROP, BeginID);
    BOOST_CHECK_EQUAL(iROP->Cryostat, BeginID.Cryostat);
    BOOST_CHECK_EQUAL(iROP->TPCset, BeginID.TPCset);
    BOOST_CHECK_EQUAL(iROP->ROP, BeginID.ROP);

    // test copy and postfix increment
    geo::ROP_id_iterator iROPI(iROP++);

    const unsigned int nReadoutPlanesInC0S0
      = geom->NROPs(readout::TPCsetID(0, 0));
    if (nReadoutPlanesInC0S0 > 1) {
      BOOST_CHECK_EQUAL(iROPI->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iROPI->TPCset,   readout::TPCsetID::TPCsetID_t(0));
      BOOST_CHECK_EQUAL(iROPI->ROP,            readout::ROPID::ROPID_t(0));
      BOOST_CHECK_EQUAL(iROP->Cryostat, geo::CryostatID::CryostatID_t(0));
      BOOST_CHECK_EQUAL(iROP->TPCset,   readout::TPCsetID::TPCsetID_t(0));
      BOOST_CHECK_EQUAL(iROP->ROP,            readout::ROPID::ROPID_t(1));
    }
    BOOST_CHECK_NE(iROPI, iROP);

    // test copy and prefix increment
    ++iROPI;
    BOOST_CHECK_EQUAL(iROPI, iROP); // arguable if both are end-iterators by now

  }

  //
  // constructed from starting point
  //
  {
    // test increment flipping TPC
    readout::ROPID ID(0, 0, 0);
    ID.ROP = geom->NROPs(ID) - 1; // last plane of first TPC set

    geo::ROP_id_iterator iROP(geom, ID);

    // check that the iterator tests true
    BOOST_CHECK(bool(iROP));
    BOOST_CHECK(!!iROP);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iROP, ID);
    BOOST_CHECK_EQUAL(iROP->Cryostat, ID.Cryostat);
    BOOST_CHECK_EQUAL(iROP->TPCset, ID.TPCset);
    BOOST_CHECK_EQUAL(iROP->ROP, ID.ROP);

    // check that the pointed ID is as expected
    ++iROP;
    if (ID.TPCset + 1 < (int) geom->NTPCsets(ID)) {
      BOOST_CHECK_EQUAL(iROP->Cryostat, ID.Cryostat);
      BOOST_CHECK_EQUAL(iROP->TPCset,   ID.TPCset + 1);
      BOOST_CHECK_EQUAL(iROP->ROP,      readout::ROPID::ROPID_t(0));
    }
    else {
      BOOST_CHECK_EQUAL(iROP->Cryostat, ID.Cryostat + 1);
      BOOST_CHECK_EQUAL(iROP->TPCset,   readout::TPCsetID::TPCsetID_t(0));
      BOOST_CHECK_EQUAL(iROP->ROP,      readout::ROPID::ROPID_t(0));
    }

    // test iterator to last plane
    readout::ROPID LastID(geom->Ncryostats() - 1, 0, 0);
    LastID.TPCset = geom->NTPCsets(LastID) - 1; // last TPC set of last cryostat
    LastID.ROP = geom->NROPs(LastID) - 1; // last readout plane of last TPC set
    geo::ROP_id_iterator iLastROP(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last readout plane ID: "
      << std::string(*iLastROP));

    // check that the iterator tests true
    BOOST_CHECK(bool(iLastROP));
    BOOST_CHECK(!!iLastROP);

    // check that the pointed ID is as expected
    BOOST_CHECK_EQUAL(*iLastROP, LastID);
    BOOST_CHECK_EQUAL(iLastROP->Cryostat, LastID.Cryostat);
    BOOST_CHECK_EQUAL(iLastROP->TPCset,   LastID.TPCset);
    BOOST_CHECK_EQUAL(iLastROP->ROP,      LastID.ROP);

    // test increment to past-the-end
    geo::ROP_id_iterator iEndROP = iLastROP;
    ++iEndROP;

    // check that the iterator tests false
    BOOST_CHECK(!bool(iEndROP));
    BOOST_CHECK(!iEndROP);

    BOOST_CHECK_EQUAL
      (iEndROP->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iEndROP->TPCset,   readout::TPCsetID::TPCsetID_t(0));
    BOOST_CHECK_EQUAL(iEndROP->ROP,      readout::ROPID::ROPID_t(0));
    BOOST_CHECK_EQUAL(iEndROP, geom->end_ROP_id());

  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::ROP_id_iterator iROP(geom, geo::iterators::end_pos);
    BOOST_TEST_CHECKPOINT
      ("End-created readout plane ID iterator: " << std::string(*iROP));

    BOOST_CHECK_EQUAL
      (iROP->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iROP->TPCset,   readout::TPCsetID::TPCsetID_t(0));
    BOOST_CHECK_EQUAL(iROP->ROP,            readout::ROPID::ROPID_t(0));

    // check that the iterator tests false
    BOOST_CHECK(!bool(iROP));
    BOOST_CHECK(!iROP);

    // construct at end position by geometry
    geo::ROP_id_iterator iROPGE = geom->end_ROP_id();
    BOOST_CHECK_EQUAL(iROPGE, iROP);

    // initialize to the end directly; this has probably ID's isValid true
    geo::ROP_id_iterator iROP2
      (geom, readout::ROPID(geom->Ncryostats(), 0, 0));
    BOOST_CHECK_EQUAL
      (iROP->Cryostat, geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_CHECK_EQUAL(iROP->TPCset,   readout::TPCsetID::TPCsetID_t(0));
    BOOST_CHECK_EQUAL(iROP->ROP,            readout::ROPID::ROPID_t(0));
    BOOST_CHECK_EQUAL(iROP2, iROP);

  }
} // GeometryIteratorTestAlg::ROPIDIteratorsTest()
