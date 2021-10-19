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
      BOOST_TEST(iter == id_iter);
      BOOST_TEST(! (iter != id_iter));

      // ID comparison
      BOOST_TEST(iter.ID() == *id_iter);

      auto pGeoElement = id_iter.get();

      // dereference
      BOOST_TEST(iter.get() == pGeoElement);
      BOOST_TEST(iter.operator->() == pGeoElement);

      if (pGeoElement) BOOST_TEST(&*iter == pGeoElement);
      else             BOOST_CHECK_THROW(*iter, cet::exception);

      // boolean conversions
      BOOST_TEST(bool(iter) == bool(id_iter));

      // check copy assignment
      ITER iter_copy(iter);
      ITERID id_iter_copy(id_iter);

      // check comparisons too
      BOOST_TEST(iter == iter_copy);
      BOOST_TEST(iter_copy == iter);
      BOOST_TEST(!(iter != iter_copy));
      BOOST_TEST(!(iter_copy != iter));

      // check increment operator
      BOOST_TEST(iter++ == id_iter++);
      BOOST_TEST(++iter_copy == ++id_iter_copy);

      BOOST_TEST(iter == iter_copy);

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

    BOOST_TEST(iCryo->Cryostat == geo::CryostatID::getInvalidID());

    // check that the iterator tests false
    BOOST_TEST(!iCryo);
    BOOST_TEST(!(bool(iCryo)));

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

    BOOST_TEST(iCryo->Cryostat == geo::CryostatID::CryostatID_t(0));

    // check that the iterator tests true
    BOOST_TEST(bool(iCryo));
    BOOST_TEST(!!iCryo);

    // initialize to the beginning directly; this has probably ID's isValid true
    geo::CryostatID BeginID;
    geom->GetBeginID(BeginID);
    geo::GeometryCore::cryostat_id_iterator iCryoD(geom, BeginID);
    BOOST_TEST(iCryoD->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iCryoD == iCryo);

    // construct from explicit begin position
    geo::GeometryCore::cryostat_id_iterator iCryoBC
      (geom, geo::GeometryCore::cryostat_id_iterator::begin_pos);
    BOOST_TEST(iCryoBC == iCryo);

    // construct at begin position by geometry
    geo::GeometryCore::cryostat_id_iterator iCryoGB = geom->begin_cryostat_id();
    BOOST_TEST(iCryoGB == iCryo);

    // check access to ID
    BOOST_TEST(*iCryo == BeginID);
    BOOST_TEST(iCryo->Cryostat == BeginID.Cryostat);

    // check access to geometry element
    geo::CryostatGeo const* pCryo = geom->CryostatPtr(BeginID);
    BOOST_TEST(iCryo.get() == pCryo);

    // test copy and postfix increment
    geo::GeometryCore::cryostat_id_iterator iCryoI(iCryo++);

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

    geo::GeometryCore::cryostat_id_iterator iLastCryo(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last cryostat ID: "
      << std::string(*iLastCryo));

    // check that the iterator tests true
    BOOST_TEST(bool(iLastCryo));
    BOOST_TEST(!!iLastCryo);

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastCryo == LastID);
    BOOST_TEST(iLastCryo->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastCryo.get() == geom->CryostatPtr(LastID));

    // test increment to past-the-end
    geo::GeometryCore::cryostat_id_iterator iEndCryo = iLastCryo;
    ++iEndCryo;

    // check that the iterator tests false
    BOOST_TEST(!bool(iEndCryo));
    BOOST_TEST(!iEndCryo);

    BOOST_TEST(iEndCryo->Cryostat == geom->Ncryostats());
    BOOST_TEST(iEndCryo == geom->end_cryostat_id());
    BOOST_TEST(!iEndCryo.get());

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

    BOOST_TEST(iCryo->Cryostat == geom->Ncryostats());

    // check that the iterator tests false
    BOOST_TEST(!bool(iCryo));
    BOOST_TEST(!iCryo);

    // check access to geometry element (result of operator* is not defined)
    BOOST_TEST(!(iCryo.get())); // should get nullptr

    // construct at end position by geometry
    geo::GeometryCore::cryostat_id_iterator iCryoGE = geom->end_cryostat_id();
    BOOST_TEST(iCryoGE == iCryo);

    // initialize to the end directly; this has probably ID's isValid true
    geo::GeometryCore::cryostat_id_iterator iCryo2
      (geom, geo::CryostatID(geom->Ncryostats()));
    BOOST_TEST(iCryo2->Cryostat == geom->Ncryostats());
    BOOST_TEST(iCryo2 == iCryo);
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
    BOOST_TEST(iCryo == iCryoID);
    BOOST_TEST(!(iCryo != iCryoID));

    // ID comparison
    BOOST_TEST(iCryo.ID() == *iCryoID);

    // check copy assignment
    geo::GeometryCore::cryostat_iterator iCryo_copy(iCryo);
    geo::GeometryCore::cryostat_id_iterator iCryoID_copy(iCryoID);

    // check comparisons too
    BOOST_TEST(iCryo == iCryo_copy);
    BOOST_TEST(iCryo_copy == iCryo);
    BOOST_TEST(!(iCryo != iCryo_copy));
    BOOST_TEST(!(iCryo_copy != iCryo));

    BOOST_TEST(iCryo == iCryo_copy);

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

    BOOST_TEST(iTPC->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iTPC->TPC == geo::TPCID::getInvalidID());

    // check that the iterator tests false
    BOOST_TEST(!iTPC);
    BOOST_TEST(!(bool(iTPC)));

  }

  //
  // begin-constructed
  //
  {
    geo::GeometryCore::TPC_id_iterator iTPC(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created TPC ID iterator: " << std::string(*iTPC));

    BOOST_TEST(iTPC->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iTPC->TPC ==                geo::TPCID::TPCID_t(0));

    // check that the iterator tests true
    BOOST_TEST(bool(iTPC));
    BOOST_TEST(!!iTPC);

    // initialize to the beginning directly; this has probably ID's isValid true
    geo::TPCID BeginID;
    geom->GetBeginID(BeginID);
    geo::GeometryCore::TPC_id_iterator iTPCD(geom, BeginID);
    BOOST_TEST(iTPCD->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iTPCD->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iTPCD == iTPC);

    // construct from explicit begin position
    geo::GeometryCore::TPC_id_iterator iTPCBC
      (geom, geo::GeometryCore::TPC_id_iterator::begin_pos);
    BOOST_TEST(iTPCBC == iTPC);

    // construct at begin position by geometry
    geo::GeometryCore::TPC_id_iterator iTPCGB = geom->begin_TPC_id();
    BOOST_TEST(iTPCGB == iTPC);

    // check access to ID
    BOOST_TEST(*iTPC == BeginID);
    BOOST_TEST(iTPC->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iTPC->TPC == BeginID.TPC);

    // check access to geometry element
    geo::TPCGeo const* pTPC = geom->TPCPtr(BeginID);
    BOOST_TEST(iTPC.get() == pTPC);

    // test copy and postfix increment
    geo::GeometryCore::TPC_id_iterator iTPCI(iTPC++);

    const unsigned int nTPCsInC0 = geom->NTPC(geo::CryostatID(0));
    if (nTPCsInC0 > 1) {
      BOOST_TEST(iTPCI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iTPCI->TPC ==                geo::TPCID::TPCID_t(0));
      BOOST_TEST(iTPC->Cryostat ==  geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iTPC->TPC ==                 geo::TPCID::TPCID_t(1));
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

    geo::GeometryCore::TPC_id_iterator iTPC(geom, ID);

    // check that the iterator tests true
    BOOST_TEST(bool(iTPC));
    BOOST_TEST(!!iTPC);

    // check that the pointed ID is as expected
    BOOST_TEST(*iTPC == ID);
    BOOST_TEST(iTPC->Cryostat == ID.Cryostat);
    BOOST_TEST(iTPC->TPC == ID.TPC);
    BOOST_TEST(iTPC.get() == geom->TPCPtr(ID));

    ++iTPC;
    // check that the pointed ID is as expected
    BOOST_TEST(iTPC->Cryostat == geo::CryostatID::CryostatID_t(ID.Cryostat + 1));
    BOOST_TEST(iTPC->TPC ==                geo::TPCID::TPCID_t(0));


    // test iterator to last TPC
    geo::TPCID LastID(geom->Ncryostats() - 1, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    geo::GeometryCore::TPC_id_iterator iLastTPC(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last TPC ID: "
      << std::string(*iLastTPC));

    // check that the iterator tests true
    BOOST_TEST(bool(iLastTPC));
    BOOST_TEST(!!iLastTPC);

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastTPC == LastID);
    BOOST_TEST(iLastTPC->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastTPC->TPC == LastID.TPC);
    BOOST_TEST(iLastTPC.get() == geom->TPCPtr(LastID));

    // test increment to past-the-end
    geo::GeometryCore::TPC_id_iterator iEndTPC = iLastTPC;
    ++iEndTPC;

    // check that the iterator tests false
    BOOST_TEST(!bool(iEndTPC));
    BOOST_TEST(!iEndTPC);

    BOOST_TEST(iEndTPC->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndTPC->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iEndTPC == geom->end_TPC_id());
    BOOST_TEST(!iEndTPC.get());

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

    BOOST_TEST(iTPC->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iTPC->TPC ==                geo::TPCID::TPCID_t(0));

    // check that the iterator tests false
    BOOST_TEST(!bool(iTPC));
    BOOST_TEST(!iTPC);

    // check access to geometry element (result of operator* is not defined)
    BOOST_TEST(!(iTPC.get())); // should get nullptr

    // construct at end position by geometry
    geo::GeometryCore::TPC_id_iterator iTPCGE = geom->end_TPC_id();
    BOOST_TEST(iTPCGE == iTPC);

    // initialize to the end directly; this has probably ID's isValid true
    geo::GeometryCore::TPC_id_iterator iTPC2
      (geom, geo::TPCID(geom->Ncryostats(), 0));
    BOOST_TEST(iTPC2->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iTPC2->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iTPC2 == iTPC);

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
    BOOST_TEST(iTPC == iTPCID);
    BOOST_TEST(!(iTPC != iTPCID));

    // ID comparison
    BOOST_TEST(iTPC.ID() == *iTPCID);

    // check copy assignment
    geo::GeometryCore::TPC_iterator iTPC_copy(iTPC);
    geo::GeometryCore::TPC_id_iterator iTPCID_copy(iTPCID);

    // check comparisons too
    BOOST_TEST(iTPC == iTPC_copy);
    BOOST_TEST(iTPC_copy == iTPC);
    BOOST_TEST(!(iTPC != iTPC_copy));
    BOOST_TEST(!(iTPC_copy != iTPC));

    BOOST_TEST(iTPC == iTPC_copy);

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

    BOOST_TEST(iPlane->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iPlane->TPC ==           geo::TPCID::getInvalidID());
    BOOST_TEST(iPlane->Plane ==       geo::PlaneID::getInvalidID());

    // check that the iterator tests false
    BOOST_TEST(!iPlane);
    BOOST_TEST(!(bool(iPlane)));

  }

  //
  // begin-constructed
  //
  {
    geo::GeometryCore::plane_id_iterator iPlane(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created plane ID iterator: " << std::string(*iPlane));

    BOOST_TEST(iPlane->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iPlane->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iPlane->Plane ==          geo::PlaneID::PlaneID_t(0));

    // check that the iterator tests true
    BOOST_TEST(bool(iPlane));
    BOOST_TEST(!!iPlane);

    // initialize to the beginning directly; this has probably ID's isValid true
    geo::PlaneID BeginID;
    geom->GetBeginID(BeginID);
    geo::GeometryCore::plane_id_iterator iPlaneD(geom, BeginID);
    BOOST_TEST(iPlaneD->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iPlaneD->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iPlaneD->Plane ==          geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iPlaneD == iPlane);

    // construct from explicit begin position
    geo::GeometryCore::plane_id_iterator iPlaneBC
      (geom, geo::GeometryCore::plane_id_iterator::begin_pos);
    BOOST_TEST(iPlaneBC == iPlane);

    // construct at begin position by geometry
    geo::GeometryCore::plane_id_iterator iPlaneGB = geom->begin_plane_id();
    BOOST_TEST(iPlaneGB == iPlane);

    // check access to ID
    BOOST_TEST(*iPlane == BeginID);
    BOOST_TEST(iPlane->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iPlane->TPC == BeginID.TPC);
    BOOST_TEST(iPlane->Plane == BeginID.Plane);

    // check access to geometry element
    geo::PlaneGeo const* pPlane = geom->PlanePtr(BeginID);
    BOOST_TEST(iPlane.get() == pPlane);

    // test copy and postfix increment
    geo::GeometryCore::plane_id_iterator iPlaneI(iPlane++);

    const unsigned int nPlanesInC0T0 = geom->Nplanes(geo::TPCID(0, 0));
    if (nPlanesInC0T0 > 1) {
      BOOST_TEST(iPlaneI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iPlaneI->TPC ==                geo::TPCID::TPCID_t(0));
      BOOST_TEST(iPlaneI->Plane ==          geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iPlane->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iPlane->TPC ==                geo::TPCID::TPCID_t(0));
      BOOST_TEST(iPlane->Plane ==          geo::PlaneID::PlaneID_t(1));
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

    geo::GeometryCore::plane_id_iterator iPlane(geom, ID);

    // check that the iterator tests true
    BOOST_TEST(bool(iPlane));
    BOOST_TEST(!!iPlane);

    // check that the pointed ID is as expected
    BOOST_TEST(*iPlane == ID);
    BOOST_TEST(iPlane->Cryostat == ID.Cryostat);
    BOOST_TEST(iPlane->TPC == ID.TPC);
    BOOST_TEST(iPlane->Plane == ID.Plane);
    BOOST_TEST(iPlane.get() == geom->PlanePtr(ID));

    // check that the pointed ID is as expected
    ++iPlane;
    if (ID.TPC + 1 < geom->NTPC(ID)) {
      BOOST_TEST(iPlane->Cryostat == ID.Cryostat);
      BOOST_TEST(iPlane->TPC == ID.TPC + 1);
      BOOST_TEST(iPlane->Plane == geo::PlaneID::PlaneID_t(0));
    }
    else {
      BOOST_TEST(iPlane->Cryostat == ID.Cryostat + 1);
      BOOST_TEST(iPlane->TPC ==       geo::TPCID::TPCID_t(0));
      BOOST_TEST(iPlane->Plane == geo::PlaneID::PlaneID_t(0));
    }

    // test iterator to last plane
    geo::PlaneID LastID(geom->Ncryostats() - 1, 0, 0);
    LastID.TPC = geom->NTPC(LastID) - 1; // last TPC of last cryostat
    LastID.Plane = geom->Nplanes(LastID) - 1; // last plane of last TPC
    geo::GeometryCore::plane_id_iterator iLastPlane(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last plane ID: "
      << std::string(*iLastPlane));

    // check that the iterator tests true
    BOOST_TEST(bool(iLastPlane));
    BOOST_TEST(!!iLastPlane);

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastPlane == LastID);
    BOOST_TEST(iLastPlane->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastPlane->TPC == LastID.TPC);
    BOOST_TEST(iLastPlane->Plane == LastID.Plane);
    BOOST_TEST(iLastPlane.get() == geom->PlanePtr(LastID));

    // test increment to past-the-end
    geo::GeometryCore::plane_id_iterator iEndPlane = iLastPlane;
    ++iEndPlane;

    // check that the iterator tests false
    BOOST_TEST(!bool(iEndPlane));
    BOOST_TEST(!iEndPlane);

    BOOST_TEST(iEndPlane->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndPlane->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iEndPlane->Plane ==          geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iEndPlane == geom->end_plane_id());
    BOOST_TEST(!iEndPlane.get());

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

    BOOST_TEST(iPlane->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iPlane->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iPlane->Plane ==          geo::PlaneID::PlaneID_t(0));

    // check that the iterator tests false
    BOOST_TEST(!bool(iPlane));
    BOOST_TEST(!iPlane);

    // check access to geometry element (result of operator* is not defined)
    BOOST_TEST(!(iPlane.get())); // should get nullptr

    // construct at end position by geometry
    geo::GeometryCore::plane_id_iterator iPlaneGE = geom->end_plane_id();
    BOOST_TEST(iPlaneGE == iPlane);

    // initialize to the end directly; this has probably ID's isValid true
    geo::GeometryCore::plane_id_iterator iPlane2
      (geom, geo::PlaneID(geom->Ncryostats(), 0, 0));
    BOOST_TEST(iPlane2->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iPlane2->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iPlane2->Plane ==          geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iPlane2 == iPlane);

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
    BOOST_TEST(iPlane == iPlaneID);
    BOOST_TEST(!(iPlane != iPlaneID));

    // ID comparison
    BOOST_TEST(iPlane.ID() == *iPlaneID);

    // check copy assignment
    geo::GeometryCore::plane_iterator iPlane_copy(iPlane);
    geo::GeometryCore::plane_id_iterator iPlaneID_copy(iPlaneID);

    // check comparisons too
    BOOST_TEST(iPlane == iPlane_copy);
    BOOST_TEST(iPlane_copy == iPlane);
    BOOST_TEST(!(iPlane != iPlane_copy));
    BOOST_TEST(!(iPlane_copy != iPlane));

    BOOST_TEST(iPlane == iPlane_copy);

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

    BOOST_TEST(iWire->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iWire->TPC ==           geo::TPCID::getInvalidID());
    BOOST_TEST(iWire->Plane ==       geo::PlaneID::getInvalidID());
    BOOST_TEST(iWire->Wire ==         geo::WireID::getInvalidID());

    // check that the iterator tests false
    BOOST_TEST(!iWire);
    BOOST_TEST(!(bool(iWire)));

  }

  //
  // begin-constructed
  //
  {
    geo::GeometryCore::wire_id_iterator iWire(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created wire ID iterator: " << std::string(*iWire));

    BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iWire->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iWire->Plane ==          geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iWire->Wire ==             geo::WireID::WireID_t(0));

    // check that the iterator tests true
    BOOST_TEST(bool(iWire));
    BOOST_TEST(!!iWire);

    // initialize to the beginning directly; this has probably ID's isValid true
    geo::WireID BeginID;
    geom->GetBeginID(BeginID);
    geo::GeometryCore::wire_id_iterator iWireD(geom, BeginID);
    BOOST_TEST(iWireD->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iWireD->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iWireD->Plane ==          geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iWireD->Wire ==             geo::WireID::WireID_t(0));
    BOOST_TEST(iWireD == iWire);

    // construct from explicit begin position
    geo::GeometryCore::wire_id_iterator iWireBC
      (geom, geo::GeometryCore::wire_id_iterator::begin_pos);
    BOOST_TEST(iWireBC == iWire);

    // construct at begin position by geometry
    geo::GeometryCore::wire_id_iterator iWireGB = geom->begin_wire_id();
    BOOST_TEST(iWireGB == iWire);

    // check access to ID
    BOOST_TEST(*iWire == BeginID);
    BOOST_TEST(iWire->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iWire->TPC == BeginID.TPC);
    BOOST_TEST(iWire->Plane == BeginID.Plane);
    BOOST_TEST(iWire->Wire == BeginID.Wire);

    // check access to geometry element
    geo::WireGeo const* pWire = geom->WirePtr(BeginID);
    BOOST_TEST(iWire.get() == pWire);

    // test copy and postfix increment
    geo::GeometryCore::wire_id_iterator iWireI(iWire++);

    const unsigned int nWiresInC0T0P0 = geom->Nwires(geo::PlaneID(0, 0, 0));
    if (nWiresInC0T0P0 > 1) {
      BOOST_TEST(iWireI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iWireI->TPC ==                geo::TPCID::TPCID_t(0));
      BOOST_TEST(iWireI->Plane ==          geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iWireI->Wire ==             geo::WireID::WireID_t(0));
      BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iWire->TPC ==                geo::TPCID::TPCID_t(0));
      BOOST_TEST(iWire->Plane ==          geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iWire->Wire ==             geo::WireID::WireID_t(1));
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

    geo::GeometryCore::wire_id_iterator iWire(geom, ID);

    // check that the iterator tests true
    BOOST_TEST(bool(iWire));
    BOOST_TEST(!!iWire);

    // check that the pointed ID is as expected
    BOOST_TEST(*iWire == ID);
    BOOST_TEST(iWire->Cryostat == ID.Cryostat);
    BOOST_TEST(iWire->TPC == ID.TPC);
    BOOST_TEST(iWire->Plane == ID.Plane);
    BOOST_TEST(iWire->Wire == ID.Wire);
    BOOST_TEST(iWire.get() == geom->WirePtr(ID));

    ++iWire;
    // check that the pointed ID is as expected
    if (ID.Plane + 1 < geom->Nplanes(ID)) {
      BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iWire->TPC ==                geo::TPCID::TPCID_t(0));
      BOOST_TEST(iWire->Plane ==                                  ID.Plane + 1);
      BOOST_TEST(iWire->Wire ==             geo::WireID::WireID_t(0));
    } else if (ID.TPC + 1 < geom->NTPC(ID)) {
      BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iWire->TPC ==                                    ID.TPC + 1);
      BOOST_TEST(iWire->Plane ==          geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iWire->Wire ==             geo::WireID::WireID_t(0));
    } else {
      BOOST_TEST(iWire->Cryostat ==                               ID.Cryostat + 1);
      BOOST_TEST(iWire->TPC ==                geo::TPCID::TPCID_t(0));
      BOOST_TEST(iWire->Plane ==          geo::PlaneID::PlaneID_t(0));
      BOOST_TEST(iWire->Wire ==             geo::WireID::WireID_t(0));
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
    BOOST_TEST(bool(iLastWire));
    BOOST_TEST(!!iLastWire);

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastWire == LastID);
    BOOST_TEST(iLastWire->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastWire->TPC == LastID.TPC);
    BOOST_TEST(iLastWire->Plane == LastID.Plane);
    BOOST_TEST(iLastWire->Wire == LastID.Wire);
    BOOST_TEST(iLastWire.get() == geom->WirePtr(LastID));

    // test increment to past-the-end
    geo::GeometryCore::wire_id_iterator iEndWire = iLastWire;
    ++iEndWire;

    // check that the iterator tests false
    BOOST_TEST(!bool(iEndWire));
    BOOST_TEST(!iEndWire);

    BOOST_TEST(iEndWire->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndWire->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iEndWire->Plane ==          geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iEndWire->Wire ==             geo::WireID::WireID_t(0));
    BOOST_TEST(iEndWire == geom->end_wire_id());
    BOOST_TEST(!iEndWire.get());

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

    BOOST_TEST(iWire->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iWire->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iWire->Plane ==          geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iWire->Wire ==             geo::WireID::WireID_t(0));

    // check that the iterator tests false
    BOOST_TEST(!bool(iWire));
    BOOST_TEST(!iWire);

    // check access to geometry element (result of operator* is not defined)
    BOOST_TEST(!(iWire.get())); // should get nullptr

    // construct at end position by geometry
    geo::GeometryCore::wire_id_iterator iWireGE = geom->end_wire_id();
    BOOST_TEST(iWireGE == iWire);

    // initialize to the end directly; this has probably ID's isValid true
    geo::GeometryCore::wire_id_iterator iWire2
      (geom, geo::WireID(geom->Ncryostats(), 0, 0, 0));
    BOOST_TEST(iWire2->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iWire2->TPC ==                geo::TPCID::TPCID_t(0));
    BOOST_TEST(iWire2->Plane ==          geo::PlaneID::PlaneID_t(0));
    BOOST_TEST(iWire2->Wire ==             geo::WireID::WireID_t(0));
    BOOST_TEST(iWire2 == iWire);

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
    BOOST_TEST(iWire == iWireID);
    BOOST_TEST(!(iWire != iWireID));

    // ID comparison
    BOOST_TEST(iWire.ID() == *iWireID);

    // check copy assignment
    geo::GeometryCore::wire_iterator iWire_copy(iWire);
    geo::GeometryCore::wire_id_iterator iWireID_copy(iWireID);

    // check comparisons too
    BOOST_TEST(iWire == iWire_copy);
    BOOST_TEST(iWire_copy == iWire);
    BOOST_TEST(!(iWire != iWire_copy));
    BOOST_TEST(!(iWire_copy != iWire));

    BOOST_TEST(iWire == iWire_copy);

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

    BOOST_TEST(iTPCset->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iTPCset->TPCset == readout::TPCsetID::getInvalidID());

    // check that the iterator tests false
    BOOST_TEST(!iTPCset);
    BOOST_TEST(!(bool(iTPCset)));

  }

  //
  // begin-constructed
  //
  {
    geo::TPCset_id_iterator iTPCset(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created TPC set ID iterator: " << std::string(*iTPCset));

    BOOST_TEST(iTPCset->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iTPCset->TPCset ==   readout::TPCsetID::TPCsetID_t(0));

    // check that the iterator tests true
    BOOST_TEST(bool(iTPCset));
    BOOST_TEST(!!iTPCset);

    // initialize to the beginning directly; this has probably ID's isValid true
    readout::TPCsetID BeginID;
    geom->GetBeginID(BeginID);
    geo::TPCset_id_iterator iTPCsetD(geom, BeginID);
    BOOST_TEST(iTPCsetD->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iTPCsetD->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iTPCsetD == iTPCset);

    // construct from explicit begin position
    geo::TPCset_id_iterator iTPCsetBC(geom, geo::iterators::begin_pos);
    BOOST_TEST(iTPCsetBC == iTPCset);

    // construct at begin position by geometry
    geo::TPCset_id_iterator iTPCsetGB = geom->begin_TPCset_id();
    BOOST_TEST(iTPCsetGB == iTPCset);

    // check access to ID
    BOOST_TEST(*iTPCset == BeginID);
    BOOST_TEST(iTPCset->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iTPCset->TPCset == BeginID.TPCset);

    // test copy and postfix increment
    geo::TPCset_id_iterator iTPCsetI(iTPCset++);

    const unsigned int nTPCsetsInC0 = geom->NTPCsets(geo::CryostatID(0));
    if (nTPCsetsInC0 > 1) {
      BOOST_TEST(iTPCsetI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iTPCsetI->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
      BOOST_TEST(iTPCset->Cryostat ==  geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iTPCset->TPCset ==    readout::TPCsetID::TPCsetID_t(1));
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
    ID.TPCset = geom->NTPCsets(ID) - 1; // last TPC set of first cryostat

    geo::TPCset_id_iterator iTPCset(geom, ID);

    // check that the iterator tests true
    BOOST_TEST(bool(iTPCset));
    BOOST_TEST(!!iTPCset);

    // check that the pointed ID is as expected
    BOOST_TEST(*iTPCset == ID);
    BOOST_TEST(iTPCset->Cryostat == ID.Cryostat);
    BOOST_TEST(iTPCset->TPCset == ID.TPCset);

    ++iTPCset;
    // check that the pointed ID is as expected
    BOOST_TEST
      (iTPCset->Cryostat == geo::CryostatID::CryostatID_t(ID.Cryostat + 1));
    BOOST_TEST(iTPCset->TPCset == readout::TPCsetID::TPCsetID_t(0));


    // test iterator to last TPC
    readout::TPCsetID LastID(geom->Ncryostats() - 1, 0);
    LastID.TPCset = geom->NTPCsets(LastID) - 1; // last TPC set of last cryostat
    geo::TPCset_id_iterator iLastTPCset(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last TPC set ID: "
      << std::string(*iLastTPCset));

    // check that the iterator tests true
    BOOST_TEST(bool(iLastTPCset));
    BOOST_TEST(!!iLastTPCset);

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastTPCset == LastID);
    BOOST_TEST(iLastTPCset->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastTPCset->TPCset == LastID.TPCset);

    // test increment to past-the-end
    geo::TPCset_id_iterator iEndTPCset = iLastTPCset;
    ++iEndTPCset;

    // check that the iterator tests false
    BOOST_TEST(!bool(iEndTPCset));
    BOOST_TEST(!iEndTPCset);

    BOOST_TEST
      (iEndTPCset->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndTPCset->TPCset == readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iEndTPCset == geom->end_TPCset_id());

  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::TPCset_id_iterator iTPCset(geom, geo::iterators::end_pos);
    BOOST_TEST_CHECKPOINT
      ("End-created TPC set ID iterator: " << std::string(*iTPCset));

    BOOST_TEST
      (iTPCset->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iTPCset->TPCset == readout::TPCsetID::TPCsetID_t(0));

    // check that the iterator tests false
    BOOST_TEST(!bool(iTPCset));
    BOOST_TEST(!iTPCset);

    // construct at end position by geometry
    geo::TPCset_id_iterator iTPCsetGE = geom->end_TPCset_id();
    BOOST_TEST(iTPCsetGE == iTPCset);

    // initialize to the end directly; this has probably ID's isValid true
    geo::TPCset_id_iterator iTPCset2
      (geom, readout::TPCsetID(geom->Ncryostats(), 0));
    BOOST_TEST
      (iTPCset2->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iTPCset2->TPCset == readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iTPCset2 == iTPCset);

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

    BOOST_TEST(iROP->Cryostat == geo::CryostatID::getInvalidID());
    BOOST_TEST(iROP->TPCset ==   readout::TPCsetID::getInvalidID());
    BOOST_TEST(iROP->ROP ==      readout::ROPID::getInvalidID());

    // check that the iterator tests false
    BOOST_TEST(!iROP);
    BOOST_TEST(!(bool(iROP)));

  }

  //
  // begin-constructed
  //
  {
    geo::ROP_id_iterator iROP(geom);
    BOOST_TEST_CHECKPOINT
      ("Begin-created readout plane ID iterator: " << std::string(*iROP));

    BOOST_TEST(iROP->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iROP->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iROP->ROP ==      readout::ROPID::ROPID_t(0));

    // check that the iterator tests true
    BOOST_TEST(bool(iROP));
    BOOST_TEST(!!iROP);

    // initialize to the beginning directly; this has probably ID's isValid true
    readout::ROPID BeginID;
    geom->GetBeginID(BeginID);
    geo::ROP_id_iterator iROPD(geom, BeginID);
    BOOST_TEST(iROPD->Cryostat == geo::CryostatID::CryostatID_t(0));
    BOOST_TEST(iROPD->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iROPD->ROP ==      readout::ROPID::ROPID_t(0));
    BOOST_TEST(iROPD == iROP);

    // construct from explicit begin position
    geo::ROP_id_iterator iROPBC(geom, geo::iterators::begin_pos);
    BOOST_TEST(iROPBC == iROP);

    // construct at begin position by geometry
    geo::ROP_id_iterator iROPGB = geom->begin_ROP_id();
    BOOST_TEST(iROPGB == iROP);

    // check access to ID
    BOOST_TEST(*iROP == BeginID);
    BOOST_TEST(iROP->Cryostat == BeginID.Cryostat);
    BOOST_TEST(iROP->TPCset == BeginID.TPCset);
    BOOST_TEST(iROP->ROP == BeginID.ROP);

    // test copy and postfix increment
    geo::ROP_id_iterator iROPI(iROP++);

    const unsigned int nReadoutPlanesInC0S0
      = geom->NROPs(readout::TPCsetID(0, 0));
    if (nReadoutPlanesInC0S0 > 1) {
      BOOST_TEST(iROPI->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iROPI->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
      BOOST_TEST(iROPI->ROP ==            readout::ROPID::ROPID_t(0));
      BOOST_TEST(iROP->Cryostat == geo::CryostatID::CryostatID_t(0));
      BOOST_TEST(iROP->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
      BOOST_TEST(iROP->ROP ==            readout::ROPID::ROPID_t(1));
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
    ID.ROP = geom->NROPs(ID) - 1; // last plane of first TPC set

    geo::ROP_id_iterator iROP(geom, ID);

    // check that the iterator tests true
    BOOST_TEST(bool(iROP));
    BOOST_TEST(!!iROP);

    // check that the pointed ID is as expected
    BOOST_TEST(*iROP == ID);
    BOOST_TEST(iROP->Cryostat == ID.Cryostat);
    BOOST_TEST(iROP->TPCset == ID.TPCset);
    BOOST_TEST(iROP->ROP == ID.ROP);

    // check that the pointed ID is as expected
    ++iROP;
    if (ID.TPCset + 1 < (int) geom->NTPCsets(ID)) {
      BOOST_TEST(iROP->Cryostat == ID.Cryostat);
      BOOST_TEST(iROP->TPCset ==   ID.TPCset + 1);
      BOOST_TEST(iROP->ROP ==      readout::ROPID::ROPID_t(0));
    }
    else {
      BOOST_TEST(iROP->Cryostat == ID.Cryostat + 1);
      BOOST_TEST(iROP->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
      BOOST_TEST(iROP->ROP ==      readout::ROPID::ROPID_t(0));
    }

    // test iterator to last plane
    readout::ROPID LastID(geom->Ncryostats() - 1, 0, 0);
    LastID.TPCset = geom->NTPCsets(LastID) - 1; // last TPC set of last cryostat
    LastID.ROP = geom->NROPs(LastID) - 1; // last readout plane of last TPC set
    geo::ROP_id_iterator iLastROP(geom, LastID);
    BOOST_TEST_CHECKPOINT("Position-created iterator to last readout plane ID: "
      << std::string(*iLastROP));

    // check that the iterator tests true
    BOOST_TEST(bool(iLastROP));
    BOOST_TEST(!!iLastROP);

    // check that the pointed ID is as expected
    BOOST_TEST(*iLastROP == LastID);
    BOOST_TEST(iLastROP->Cryostat == LastID.Cryostat);
    BOOST_TEST(iLastROP->TPCset ==   LastID.TPCset);
    BOOST_TEST(iLastROP->ROP ==      LastID.ROP);

    // test increment to past-the-end
    geo::ROP_id_iterator iEndROP = iLastROP;
    ++iEndROP;

    // check that the iterator tests false
    BOOST_TEST(!bool(iEndROP));
    BOOST_TEST(!iEndROP);

    BOOST_TEST
      (iEndROP->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iEndROP->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iEndROP->ROP ==      readout::ROPID::ROPID_t(0));
    BOOST_TEST(iEndROP == geom->end_ROP_id());

  }

  //
  // end-constructed
  //
  {
    // construct from end position
    geo::ROP_id_iterator iROP(geom, geo::iterators::end_pos);
    BOOST_TEST_CHECKPOINT
      ("End-created readout plane ID iterator: " << std::string(*iROP));

    BOOST_TEST
      (iROP->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iROP->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iROP->ROP ==            readout::ROPID::ROPID_t(0));

    // check that the iterator tests false
    BOOST_TEST(!bool(iROP));
    BOOST_TEST(!iROP);

    // construct at end position by geometry
    geo::ROP_id_iterator iROPGE = geom->end_ROP_id();
    BOOST_TEST(iROPGE == iROP);

    // initialize to the end directly; this has probably ID's isValid true
    geo::ROP_id_iterator iROP2
      (geom, readout::ROPID(geom->Ncryostats(), 0, 0));
    BOOST_TEST
      (iROP->Cryostat == geo::CryostatID::CryostatID_t(geom->Ncryostats()));
    BOOST_TEST(iROP->TPCset ==   readout::TPCsetID::TPCsetID_t(0));
    BOOST_TEST(iROP->ROP ==            readout::ROPID::ROPID_t(0));
    BOOST_TEST(iROP2 == iROP);

  }
} // GeometryIteratorTestAlg::ROPIDIteratorsTest()
