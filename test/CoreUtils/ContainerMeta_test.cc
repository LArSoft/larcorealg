/**
 * @file    ContainerMeta_test.cc
 * @brief   Unit test for some of the utilities in `ContainerMeta.h`
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 25, 2019
 * @see     `larcorealg/CoreUtils/ContainerMeta.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE ( ContainerMeta_test )
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h"

// C/C++ standard libraries
#include <list>
#include <memory> // std::unique_ptr<>
#include <functional> // std::reference_wrapper
#include <type_traits>


//------------------------------------------------------------------------------
//--- static tests (would fail at compile time)
//------------------------------------------------------------------------------
struct O {
  int value = 0;
};

static_assert(std::is_same_v<util::collection_value_t                <std::list<O                              >      >, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<O      *                       >      >, O      *       >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<O const*                       >      >, O const*       >);

static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O                              >      >, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O      *                       >      >, O      *      &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O const*                       >      >, O const*      &>);

static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O                              >      >, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O      *                       >      >, O      * const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O const*                       >      >, O const* const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::list<O                              > const>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<O      *                       > const>, O      *       >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<O const*                       > const>, O const*       >);

static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O                              > const>, O        const&>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O      *                       > const>, O      * const&>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O const*                       > const>, O const* const&>);

static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O                              > const>, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O      *                       > const>, O      * const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O const*                       > const>, O const* const&>);

//
// std::reference_wrapper content
//
static_assert(std::is_same_v<util::collection_value_t                <std::list<std::reference_wrapper<O      >>      >, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<std::reference_wrapper<O      >> const>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<std::reference_wrapper<O const>>      >, O const        >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<std::reference_wrapper<O const>> const>, O const        >);

//
// std::reference_wrapper container
//
static_assert(std::is_same_v<util::collection_value_t                <                       std::list<O      >&      >, O              >);
static_assert(std::is_same_v<util::collection_value_t                <                       std::list<O      > const&>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::reference_wrapper<std::list<O      >>      >, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::reference_wrapper<std::list<O      >>     &>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::reference_wrapper<std::list<O      >>const&>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::reference_wrapper<std::list<O> const>const&>, O              >);

static_assert(std::is_same_v<util::collection_value_access_t         <                       std::list<O      >&      >, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <                       std::list<O      > const&>, O        const&>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::reference_wrapper<std::list<O      >>      >, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::reference_wrapper<std::list<O      >>     &>, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::reference_wrapper<std::list<O      >>const&>, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::reference_wrapper<std::list<O> const>const&>, O        const&>);

static_assert(std::is_same_v<util::collection_value_constant_access_t<                       std::list<O      >&      >, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<                       std::list<O      > const&>, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::reference_wrapper<std::list<O      >>      >, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::reference_wrapper<std::list<O      >>     &>, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::reference_wrapper<std::list<O      >>const&>, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::reference_wrapper<std::list<O> const>const&>, O        const&>);

//
// C array container
//
static_assert(std::is_same_v<util::collection_value_t                <O      *>, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <O      *>, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<O      *>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <O const*>, O const >);
static_assert(std::is_same_v<util::collection_value_access_t         <O const*>, O const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<O const*>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <O[10U]  >, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <O[10U]  >, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<O[10U]  >, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::unique_ptr<O      >>, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <std::unique_ptr<O      >>, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::unique_ptr<O      >>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::unique_ptr<O const>>, O const >);
static_assert(std::is_same_v<util::collection_value_access_t         <std::unique_ptr<O const>>, O const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::unique_ptr<O const>>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::unique_ptr<O[10U] >>, O         [10U]>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::unique_ptr<O[10U] >>, O(&)      [10U]>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::unique_ptr<O[10U] >>, O const(&)[10U]>);

static_assert(std::is_same_v<util::collection_value_t                <std::unique_ptr<O[]    >>, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <std::unique_ptr<O[]    >>, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::unique_ptr<O[]    >>, O const&>);

//
// util::collection_reference_t
//
static_assert(std::is_same_v<util::collection_reference_t<                       std::list<O>   >, std::reference_wrapper<std::list<O>      >>);
static_assert(std::is_same_v<util::collection_reference_t<                       std::list<O>&  >, std::reference_wrapper<std::list<O>      >>);
static_assert(std::is_same_v<util::collection_reference_t<std::reference_wrapper<std::list<O>>  >, std::reference_wrapper<std::list<O>      >>);
static_assert(std::is_same_v<util::collection_reference_t<                       O*             >,                                  O*       >);
static_assert(std::is_same_v<util::collection_reference_t<                       O const*       >,                                  O const* >);
static_assert(std::is_same_v<util::collection_reference_t<                       O const*      &>,                                  O const* >);
static_assert(std::is_same_v<util::collection_reference_t<                       O const* const&>,                                  O const* >);
static_assert(std::is_same_v<util::collection_reference_t<                       O[10U]         >,                                  O*       >);


//------------------------------------------------------------------------------
void make_collection_reference_test() {

  /*
   * `std::unique_ptr<T[N]>` is a strange beast which really deals with pointers
   * to whole C arrays (`T(*)[N]`).
   * It should be understood that this is not the same as `T*`
   */
  constexpr std::size_t Size = 10U;

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::list<O>    list        {{ O{} }};
  O               array       [Size];
  O*              ptr         = array;
  O const*        cptr        = ptr;
  O const* const  cptrc       = ptr;
  O const* const& cptrcr      = cptrc;
  auto            refw        = std::ref(list);
  auto            crefw       = std::cref(list);
  auto const      uptr        = std::unique_ptr<O>(new O[Size]);
  auto            array_uptr  = std::unique_ptr<O[Size]>();
  auto            garray_uptr = std::unique_ptr<O[]>(new O[Size]);

  static_assert(std::is_same_v<decltype(util::make_collection_reference(list       )), decltype(std::ref(list))>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(array      )), decltype(ptr           )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(ptr        )), decltype(ptr           )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(cptr       )), decltype(cptr          )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(cptrc      )), decltype(cptr          )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(cptrcr     )), decltype(cptr          )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(refw       )), decltype(refw          )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(crefw      )), decltype(crefw         )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(uptr       )), decltype(ptr           )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(array_uptr )), decltype(&array        )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(garray_uptr)), decltype(ptr           )>);

  BOOST_TEST(&(util::make_collection_reference(list       ).get().front()) == &list.front()      );
  BOOST_TEST(&(util::make_collection_reference(array      )[0]           ) == &array[0]          );
  BOOST_TEST(&(util::make_collection_reference(ptr        )[0]           ) ==  ptr               );
  BOOST_TEST(&(util::make_collection_reference(cptr       )[0]           ) ==  cptr              );
  BOOST_TEST(&(util::make_collection_reference(cptrc      )[0]           ) ==  cptr              );
  BOOST_TEST(&(util::make_collection_reference(cptrcr     )[0]           ) ==  cptr              );
  BOOST_TEST(&(util::make_collection_reference(refw       ).get().front()) == &list.front()      );
  BOOST_TEST(&(util::make_collection_reference(crefw      ).get().front()) == &list.front()      );
  BOOST_TEST(&(util::make_collection_reference(uptr       )[0]           ) ==  uptr.get()        );
  BOOST_TEST(&(util::make_collection_reference(array_uptr )[0]           ) ==  array_uptr.get()  );
  BOOST_TEST(&(util::make_collection_reference(garray_uptr)[0]           ) ==  garray_uptr.get() );

} // make_collection_reference_test()


//------------------------------------------------------------------------------
void collection_from_reference_test() {

  constexpr std::size_t Size = 10U;

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::list<O>        list        {{ O{} }};
  std::list<O>      & ref         = list;
  std::list<O> const& cref        = list;
  auto                refw        = std::ref(list);
  auto                crefw       = std::cref(list);

  O                   array       [Size];
  O*                  ptr         = array;
  O const*            cptr        = ptr;

  auto const          uptr        = std::unique_ptr<O>(new O[Size]);
  auto                array_uptr  = std::unique_ptr<O[Size]>();
  auto const          array_uptrc = std::unique_ptr<O[Size]>();
  auto                garray_uptr = std::unique_ptr<O[]>(new O[Size]);

  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(list       )>, decltype(ref )>);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(ref        )>, decltype(ref )>);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(cref       )>, decltype(cref)>);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(refw       )>, decltype(ref )>);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(crefw      )>, decltype(cref)>);

  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(array      )>, O      *      >);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(ptr        )>, O      *      >);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(cptr       )>, O const*      >);

  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(uptr       )>, O      *     >);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(array_uptr )>, O(*)[Size]   >);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(array_uptrc)>, O(*)[Size]   >);
  static_assert(std::is_same_v<util::collection_from_reference_t<decltype(garray_uptr)>, O      *     >);

  static_assert(std::is_same_v<decltype(util::collection_from_reference(list       )), decltype(ref )>);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(ref        )), decltype(ref )>);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(cref       )), decltype(cref)>);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(refw       )), decltype(ref )>);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(crefw      )), decltype(cref)>);

  static_assert(std::is_same_v<decltype(util::collection_from_reference(array      )), O      *      >);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(ptr        )), O      *      >);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(cptr       )), O const*      >);

  static_assert(std::is_same_v<decltype(util::collection_from_reference(uptr       )), O      *     >);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(array_uptr )), O(*)[Size]   >);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(array_uptrc)), O(*)[Size]   >);
  static_assert(std::is_same_v<decltype(util::collection_from_reference(garray_uptr)), O      *     >);

  BOOST_TEST(&util::collection_from_reference(list       ) == &list            );
  BOOST_TEST(&util::collection_from_reference(ref        ) == &list            );
  BOOST_TEST(&util::collection_from_reference(cref       ) == &list            );
  BOOST_TEST(&util::collection_from_reference(refw       ) == &list            );
  BOOST_TEST(&util::collection_from_reference(crefw      ) == &list            );

  BOOST_TEST( util::collection_from_reference(array      ) == array            );
  BOOST_TEST( util::collection_from_reference(ptr        ) == ptr              );
  BOOST_TEST( util::collection_from_reference(cptr       ) == cptr             );

  BOOST_TEST( util::collection_from_reference(uptr       ) == uptr.get()       );
  BOOST_TEST( util::collection_from_reference(array_uptr ) == array_uptr.get() );
  BOOST_TEST( util::collection_from_reference(array_uptrc) == array_uptrc.get());
  BOOST_TEST( util::collection_from_reference(garray_uptr) == garray_uptr.get());

} // collection_from_reference_test()


//------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(CollectionReferenceTestCase) {
  make_collection_reference_test();
  collection_from_reference_test();
} // BOOST_AUTO_TEST_CASE(CollectionReferenceTestCase)

//------------------------------------------------------------------------------
