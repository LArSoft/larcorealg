/**
 * @addtogroup Utilities General utilities
 * @brief General programming utilities.
 *
 * @{
 */
/**
 * @defgroup Metaprogramming General utilities for metaprogramming
 * @brief General utilities for use with templates and metaprogramming.
 */
/**
 * @}
 */

/** ****************************************************************************
 * @file   larcorealg/CoreUtils/MetaUtils.h
 * @brief  Basic C++ metaprogramming utilities.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 12, 2017
 * @ingroup Metaprogramming
 */

#ifndef LARCOREALG_COREUTILS_METAUTILS_H
#define LARCOREALG_COREUTILS_METAUTILS_H

// C/C++ standard libraries
#include <array>
#include <functional> // std::reference_wrapper<>
#include <memory>     // std::addressof()
#include <string>
#include <string_view>
#include <type_traits>

/**
 * @namespace util
 * @brief Namespace for general, non-LArSoft-specific utilities.
 * @ingroup Utilities
 */
namespace util {

  //--- BEGIN MetaprogrammingBase ----------------------------------------------
  /**
   * @defgroup MetaprogrammingBase Simple utility traits
   * @brief Simple traits for the implementation of other traits.
   * @ingroup Metaprogramming
   */
  /// @{

  namespace details {

    /// Implementation detail of `staticDumpClassName()`.
    template <typename T>
    struct ClassNameStaticDumper;

  } // namespace details

  //----------------------------------------------------------------------------
  /// Trait returning the very same type as in the template argument.
  template <typename T>
  struct self_type {
    using type = T;
  };

  /// The very same type as in the template argument.
  template <typename T>
  using self_t = typename self_type<T>::type;

  //----------------------------------------------------------------------------
  /**
   * @brief A `std::false_type` with a template argument.
   * @see util::always_true_type, util::always_false_v
   *
   * This type allows a `static_assert` to fail only when the template type it's
   * in is being instantiated:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct MandatoryCustomizationPoint {
   *   static_assert(util::always_false_type<T>(), "You have to specialize!");
   * };
   *
   * template <typename T>
   * struct MandatoryCustomizationPoint<std::vector<T>> {
   *   using type = typename std::vector<T>::reference;
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this example, using `std::false_type` instead of
   * `util::always_false_type` might have tripped the compiler to trigger the
   * assertion failure even if the class is not instantiated.
   */
  template <typename>
  struct always_false_type : public std::false_type {};

  //----------------------------------------------------------------------------
  /**
   * @brief A templated constant, always false.
   * @see util::always_false_type, util::always_true_v
   *
   * This constant allows a `static_assert` to fail only when the template type
   * it's in is being instantiated:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct MandatoryCustomizationPoint {
   *   static_assert(util::always_false_v<T>, "You have to specialize!");
   * };
   *
   * template <typename T>
   * struct MandatoryCustomizationPoint<std::vector<T>> {
   *   using type = typename std::vector<T>::reference;
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this example, using `std::false_type` might have tipped the compiler to
   * trigger the assertion failure even if the class is not instantiated.
   */
  template <typename>
  constexpr bool always_false_v = false;

  /**
   * @brief A `std::true_type` with a template argument.
   * @see util::always_false_type, util::always_true_v
   *
   * This is one way to allow to specialize for classes with a certain type:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T, typename = void>
   * class ReferenceTypeExtractor {
   *   static_assert(util::always_false_type<T>(), "Type has no reference!");
   * };
   *
   * template <typename Cont>
   * struct ReferenceTypeExtractor<
   *   Cont,
   *   std::enable_if_t
   *     <util::always_true_type<typename Cont::value_type>::value>
   *   >
   * {
   *   using type = typename Cont::reference;
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename>
  struct always_true_type : public std::true_type {};

  /**
   * @brief A template constant always true.
   * @see util::always_true_type, util::always_false_v
   *
   * This is one way to allow to specialize for classes with a certain type:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T, typename = void>
   * class ReferenceTypeExtractor {
   *   static_assert(util::always_false_v<T>, "Type has no reference!");
   * };
   *
   * template <typename Cont>
   * struct ReferenceTypeExtractor
   *   <Cont, std::enable_if_t<util::always_true_v<typename Cont::value_type>>>
   * {
   *   using type = typename Cont::reference;
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename>
  constexpr bool always_true_v = true;

  /// The negation of `std::is_same`.
  template <typename A, typename B>
  using is_not_same = std::negation<std::is_same<A, B>>;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait: index of the first occurrence of `T` among the specified
   *        `Types`, starting from the one with index `StartFrom`.
   * @tparam T the type of check the presence of
   * @tparam StartFrom number of `Types` that will be ignored
   * @tparam Types the possible types `T` can match.
   * @see    `util::find_next_type`
   *
   * The value of the trait is the index of `T` within the specified list of
   * `Types` (first type as index `0`).
   * The first `StartFrom` `Types` are ignored, but still counted.
   * The match is exact, as in `std::is_same`.
   * If none of the `Types` exactly matches `T`, the trait value will be the
   * number of types (i.e. `sizeof...(Types)`), which is the index after the
   * last of the types.
   *
   * This is a integral trait (type `std::size_t`): use it as
   * `std::integer_constant`.
   */
  template <typename T, std::size_t StartFrom, typename... Types>
  struct find_next_type;

  template <typename T, std::size_t StartFrom, typename... Types>
  constexpr std::size_t find_next_type_v = find_next_type<T, StartFrom, Types...>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait: index of the first occurrence of `T` among the specified
   *        `Types`.
   * @tparam T the type of check the presence of
   * @tparam Types the possible types `T` can match.
   * @see    `util::find_next_type`
   *
   * The value of the trait is the index of `T` within the specified list of
   * `Types` (first type as index `0`). The match is exact, as in
   * `std::is_same`. If none of the `Types` exactly matches `T`, the trait value
   * will be the number of types (i.e. `sizeof...(Types)`), which is the index
   * after the last of the types.
   *
   * This is a integral trait (type `std::size_t`): use it as
   * `std::integer_constant`.
   */
  template <typename T, typename... Types>
  using find_type = find_next_type<T, 0U, Types...>;

  /// Index of the first occurrence of `T` among the specified `Types`.
  /// @see `util::find_type`
  template <typename T, typename... Types>
  constexpr std::size_t find_type_v = find_type<T, Types...>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait: whether `T` is among the specified `Types`.
   * @tparam T the type of check the presence of
   * @tparam Types the possible types `T` can match.
   *
   * Matching is for the exact type, as in `std::is_same`.
   *
   * This is a boolean trait: use it as `std::bool_constant`.
   */
  template <typename T, typename... Types>
  struct is_any_of;

  /// Whether `T` is among the specified `Types` (see `util::is_any_of`).
  template <typename T, typename... Types>
  constexpr bool is_any_of_v = is_any_of<T, Types...>::value;

  //----------------------------------------------------------------------------
  /// Whether `T` and `U` are the same type, after being applied `std::decay`.
  template <typename T, typename U>
  constexpr auto is_same_decay_v = std::is_same_v<std::decay_t<T>, std::decay_t<U>>;

  //----------------------------------------------------------------------------
  // @{
  /**
   * @brief Helper to determine the type of a variable at compilation time.
   * @tparam T type to be investigated
   *
   * It may be difficult to understand which type is being used in a failing
   * static assertion or in some complicate metaprogramming code (is there any
   * other kind?), especially when due to a compilation failure the code can't
   * be run.
   * This class is supposed to help by forcing the compiler to a halt, and it is
   * devised so that the compiler should print in the error message what it
   * thinks the type `T` is.
   *
   * An example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * #include "larcorealg/CoreUtils/MetaUtils.h"
   *
   * void f() {
   *   constexpr auto v = 5U - 6U; // which type is `v` of?
   *   util::staticDumpClassName(v);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * For example, Clang 5.0.1 emits these errors:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In file included from metatest.cpp:1:
   * ./MetaUtils.h:217:7: error: static_assert failed "ClassNameStaticDumper<T>: look for T in the error message context"
   *       static_assert(
   *       ^
   * ./MetaUtils.h:228:39: note: in instantiation of template class 'util::details::ClassNameStaticDumper<unsigned int>' requested here
   *   void staticDumpClassName() { (void) details::ClassNameStaticDumper<T>(); }
   *                                       ^
   * ./MetaUtils.h:195:46: note: in instantiation of function template specialization 'util::staticDumpClassName<unsigned int>' requested here
   *   [[noreturn]] void staticDumpClassName(T) { staticDumpClassName<T>(); }
   *                                              ^
   * metatest.cpp:5:16: note: in instantiation of function template specialization 'util::staticDumpClassName<unsigned int>' requested here
   *   (void) util::staticDumpClassName(v);
   *                ^
   * 1 error generated.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * From the first "note" we can see that the type of `v` is `unsigned int`.
   * Note that if `v` is not copiable, an additional error will be emitted.
   * To avoid that, the type can be specified as template parameter, as in
   * `util::staticDumpClassName<decltype(v)>()`.
   * The same program is processed by GNU GCC with an error message:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In file included from metatest.cpp:1:0:
   * MetaUtils.h: In instantiation of ‘struct util::details::ClassNameStaticDumper<unsigned int>’:
   * MetaUtils.h:248:48:   required from ‘void util::staticDumpClassName() [with T = unsigned int]’
   * MetaUtils.h:215:68:   required from ‘void util::staticDumpClassName(T) [with T = unsigned int]’
   * metatest.cpp:5:30:   required from here
   * MetaUtils.h:237:7: error: static assertion failed: ClassNameStaticDumper<T>: look for T in the error message context
   *        static_assert(
   *        ^~~~~~~~~~~~~
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * where the type is mentioned in all the three context lines.
   */
  template <typename T>
  void staticDumpClassName();

  template <typename T>
  void staticDumpClassName(T)
  {
    staticDumpClassName<T>();
  }
  // @}

  /// @}
  //--- END MetaprogrammingBase ------------------------------------------------

  //--- BEGIN Type identification ----------------------------------------------
  /**
   * @defgroup MetaprogrammingTypeIdentification Determination of specific types
   * @brief Traits to identify specific types.
   * @ingroup Metaprogramming
   */
  /// @{

  //----------------------------------------------------------------------------

  /**
   * @brief Trait describing whether `T` is a template instance of `Template`.
   * @tparam Template template class to be detected
   * @tparam T type to be tested
   *
   * This trait is true if the type `T` is an instance of template class
   * `Template`, that is, if `T` is `Template<...>`, with the ellipsis
   * represents any template argument. Before being tested, `T` is stripped
   * of reference and constantness qualifiers, so that for example the answer
   * will be the same for `std::vector<int>` as for `std::vector<int> const`,
   * `std::vector<int> volatile&`, etc.
   *
   * @bug The limitation of this implementation is that only `Template` types
   *      taking only type arguments can be used. For example, attempting to
   *      use it with `std::array`, which contains a non-type argument (of type
   *      `std::size_t`), will cause a compilation error. For example, GCC 7.2
   *      reports:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * error: type/value mismatch at argument 1 in template parameter list for ‘template<template<class ...> class Template, class T> constexpr bool const util::is_instance_of_v<Template, T>’
   * static_assert(util::is_instance_of_v<std::array, std::array<int, 2U>>);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  template <template <typename...> typename Template, typename T>
  struct is_instance_of;

  /// A constant describing whether `T` is a template instance of `Template`.
  /// @see `util::is_instance_of`
  template <template <typename...> typename Template, typename T>
  constexpr bool is_instance_of_v = is_instance_of<Template, T>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Identifies whether the specified type is a STL array.
   * @tparam T the type to be tested
   * @see `util::is_STLarray_v`
   */
  template <typename T>
  struct is_STLarray;

  /**
   * @brief A constant describing whether the specified type is a STL array.
   * @tparam T the type to be tested
   * @see `util::is_STLarray`
   */
  template <typename T>
  constexpr bool is_STLarray_v = is_STLarray<T>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Identifies whether the specified type is a `std::reference_wrapper`.
   * @tparam T the type to be tested
   * @see `util::is_reference_wrapper_v`
   */
  template <typename T>
  using is_reference_wrapper = is_instance_of<std::reference_wrapper, T>;

  /**
   * @brief A constant describing whether the specified type is a
   *        `std::reference_wrapper`.
   * @tparam T the type to be tested
   * @see `util::is_reference_wrapper`
   */
  template <typename T>
  constexpr bool is_reference_wrapper_v = is_reference_wrapper<T>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Identifies whether the specified type is a `std::unique_ptr`.
   * @tparam T the type to be tested
   * @see `util::is_unique_ptr_v`
   */
  template <typename T>
  using is_unique_ptr = is_instance_of<std::unique_ptr, T>;

  /**
   * @brief A constant describing whether the specified type is a
   *        `std::unique_ptr`.
   * @tparam T the type to be tested
   * @see `util::is_unique_ptr`
   */
  template <typename T>
  constexpr bool is_unique_ptr_v = is_unique_ptr<T>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait: whether type `T` is a character type.
   * @tparam T the type to be tested
   *
   * Character types are `char` (in all its sign options), `wchar_t`, `char32_t`
   * and `char16_t`, in any combination of constantness and volatility.
   * References to types yield the same value as the types they reference.
   */
  template <typename T>
  struct is_character_type;

  /// Whether type `T` is a character type (see `util::is_character_type`).
  template <typename T>
  constexpr bool is_character_type_v = is_character_type<T>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait: whether type `T` is a character string type.
   * @tparam T the type to be tested
   * @see `util::is_character_type`
   *
   * In this definition, any container of character types is a string.
   * A container is defined as a type having a `value_type` member.
   * Also, C-style arrays and pointers to characters are considered strings.
   * Reference types yield the same value as their referenced type.
   */
  template <typename T>
  struct is_string_type;

  /// Whether type `T` is a character string type (see `util::is_string_type`).
  template <typename T>
  constexpr bool is_string_type_v = is_string_type<T>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait: whether type `T` is a STL string type.
   * @tparam T the type to be tested
   *
   * This trait has `value` `true` if `T` is an instance of `std::basic_string`
   * template, `false` otherwise.
   */
  template <typename T>
  struct is_basic_string_type;

  /// Whether type `T` is a character string type
  /// (see `util::is_basic_string_type`).
  template <typename T>
  constexpr bool is_basic_string_type_v = is_basic_string_type<T>::value;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait: whether type `T` is a `std::string_view` type.
   * @tparam T the type to be tested
   *
   * This trait has `value` `true` if `T` is an instance of
   * `std::basic_string_view` template, `false` otherwise.
   */
  template <typename T>
  struct is_basic_string_view_type;

  /// Whether type `T` is a character string type
  /// (see `util::is_basic_string_view_type`).
  template <typename T>
  constexpr bool is_basic_string_view_type_v = is_basic_string_view_type<T>::value;

  /// @}
  //--- END Type identification ------------------------------------------------

  //--- BEGIN Type manipulation ------------------------------------------------
  /**
   * @defgroup MetaprogrammingTypeManipulation Manipulation of types
   * @brief Traits to change types.
   * @ingroup Metaprogramming
   */
  /// @{

  //----------------------------------------------------------------------------
  /**
   * @brief Trait with type `Base`, plus the constantness as in `Key`.
   * @tparam Base the basic type being returned
   * @tparam Key a type expressing the constantness wanted for `Base`
   *
   * The `type` member of this trait is:
   * * the type `Base` with constantness removed, if `Key` is non-constant type
   * * the type `Base` with constantness added, if `Key` is a constant type
   *
   * This trait passes through references. Both `Base` and `Key` are treated as
   * they were no references (e.g. a `int const&` is treated as a `int const`).
   * The referenceness of the type in the trait is the same as the one of
   * `Base`.
   *
   * Therefore, for example:
   * * `with_const_as<int const, double>` yields `int`;
   * * `with_const_as<int const, double const>` yields `int const`;
   * * `with_const_as<int const&, double>` yields `int&`;
   * * `with_const_as<int, double const&>` yields `int const`.
   *
   * R-value references are likewise taken into account.
   */
  template <typename Base, typename Key>
  struct with_const_as;

  /**
   * @brief The type `Base`, plus the constantness as in `Key`.
   * @tparam Base the basic type being returned
   * @tparam Key a type expressing the constantness wanted for `Base`
   * @see `util::with_const_as`
   */
  template <typename Base, typename Key>
  using with_const_as_t = typename with_const_as<Base, Key>::type;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait with type `T` stripped of all known reference types.
   * @tparam T type to remove referenceness from
   *
   * In addition of the standard C++ references, this trait also removes all
   * pseudo-referenceness known (to it).
   * That currently includes:
   * * C++ l-value and r-value reference
   * * `std::reference_wrapper`
   *
   * The implementation removes all references recursively.
   */
  template <typename T>
  struct strip_referenceness_type;

  /**
   * @brief The type `T` stripped of all known reference types.
   * @tparam T type to remove referenceness from
   * @see `util::strip_referenceness_type`
   */
  template <typename T>
  using strip_referenceness_t = typename strip_referenceness_type<T>::type;

  //----------------------------------------------------------------------------
  /**
   * @brief Returns the address of the referenced object.
   * @tparam Ref type of reference
   * @param ref reference
   * @return a pointer to the referenced object
   *
   * This function also manages `std::reference_wrapper` arguments, by returning
   * the address of the object they reference.
   * In all other cases, the return value is simply `std::addressof(obj)`.
   */
  template <typename Ref>
  auto referenced_address(Ref&& ref);

  //----------------------------------------------------------------------------
  /**
   * @brief Functor applying the proper `referenced_address()` function.
   * @see `referenced_address()`
   *
   * This class operates in the same way as `util::referenced_address()`, but it
   * is easier to use in STL algorithms since it's not a template:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> data(4U, 0);
   * std::vector<int const*> dataPtr;
   * std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtr),
   *   util::reference_addresser());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will fill `dataPtr` with the pointers to the elements in `data`.
   * This is easier than taking the address of the correct template instance of
   * `util::referenced_address()`, or than writing a lambda function for that.
   */
  struct reference_addresser {
    template <typename Ref>
    decltype(auto) operator()(Ref&& ref) const
    {
      return addressof(std::forward<Ref>(ref));
    }

    template <typename Ref>
    static decltype(auto) addressof(Ref&& ref)
    {
      return referenced_address(std::forward<Ref>(ref));
    }

  }; // struct reference_addresser

  //----------------------------------------------------------------------------
  /**
   * @brief Trait with type `T` into `std::reference_wrapper` if reference.
   * @tparam T type to be wrapped
   * @see `util::lvalue_reference_into_wrapper()`
   *
   * If the argument type `T` is a l-value reference, the corresponding `type`
   * will be a `std::reference_wrapper` object wrapping the type `T` references
   * (constantness will be preserved).
   * If the argument is already a `std::reference_wrapper` (possibly constant,
   * possibly any reference), no action is taken and `type` is the same as `T`
   * except for no reference.
   * Otherwise, `type` will be `T` itself with any reference removed.
   */
  template <typename T>
  struct lvalue_reference_into_wrapper_type;

  /**
   * @brief The type `T` stripped of all known reference types.
   * @tparam T type to remove referenceness from
   * @see `util::lvalue_reference_into_wrapper_type`
   */
  template <typename T>
  using lvalue_reference_into_wrapper_t = typename lvalue_reference_into_wrapper_type<T>::type;

  /**
   * @brief Converts a l-value reference object into a `std::reference_wrapper`.
   * @tparam T type of the object to be converted
   * @param obj object to be converted
   * @return either `obj` or a `std::reference_wrapper` around it
   * @see `util::lvalue_reference_into_wrapper_type`
   *
   * This function operates on a fashion similar to
   * `util::lvalue_reference_into_wrapper_type`, but it performs the conversion
   * of an instantiated object rather than just reporting a type.
   */
  template <typename T>
  auto lvalue_reference_into_wrapper(T&& obj)
  {
    return util::lvalue_reference_into_wrapper_t<T>(obj);
  }

  //----------------------------------------------------------------------------

  /// @}
  //--- END Type manipulation --------------------------------------------------

} // namespace util

//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
namespace util {

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    template <std::size_t Index, std::size_t Skip, typename T, typename... Types>
    struct find_type_impl;

    template <std::size_t Index, std::size_t Skip, typename T, typename Type, typename... Others>
    struct find_type_impl<Index, Skip, T, Type, Others...>
      : std::integral_constant<
          std::size_t,
          (Skip == 0) && std::is_same_v<T, Type> ?
            Index :
            find_type_impl<Index + 1U, ((Skip > 0U) ? Skip - 1U : 0U), T, Others...>::value> {};

    template <std::size_t Index, std::size_t Skip, typename T>
    struct find_type_impl<Index, Skip, T> : std::integral_constant<std::size_t, Index> {};

    //--------------------------------------------------------------------------
    template <typename T, typename = void>
    struct is_character_type_impl
      : std::bool_constant<util::is_any_of_v<std::decay_t<T>,
                                             signed char,
                                             unsigned char,
                                             char,
                                             wchar_t,
#ifdef __cpp_char8_t // C++20
                                             char8_t,
#endif                                                 // __cpp_char8_t
                                             char16_t, // this is defined unsigned
                                             char32_t  // this is defined unsigned
                                             >> {
    };

    //--------------------------------------------------------------------------
    template <typename T, typename = void>
    struct is_string_type_impl : std::false_type {};

    template <typename T>
    struct is_string_type_impl<
      T,
      std::enable_if_t<is_character_type_impl<typename T::value_type>::value>> : std::true_type {};

    template <typename T>
    struct is_string_type_impl<
      T,
      std::enable_if_t<std::is_pointer_v<std::decay_t<T>> &&
                       is_character_type_impl<std::remove_pointer_t<std::decay_t<T>>>::value>>
      : std::true_type {};

    template <typename T>
    struct is_string_type_impl<
      T,
      std::enable_if_t<std::is_array_v<std::decay_t<T>> &&
                       is_character_type_impl<std::remove_extent_t<std::decay_t<T>>>::value>>
      : std::true_type {};

    //--------------------------------------------------------------------------
    template <typename T>
    struct is_basic_string_type_impl : std::false_type {};

    template <typename... Args>
    struct is_basic_string_type_impl<std::basic_string<Args...>> : std::true_type {};

    //--------------------------------------------------------------------------
    template <typename T>
    struct is_basic_string_view_type_impl : std::false_type {};

    template <typename... Args>
    struct is_basic_string_view_type_impl<std::basic_string_view<Args...>> : std::true_type {};

    //--------------------------------------------------------------------------
    /// Implementation detail of `staticDumpClassName()`.
    template <typename T>
    struct ClassNameStaticDumper {
      static_assert(util::always_false_type<T>(),
                    "ClassNameStaticDumper<T>: look for T in the error message context");
    }; // struct ClassNameStaticDumper

    //--------------------------------------------------------------------------
    // implementation for `with_const_as`

    // - final implementation:
    template <typename Base, typename /* Key */, typename = void>
    struct with_const_as_impl {
      using type = std::remove_const_t<Base>;
    };

    template <typename Base, typename Key>
    struct with_const_as_impl<Base, Key, std::enable_if_t<std::is_const_v<Key>>> {
      using type = std::add_const_t<Base>;
    };

    // - implementation dispatcher for reference types
    //   - pass through for not-reference types
    template <typename Base, typename Key, typename = void>
    struct with_const_as_dispatch_ref : with_const_as_impl<Base, Key> {};
    //   - lvalue reference
    template <typename Base, typename Key>
    struct with_const_as_dispatch_ref<Base,
                                      Key,
                                      std::enable_if_t<std::is_lvalue_reference_v<Base>>> {
      using type = std::add_lvalue_reference_t<
        typename with_const_as_impl<std::remove_reference_t<Base>, Key>::type>;
    };
    //   - rvalue reference
    template <typename Base, typename Key>
    struct with_const_as_dispatch_ref<Base,
                                      Key,
                                      std::enable_if_t<std::is_rvalue_reference_v<Base>>> {
      using type = std::add_rvalue_reference_t<
        typename with_const_as_impl<std::remove_reference_t<Base>, Key>::type>;
    };

    // - key management
    template <typename Base, typename Key>
    struct with_const_as_dispatch_keyref
      : with_const_as_dispatch_ref<Base, std::remove_reference_t<Key>> {};

    // - top level implementation dispatcher
    template <typename Base, typename Key>
    struct with_const_as_dispatcher : with_const_as_dispatch_keyref<Base, Key> {};

    //--------------------------------------------------------------------------
    //--- implementation of `is_instance_of`
    template <template <typename...> typename Template, typename T>
    struct is_instance_of_impl : std::false_type {};

    template <template <typename...> typename Template, typename... Args>
    struct is_instance_of_impl<Template, Template<Args...>> : std::true_type {};

    //--------------------------------------------------------------------------
    //--- implementation of `strip_referenceness_type`

    template <typename T>
    struct strip_referenceness_type_impl;

    // implementation layer dealing with `std::reference_wrapper`
    template <typename T, typename = void>
    struct strip_referenceness_type_impl_wrapref {
      using type = T;
    }; // exit here

    // - handle any constantness and volatility
    template <typename T>
    struct strip_referenceness_type_impl_wrapref<
      T,
      std::enable_if_t<util::is_reference_wrapper_v<std::remove_cv_t<T>>>>
      : strip_referenceness_type_impl<typename T::type> // back to square one
    {};

    // implementation layer dealing with C++ references
    template <typename T>
    struct strip_referenceness_type_impl_ref : strip_referenceness_type_impl_wrapref<T> {};

    template <typename T>
    struct strip_referenceness_type_impl_ref<T&> : strip_referenceness_type_impl_wrapref<T> {};

    template <typename T>
    struct strip_referenceness_type_impl_ref<T&&> : strip_referenceness_type_impl_wrapref<T> {};

    // entry point: start by dealing with C++ references
    template <typename T>
    struct strip_referenceness_type_impl : strip_referenceness_type_impl_ref<T> {};

    //--------------------------------------------------------------------------
    //--- implementation of `referenced_address()`

    template <typename T, typename = void>
    struct referenced_address_impl {
      static auto addressof(T& obj) { return std::addressof(obj); }
    };

    template <typename T>
    struct referenced_address_impl<T, std::enable_if_t<util::is_reference_wrapper_v<T>>> {
      static auto addressof(T& obj) { return std::addressof(obj.get()); }
    };

    //--------------------------------------------------------------------------
    //--- implementation of `lvalue_reference_into_wrapper_type`
    /*
     * The implementation develops across levels:
     *  1. any kind of `std::reference_wrapper` is handled
     *  2. l-value references are handled
     */
    template <typename T>
    struct lvalue_reference_into_wrapper_type_impl_final {
      using type = std::remove_reference_t<T>;
    };

    template <typename T>
    struct lvalue_reference_into_wrapper_type_impl_final<T&> {
      using type = std::reference_wrapper<T>;
    };

    template <typename T, typename = void>
    struct lvalue_reference_into_wrapper_type_impl_wrapper {
      using type = typename lvalue_reference_into_wrapper_type_impl_final<T>::type;
    };

    template <typename T>
    struct lvalue_reference_into_wrapper_type_impl_wrapper<
      T,
      std::enable_if_t<util::is_reference_wrapper_v<T>>> {
      using type = std::remove_reference_t<T>;
    };

    template <typename T>
    struct lvalue_reference_into_wrapper_type_impl
      : lvalue_reference_into_wrapper_type_impl_wrapper<T> {};

    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  template <template <typename...> typename Template, typename T>
  struct is_instance_of : details::is_instance_of_impl<Template, std::decay_t<T>> {};

  //----------------------------------------------------------------------------
  template <typename T, std::size_t StartFrom, typename... Types>
  struct find_next_type : details::find_type_impl<0U, StartFrom, T, Types...> {};

  //----------------------------------------------------------------------------
  template <typename T, typename... Types>
  struct is_any_of : std::bool_constant<((find_type_v<T, Types...>) < sizeof...(Types))> {};

  //----------------------------------------------------------------------------
  template <typename T>
  struct is_character_type : details::is_character_type_impl<T> {};

  //----------------------------------------------------------------------------
  template <typename T>
  struct is_string_type : details::is_string_type_impl<T> {};

  //----------------------------------------------------------------------------
  template <typename T>
  struct is_basic_string_type : details::is_basic_string_type_impl<std::decay_t<T>> {};

  //----------------------------------------------------------------------------
  template <typename T>
  struct is_basic_string_view_type : details::is_basic_string_view_type_impl<std::decay_t<T>> {};

  //----------------------------------------------------------------------------
  template <typename T>
  void staticDumpClassName()
  {
    (void)details::ClassNameStaticDumper<T>();
  }

  //----------------------------------------------------------------------------
  template <typename>
  struct is_STLarray : public std::false_type {};

  template <typename T, std::size_t N>
  struct is_STLarray<std::array<T, N>> : public std::true_type {};

  //----------------------------------------------------------------------------
  template <typename Base, typename Key>
  struct with_const_as : public details::with_const_as_dispatcher<Base, Key> {};

  //----------------------------------------------------------------------------
  template <typename T>
  struct strip_referenceness_type : public details::strip_referenceness_type_impl<T> {};

  //----------------------------------------------------------------------------
  template <typename T>
  struct lvalue_reference_into_wrapper_type
    : public details::lvalue_reference_into_wrapper_type_impl<T> {};

  //----------------------------------------------------------------------------
  template <typename Ref>
  auto referenced_address(Ref&& ref)
  {
    return details::referenced_address_impl<Ref>::addressof(ref);
  }

  //----------------------------------------------------------------------------

} // namespace util

//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_METAUTILS_H
