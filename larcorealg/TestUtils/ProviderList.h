/**
 * @file   ProviderList.h
 * @brief  Container for service providers used in a test
 * @date   April 22nd, 2016
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @see    ProviderTesterHelpers.h
 *
 * This is a header-only library.
 * It depends only on standard C++ and on LArSoft pure headers with that same
 * feature, and therefore does not require additional linkage.
 *
 * This library provides:
 *
 * - ProviderList, an object managing service providers
 *
 */

#ifndef LARCORE_TESTUTILS_PROVIDERLIST_H
#define LARCORE_TESTUTILS_PROVIDERLIST_H 1

// LArSoft libraries
#include "larcorealg/TestUtils/ProviderTestHelpers.h"

// art libraries
#include "cetlib_except/demangle.h"

// C/C++ standard libraries
#include <memory>    // std::unique_ptr()
#include <stdexcept> // std::runtime_error
#include <typeinfo>
#include <unordered_map>
#include <utility> // std::forward()

namespace testing {

  namespace details {

    /// A base class with a virtual table
    struct MovableClassWrapperBase {
      virtual ~MovableClassWrapperBase() = default;

    }; // struct MovableClassWrapperBase

    /**
       * @brief A class containing an owned object
       * @tparam T type of the contained object
       *
       * Differently from std::any, `T` does not have to be copiable nor
       * movable. The price is that the object is dynamically allocated.
       *
       * Given that we don't require T to be movable and we don't require it to
       * be polymorphic, we need some way to move T and have T recognized as
       * such when put in a mixed container. Hence this wrapper.
       * The price is quite high: two chained dynamic allocations per object
       * (one in the wrapper to provide mobility, the other in the container
       * not to lose polymorphism).
       *
       */
    template <typename T>
    class MovableClassWrapper : public details::MovableClassWrapperBase {
      using datum_t = T;                          ///< contained type
      using this_t = MovableClassWrapper<T>;      ///< this type
      using pointer_t = std::shared_ptr<datum_t>; ///< pointer storing datum

      template <typename U>
      friend class MovableClassWrapper;

    public:
      struct share_t {}; /// type to specify a share constructor

      static constexpr share_t share = {};

      /// Default constructor: no datum present (move one in later on)
      MovableClassWrapper() : ptr(std::make_unique<datum_t>()) {}

      template <typename U>
      MovableClassWrapper(std::unique_ptr<U>&& from) : ptr(std::move(from))
      {}

      //@{
      /// Constructor and assignment from a unique pointer: steal the content
      MovableClassWrapper(pointer_t&& from) : ptr(std::move(from)) {}
      MovableClassWrapper& operator=(pointer_t&& from)
      {
        ptr = std::move(from);
        return *this;
      }
      //@}

      /// Share constructor (kind of copy)
      template <typename U>
      MovableClassWrapper(MovableClassWrapper<U> const& from, share_t) : ptr(from.ptr)
      {}

      /// Returns a reference to the datum with the correct type
      datum_t& ref() { return *ptr; }

      /// Returns a constant reference to the datum with the correct type
      datum_t const& ref() const { return *ptr; }

      /// Returns a pointer to the datum with the correct type
      datum_t* get() { return ptr.get(); }

      /// Returns a constant pointer to the datum with the correct type
      datum_t const* get() const { return ptr.get(); }

      //@{
      /// Returns whether there is a valid pointer
      bool valid() const { return bool(ptr); }

      operator bool() const { return valid(); }
      //@}

    private:
      pointer_t ptr;

    }; // namespace testing

  } // namespace details

  /** *************************************************************************
    * @brief Container of service providers accessed by type and optional label
    *
    * This container is expected to contain elements that are service providers
    * of different types. Each provider is accessed by its class type and
    * an optional instance label to discriminate between providers of the same
    * type.
    *
    * The list owns the providers. A provider is created with `setup()`
    * (or `setup_instance()` if a instance label is needed). This method
    * relies on a class `testing::ProvideSetupClass` to create and
    * correctly set up the provider. For example, to set up the provider
    * `LArPropertiesStandard`:
    *
    *     provList.setup<LArPropertiesStandard>(pset);
    *
    * assuming that `LArPropertiesStandard` provider has a constructor with
    * as only argument `pset` (supposedly, a `fhicl::ParameterSet`).
    * If a custom setup is needed, the methods `custom_setup_instance()` and
    * `custom_setup()`  take as argument a setup function, which can do
    * whatever it takes to perform the set up.
    *
    * After a provider is set up, a reference to it can be obtained by `get()`:
    *
    *     auto& larp = provList.get<LArPropertiesStandard>();
    *
    * If no such class is available, an exception will be thrown.
    * The presence of a provider can be checked beforehand with `has()`.
    *
    * @note The presence of multiple providers of the same type, supported via
    * instance names, is not useful in the current art/LArSoft.
    *
    *
    * How can a provider support the `setup()` construction method?
    * --------------------------------------------------------------
    *
    * A provider can specify its own setup by specialising the class
    * `testing::ProvideSetupClass`. A default implementation is provided,
    * that constructs the provider with a parameter set.
    *
    */
  class ProviderList {
    // Sparse implementation notes:
    // - we use MovableClassWrapperBase in place of std::any because our
    //   providers are recommended to be not copiable

    /// type of smart pointer we use to store elements
    template <typename T>
    using smart_pointer_t = std::unique_ptr<T>;

    /// Type of objects contained in the list
    using pointer_t = smart_pointer_t<details::MovableClassWrapperBase>;

    /// Type of list element with explicit element type memory
    template <typename T>
    using concrete_type_t = details::MovableClassWrapper<std::decay_t<T>>;
    /// Type of smart pointer to typed list element
    template <typename T>
    using concrete_pointer_t = smart_pointer_t<concrete_type_t<T>>;

  public:
    /// base exception class for ProviderList
    struct exception : public std::runtime_error {
      using std::runtime_error::runtime_error;
    };
    /// Exception thrown on a request about an unregistered type
    struct provider_not_available : public exception {
      using exception::exception;
    };
    /// Exception thrown on when object is not available any more
    struct provider_deleted : public exception {
      using exception::exception;
    };
    /// Exception thrown on a invalid type request
    struct provider_wrong : public exception {
      using exception::exception;
    };

    /**
       * @brief Construct and register an object of type T
       * @tparam T type of the object to be constructed (caller specifies it)
       * @tparam SetupProc type of functor performing the actual setup
       * @tparam Args type of constructor arguments (compiler fills them in)
       * @param label name of this instance of object type T (can be empty)
       * @param provSetup functor performing the setup
       * @param args arguments to provSetup for the construction of T
       *
       * An object is instantiated and associates it with the specified instance
       * label. It can then be accessed with a `get<T>(label)` call.
       *
       * The functor `provSetup` is expected to return a unique pointer to the
       * newly created provider, `std::unique_ptr<T>`.
       */
    template <typename T, typename SetupProc, typename... Args>
    bool custom_setup_instance(std::string label, SetupProc&& provSetup, Args&&... args)
    {
      auto k = key<T>(label); // key
      auto it = data.find(k);
      if (it != data.end()) return false;

      pointer_t ptr = std::make_unique<concrete_type_t<T>>(provSetup(std::forward<Args>(args)...));
      data.emplace_hint(it, std::move(k), std::move(ptr));
      return true;
    } // custom_setup_instance()

    /// Construct and register an object of type T with specified arguments
    template <typename T, typename SetupProc, typename... Args>
    bool custom_setup(SetupProc&& provSetup, Args&&... args)
    {
      return custom_setup_instance<T, SetupProc, Args...>(
        "", std::forward<SetupProc>(provSetup), std::forward<Args>(args)...);
    } // custom_setup()

    template <typename T, typename... Args>
    bool setup_instance(std::string label, Args&&... args)
    {
      auto k = key<T>(label); // key
      auto it = data.find(k);
      if (it != data.end()) return false;

      pointer_t ptr =
        std::make_unique<concrete_type_t<T>>(setupProvider<T>(std::forward<Args>(args)...));
      data.emplace_hint(it, std::move(k), std::move(ptr));
      return true;
    } // setup_instance()

    /// Construct and register an object of type T with specified arguments
    template <typename T, typename... Args>
    bool setup(Args&&... args)
    {
      return setup_instance<T>("", std::forward<Args>(args)...);
    }

    /**
       * @brief Registers and gets ownership of the specified object
       * @tparam T type of object being acquired
       * @param obj_ptr pointer to the object to be acquired
       * @param label name of the object instance
       * @return whether the object was acquired or not
       *
       * The ProviderList takes ownership of the specified provider.
       * If an object of type T is already registered, the pointer is left
       * untouched and `false` is returned.
       */
    template <typename T>
    bool acquire(std::unique_ptr<T>&& obj_ptr, std::string label = "")
    {
      auto k = key<T>(label); // key
      auto it = data.find(k);
      if (it != data.end()) return false;

      pointer_t ptr = std::make_unique<concrete_type_t<T>>(std::move(obj_ptr));
      data.emplace_hint(it, std::move(k), std::move(ptr));
      return true;
    } // acquire()

    /**
       * @brief Drops the object with the specified type and label
       * @tparam T type of object being acquired
       * @param label name of the object instance
       * @return whether the object was present or not
       *
       * If present, the object is destroyed
       */
    template <typename T>
    bool erase(std::string label = "")
    {
      auto k = key<T>(label); // key
      auto target_it = data.find(k);
      if (target_it == data.end()) return false;

      // erase this and all the aliases pointing to it
      auto const* target_ptr = target_it->second.get();
      auto it = data.begin();
      while (it != data.end()) {
        if (it->second.get() == target_ptr)
          it = data.erase(it);
        else
          ++it;
      } // while
      return true;
    } // erase()

    /// Sets the Alias type as an alias of the Prov provider (with labels)
    template <typename Prov, typename Alias>
    bool set_alias(std::string alias_label = "", std::string prov_label = "")
    {
      // find the alias location
      auto alias_k = key<Alias>(alias_label); // key
      auto alias_it = data.find(alias_k);
      if (alias_it != data.end()) return false;

      // find the original provider location
      auto prov_elem = get_elem<Prov>(prov_label);

      // register the shared object to the alias
      data.emplace_hint(alias_it,
                        std::move(alias_k),
                        std::make_unique<concrete_type_t<Alias>>(
                          prov_elem, typename concrete_type_t<Alias>::share_t()));
      return true;
    } // set_alias()

    /// @{
    /**
       * @brief Retrieve the object of type T stored with the specified label
       * @tparam T type of the object to be retrieved
       * @param label optional label used when the object was inserted
       * @return the specified object as a reference to type T
       * @throw provider_not_available no type T class was stored with label
       * @throw provider_deleted the object that was stored is not present
       * @throw provider_wrong the object is not compatible with the type T
       */
    template <typename T>
    T const& get(std::string label = "") const
    {
      return get_elem<T>(label).ref();
    }

    template <typename T>
    T& get(std::string label = "")
    {
      return get_elem<T>(label).ref();
    }
    /// @}

    /// @{
    /**
       * @brief Retrieve the object of type T stored with the specified label
       * @tparam T type of the object to be retrieved
       * @param label optional label used when the object was inserted
       * @return the specified object as a pointer to type T
       * @throw provider_not_available no type T class was stored with label
       * @throw provider_deleted the object that was stored is not present
       * @throw provider_wrong the object is not compatible with the type T
       */
    template <typename T>
    T const* getPointer(std::string label = "") const
    {
      return get_elem<T>(label).get();
    }

    template <typename T>
    T* getPointer(std::string label = "")
    {
      return get_elem<T>(label).get();
    }
    /// @}

    /// Returns whether we have a slot for this object
    template <typename T>
    bool known(std::string label = "") const
    {
      return find<T>(label) != data.end();
    }

    /// Returns whether the specified object is available
    template <typename T>
    bool valid(std::string label = "") const
    {
      auto it = find<T>(label);
      return (it != data.end()) && bool(it->second);
    }

  private:
    using key_type = size_t; ///< type used for key in the internal registry

    std::unordered_map<key_type, pointer_t> data; ///< all our singletons

    /// Convert a type into a demangled type name
    template <typename T>
    static std::string type_name()
    {
      return cet::demangle_symbol(typeid(T).name());
    }

    /// Convert a pointer to object into a demangled type name
    template <typename T>
    static std::string type_name(T const*)
    {
      return type_name<T>();
    }

    /// @{
    /// Returns an iterator pointing to the requested key, or data.end()
    template <typename T>
    auto find(std::string label = "") const
    {
      return data.find(key<T>(label));
    }

    template <typename T>
    auto find(std::string label = "")
    {
      return data.find(key<T>(label));
    }
    /// @}

    template <typename T>
    concrete_type_t<T> const& get_elem(std::string label = "") const
    {
      auto it = find<T>(label);
      if (it == data.end()) throw provider_not_available("Not available: " + type_name<T>());
      if (!(it->second)) throw provider_deleted("Deleted: " + type_name<T>());
      auto* ptr = dynamic_cast<details::MovableClassWrapper<T>*>(it->second.get());
      if (!ptr) {
        throw provider_wrong("Wrong: " + type_name(it->second.get()) +
                             " [requested: " + type_name<T>() + "]");
      }
      return *ptr;
    } // get_elem()

    template <typename T>
    concrete_type_t<T>& get_elem(std::string label = "")
    {
      auto it = find<T>(label);
      if (it == data.end()) throw provider_not_available("Not available: " + type_name<T>());
      if (!(it->second)) throw provider_deleted("Deleted: " + type_name<T>());
      auto* ptr = dynamic_cast<details::MovableClassWrapper<T>*>(it->second.get());
      if (!ptr) {
        throw provider_wrong("Wrong: " + type_name(it->second.get()) +
                             " [requested: " + type_name<T>() + "]");
      }
      return *ptr;
    } // get_elem()

    /// Extracts and returns the key out of a type and label
    template <typename T>
    static key_type key(std::string label = "")
    {
      return typeid(std::decay_t<T>).hash_code() ^ std::hash<std::string>()(label);
    }

  }; // class ProviderList

} // namespace testing

#endif // LARCORE_TESTUTILS_PROVIDERLIST_H
