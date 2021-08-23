/**
 * @file   larcorealg/CoreUtils/DebugUtils.h
 * @brief  Functions to help debugging by instrumenting code.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   April 8, 2016
 * @see    `larcorealg/CoreUtils/DebugUtils.cxx`
 *
 * This library contains:
 *  - a function to return the name of the type of a variable
 *  - a function printing into a stream the current call stack
 *
 */
#ifndef LARCOREALG_COREUTILS_DEBUGUTILS_H
#define LARCOREALG_COREUTILS_DEBUGUTILS_H

// LArSoft includes
#include "larcorealg/CoreUtils/MetaUtils.h"

// framework and support libraries
#include "cetlib_except/demangle.h"

// C/C++ standard libraries
#include <cstddef> // std::ptrdiff_t
#include <cstdlib> // std::free()
#include <utility> // std::pair<>
#include <vector>
#include <string>
#include <bitset>
#include <typeinfo>
#include <ostream>
#include <utility> // std::forward()

// non-standard libraries
#include <execinfo.h> // backtrace()...
// #include <experimental/filesystem> // std::experimental::filesystem::path


namespace lar::debug {

  /** ***********************************************************************
   * @brief Outputs a demangled name for type T.
   * @tparam T type whose name must be demangled (optional)
   * @return a string with demangled name
   *
   * It relies on cetlib.
   * The type to be demangled can be specified either as template argument:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto name = lar::debug::demangle<std::string>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * or via a argument pointer:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto name = lar::debug::demangle(this);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  template <typename T>
  inline std::string demangle(T const* = nullptr);


  //----------------------------------------------------------------------------
  /// Structure with information about a single call, parsed.
  struct CallInfo_t {
      private:
    using range_t = std::pair<size_t, size_t>;

      public:
    CallInfo_t(std::string const& s) { ParseString(s); }
    CallInfo_t(const char* s) { ParseString(std::string(s)); }

    /// Returns whether there is some information parsed.
    operator bool() const
       { return !libraryName.empty() || !mangledFunctionName.empty(); }
    /// Returns whether no information was parsed out of the original.
    bool operator! () const
       { return libraryName.empty() && mangledFunctionName.empty(); }

    /// Returns whether the translation was complete (offset is optional!).
    bool ParseString(std::string const& s);

    /// Returns the function name (mangled if nothing better).
    std::string const& function() const
      { return functionName.empty()? mangledFunctionName: functionName; }

    /// Returns only the library name (with suffix).
    std::string shortLibrary() const
      {
        size_t sep = libraryName.rfind('/');
        return (sep == std::string::npos)
          ? libraryName: libraryName.substr(sep + 1);
      //  return std::experimental::filesystem::path(libraryName).filename();
      }

    std::string original;            ///< String from the backtrace, unparsed.
    std::string libraryName;         ///< Parsed library name.
    std::string functionName;        ///< Parsed function name, demangled.
    std::string mangledFunctionName; ///< Parsed function name, unprocessed.
    void* address = nullptr;         ///< Function address.
    std::ptrdiff_t offset = 0;       ///< Instruction pointer offset.

       private:

    /// Returns whether the range is empty or invalid.
    static bool emptyRange(range_t const& r) { return r.first >= r.second; }

    /// Translates a range into a string.
    static std::string extract(std::string const& s, range_t const& r)
      { return emptyRange(r)? "": s.substr(r.first, r.second - r.first); }

    /// Runs the demangler and stores the result.
    void demangleFunction()
      { functionName = cet::demangle_symbol(mangledFunctionName); }

    /// Fills the information from an original string and parsed ranges.
    void setAll(
      std::string const& s,
      range_t addressStr, range_t libraryStr,
      range_t functionStr, range_t offsetStr
      );

  }; // CallInfo_t


  //----------------------------------------------------------------------------
  /**
   * @brief Class handling the output of information in a CallInfo_t object.
   *
   * This class has a "default" print function (also replicated as a call
   * operator), and a set of options that can be tweaked to change the amount
   * of information and format to be printed.
   *
   */
  class CallInfoPrinter {
      public:
    /// Set of options for printing.
    struct opt {
      /// List of available options.
      enum option_t {
        address,       ///< Print the instruction pointer memory address.
        demangled,     ///< Use demangled function names when possible.
        library,       ///< Print the library name the function lives in.
        shortLibrary,  ///< Print a shorter library name (base name only).
        offset,        ///< Print the offset from the beginning of function.
        NOptions       ///< Number of available options.
      }; // option_t

      std::bitset<NOptions> options; ///< Value of current options.

      /// Set one option `o` to the specified set value (true by default).
      opt& set(option_t o, bool set = true)
        { options.set(o, set); return *this; }

      /// Returns whether the specified option is set.
      bool has(option_t o) const { return options.test(o); }

    }; // opt

    opt options; ///< Set of current options.

    /// Default constructor: use default options.
    CallInfoPrinter() { setDefaultOptions(); }

    /// Constructor: use specified options.
    CallInfoPrinter(opt opts): options(opts) {}

    /// Override all the options.
    void setOptions(opt opts) { options = opts; }

    /// Print the content of info into the stream out, using the current options.
    template <typename Stream>
    void print(Stream&& out, CallInfo_t const& info) const;

    /// Print the content of info into the stream out, using the current options
    template <typename Stream>
    void operator() (Stream&& out, CallInfo_t const& info) const
      { print(std::forward<Stream>(out), info); }

    /// Sets this object to use a set of default options
    void setDefaultOptions() { options = defaultOptions(); }

    /// Returns a set of default options
    static opt defaultOptions()
      {
        opt options;
        options.set(opt::demangled);
        options.set(opt::library);
        options.set(opt::shortLibrary);
        options.set(opt::address);
        return options;
      }

  }; // CallInfoPrinter


  //----------------------------------------------------------------------------
  /// Helper operator to insert a call information in a stream with default options.
  template <typename Stream>
  inline Stream& operator<< (Stream&& out, CallInfo_t const& info);


  //----------------------------------------------------------------------------
  /// Backtrace printing options
  struct BacktracePrintOptions {

    unsigned int maxLines = 5; ///< Total number of lines to print.
    unsigned int skipLines = 1; ///< Number of lines to skip.
    bool countOthers = true; ///< Whether to print number of omitted lines.
    std::string indent; ///< Indentation string for all lines.
    std::string firstIndent; ///< Special indentation for the first line.

    /// Options for each single backtrace call information line.
    CallInfoPrinter::opt callInfoOptions = CallInfoPrinter::defaultOptions();

    /// Sets all indentation to the same specified `uniformIndent` string.
    void setUniformIndent(std::string uniformIndent)
      { indent = firstIndent = uniformIndent; }

  }; // struct BacktracePrintOptions


  /**
   * @brief Prints the full backtrace into a stream.
   * @tparam Stream type of output stream
   * @param out the output stream to insert output into
   * @param options printing options (see BacktracePrintOptions)
   *
   */
  template <typename Stream>
  void printBacktrace(Stream&& out, BacktracePrintOptions options);

  /**
   * @brief Prints the full backtrace into a stream with default options.
   * @tparam Stream type of output stream
   * @param out the output stream to insert output into
   */
  template <typename Stream>
  void printBacktrace(Stream&& out)
    { printBacktrace(std::forward<Stream>(out), BacktracePrintOptions()); }

  /**
   * @brief Prints the full backtrace into a stream.
   * @tparam Stream type of output stream
   * @param out the output stream to insert output into
   * @param maxLines print at most this many lines in the output (default: 5)
   * @param indent prepend a string in front of any new line (default: "  ")
   * @param callInfoOptions use these output options (default ones if null)
   *
   * The call information output options are described in
   * `CallInfoPrinter::opt` structure.
   *
   */
  template <typename Stream>
  void printBacktrace(
    Stream&& out,
    unsigned int maxLines, std::string indent = "  ",
    CallInfoPrinter::opt const* callInfoOptions = nullptr
    );

  //----------------------------------------------------------------------------
  /**
   * @brief Class triggering a `static_assert` failure.
   * @tparam T type accompanying the assertion
   * @tparam Enable assertion will fail only if `Enable` expands to `true`
   * @ingroup MetaprogrammingBase
   *
   * Instantiating this class anywhere (where it's legit) will trigger a static
   * assertion failure. Since the error message emitted by the compiler usually
   * contains an expansion of the template parameters, it is then possible to
   * see the "value" of type `T` that was used when the assertion failed.
   * The argument `Enable` allows to tune when the assertion should fail.
   *
   * For the following example, we want to investigate the value of the type
   * `element_type`, which is provided, among others, by `std::unique_ptr`.
   * We want to find out the exact type `element_type` of the collection type
   * passed to `OurClass`, but only when the collection type is, say, not
   * constant:
   * 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename Coll>
   * struct OurClass {
   *
   *   using Collection_t = Coll;
   *
   *   using value_type = typename Collection_t::element_type;
   *
   *   // DEBUG: have the compiler print `value_type`
   *   lar::debug::static_assert_on
   *     <value_type, std::is_const_v<std::remove_reference_t<Coll>>>
   *     debugVar;
   *
   * }; // struct OurClass
   *
   *
   * // this should never trigger a static assertion failure:
   * OurClass<std::unique_ptr<double>> doubleData;
   *
   * // this triggers a static assertion failure:
   * OurClass<std::unique_ptr<double[4]> const fourVectorData;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * (a working example is provided in `DebugUtils_test.h`).
   * The output with GCC 7.2 is similar to the following:
   * 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In file included from larcorealg/test/CoreUtils/DebugUtils_test.cc:17:0:
   * larcorealg/larcorealg/CoreUtils/DebugUtils.h: In instantiation of ‘struct lar::debug::details::THE_TYPE_IS<int [10]>’:
   * larcorealg/larcorealg/CoreUtils/DebugUtils.h:476:29:   required from ‘struct lar::debug::static_assert_on<int [10], true>’
   * larcorealg/test/CoreUtils/DebugUtils_test.cc:49:5:   required from ‘struct OurClass<const std::unique_ptr<int [10]> >’
   * larcorealg/test/CoreUtils/DebugUtils_test.cc:61:51:   required from here
   * larcorealg/larcorealg/CoreUtils/DebugUtils.h:467:7: error: static assertion failed: static_assert_on<T>: check the error message ("THE_TYPE_IS<>") for expansion of type `T`.
   *        static_assert(::util::always_false_v<T>,
   *        ^~~~~~~~~~~~~
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * The message of the assertion points to the key string ("THE_TYPE_IS"), and
   * it can be seen in the second line of this excerpt that the information is
   * printed as `struct lar::debug::details::THE_TYPE_IS<int [10]>`.
   * This is Clang 5.0:
   * 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In file included from larcorealg/test/CoreUtils/DebugUtils_test.cc:17:
   * larcorealg/larcorealg/CoreUtils/DebugUtils.h:451:7: error: static_assert failed "static_assert_on<T>: check the error message (\"THE_TYPE_IS<>\") for expansion of type `T`."
   *       static_assert(::util::always_false_v<T>,
   *       ^             ~~~~~~~~~~~~~~~~~~~~~~~~~
   * larcorealg/larcorealg/CoreUtils/DebugUtils.h:460:29: note: in instantiation of template class 'lar::debug::details::THE_TYPE_IS<int [10]>' requested here
   *     details::THE_TYPE_IS<T> _;
   *                             ^
   * larcorealg/test/CoreUtils/DebugUtils_test.cc:49:5: note: in instantiation of template class 'lar::debug::static_assert_on<int [10], true>' requested here
   *     debugVar;
   *     ^
   * larcorealg/test/CoreUtils/DebugUtils_test.cc:61:10: note: in instantiation of template class 'OurClass<const std::__1::unique_ptr<int [10], std::__1::default_delete<int [10]> > >' requested here
   *   (void) OurClass<std::unique_ptr<int[10]> const>();
   *          ^
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * where the type can be read in the message of the first note.
   */
  template <typename T, bool Enable /* = true */>
  struct static_assert_on;


  //----------------------------------------------------------------------------


} // namespace lar::debug


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace lar::debug {

  //----------------------------------------------------------------------------
  template <typename T>
  inline std::string demangle(T const* /* = nullptr */)
    { return cet::demangle_symbol(typeid(std::decay_t<T>).name()); }


  //----------------------------------------------------------------------------
  template <typename Stream>
  void CallInfoPrinter::print(Stream&& out, CallInfo_t const& info) const {
    if (!info) {
      out << info.original << " (?)";
      return;
    }

    if (info.mangledFunctionName.empty()) {
      if (options.has(opt::library)) {
        out << "in "
          << (options.has(opt::shortLibrary)? info.shortLibrary(): info.libraryName);
      }
      else out << "unknown";
    }
    else {
    //  auto flags = out.flags(); // not available in message facility streams
      out << info.function();
      auto offset = info.offset;
      if (offset && options.has(opt::offset)) {
        out << " [";
        if (offset > 0) out << '+';
        else {
          out << '-';
          offset = -offset;
        }
        out << std::showbase << std::hex << offset << "]";
      } // if has offset
    //  out << std::ios::setiosflags(flags);
      if (!info.libraryName.empty() && options.has(opt::library)) {
        out << " in "
          << (options.has(opt::shortLibrary)? info.shortLibrary(): info.libraryName);
      }
    }
    if (info.address && options.has(opt::address))
       out << " at " << ((void*) info.address);
    out << std::flush;
  } // CallInfoPrinter::print()


  //----------------------------------------------------------------------------
  template <typename Stream>
  inline Stream& operator<< (Stream&& out, CallInfo_t const& info) {
    CallInfoPrinter print;
    print(std::forward<Stream>(out), info);
    return out;
  }


  //----------------------------------------------------------------------------
  template <typename Stream>
  void printBacktrace(Stream&& out, BacktracePrintOptions options) {
    unsigned int nSkip = std::max(options.skipLines, 0U);
    std::vector<void*> buffer
      (nSkip + std::max(options.maxLines, 200U), nullptr);

    unsigned int const nItems
      = (unsigned int) backtrace(buffer.data(), buffer.size());

    // convert the calls in the buffer into a vector of strings
    char** symbols = backtrace_symbols(buffer.data(), buffer.size());
    if (!symbols) {
      out << options.firstIndent << "<failed to get the call stack>\n"
        << std::flush;
      return;
    }
    std::vector<CallInfo_t> callStack;
    for (size_t i = 0; i < buffer.size(); ++i)
      callStack.push_back((const char*) symbols[i]);
    std::free(symbols);

    size_t lastItem = nSkip + options.maxLines;
    if (lastItem > nItems) lastItem = nItems;
    if (lastItem >= buffer.size()) --lastItem;

    CallInfoPrinter print(options.callInfoOptions);
    for (size_t i = nSkip; i < lastItem; ++i) {
      out << (i == 0? options.firstIndent: options.indent);
      print(std::forward<Stream>(out), callStack[i]);
      out << "\n";
    }
    if ((lastItem < nItems) && options.countOthers) {
      out << options.indent << " ... and other " << (nItems - lastItem);
      if (nItems == buffer.size()) out << " (or more)";
      out << " levels\n";
    }
    out << std::flush;

  } // printBacktrace()


  //----------------------------------------------------------------------------
  template <typename Stream>
  void printBacktrace(
    Stream&& out,
    unsigned int maxLines, std::string indent /* = "  " */,
    CallInfoPrinter::opt const* callInfoOptions /* = nullptr */
    )
  {
    BacktracePrintOptions options;
    options.maxLines = maxLines;
    options.indent = options.firstIndent = indent;
    if (callInfoOptions) options.callInfoOptions = *callInfoOptions;
    printBacktrace(std::forward<Stream>(out), options);
  }


  //----------------------------------------------------------------------------
  namespace details {

    template <typename T>
    struct THE_TYPE_IS {
      // if the assertion condition didn't depend on a template parameters,
      // the assertion failure error would *always* be triggered
      static_assert(::util::always_false_v<T>,
        "static_assert_on<T>: check the error message (\"THE_TYPE_IS<>\") for expansion of type `T`."
        );
    }; // THE_TYPE_IS<>

  } // namespace details

  template <typename T, bool Enable = true>
  struct static_assert_on {
    static_assert_on() = default;
    static_assert_on(T) {}
    details::THE_TYPE_IS<T> _;
  }; // static_assert_on<>

  template <typename T>
  struct static_assert_on<T, false> {};


  //----------------------------------------------------------------------------

} // namespace lar::debug


#endif // LARCOREALG_COREUTILS_DEBUGUTILS_H
