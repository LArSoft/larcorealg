/**
 * @file    larcorealg/Geometry/WireReadoutDumper.h
 * @brief   Algorithm dumping the full detector geometry down to wires.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    January 31, 2025
 * @file    larcorealg/Geometry/WireReadoutDumper.cxx
 */

#ifndef LARCOREALG_GEOMETRY_WIREREADOUTDUMPER_H
#define LARCOREALG_GEOMETRY_WIREREADOUTDUMPER_H

// LArSoft libraries
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

// C++ libraries
#include <iosfwd> // std::ostream
#include <string>
#include <utility> // std::move()

// -----------------------------------------------------------------------------
namespace geo {
  class WireReadoutDumper;
}
/**
 * @brief Algorithm dumping the full detector geometry down to wires.
 *
 * This algorithm dumps all the available information from the geometry,
 * including the chamber information (from `geo::GeometryCore`), the wire and
 * readout information (from `geo::WireReadoutGeom`) and the auxiliary detectors
 * (from `geo::AuxDetGeometryCore`).
 * Any of the service providers can be omitted (`nullptr`), in which case the
 * relative information will not be printed.
 *
 * A typical usage example is to create the algorithm and then call its `dump()`
 * method, like in:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto const* geom = lar::providerFrom<geo::Geometry>();
 * auto const* wireGeom = &(art::ServiceHandle<geo::WireReadout>()->Get());
 * auto const* auxDetGeom = art::ServiceHandle<geo::AuxDetGeometry>()->GetProviderPtr();
 * geo::WireReadoutDumper const dumper{ geom, wireGeom, auxDetGeom };
 * dumper.dump(std::cout);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * The interface supports only dumping to C++ output streams (including
 * `std::ostringstream` to bridge to libraries which can us strings).
 * However, for streams which support inserter operators to `std::ostream`,
 * the following is also possible:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * mf::LogInfo{ "DumpGeometry" } << dumper.toStream();
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * The messagefacility library streams (like in this example) fall into this
 * category.
 *
 * The output is documented in the `dump()` method.
 *
 */
class geo::WireReadoutDumper {

  /// Cached geometry service provider.
  geo::GeometryCore const* fGeom = nullptr;

  /// Cached wire readout service provider.
  geo::WireReadoutGeom const* fWireGeom = nullptr;

  /// Cached auxiliary detector geometry service provider.
  geo::AuxDetGeometryCore const* fAuxGeom = nullptr;

  /// Helper for insertion into C++ output streams.
  struct StreamAdapter {
    WireReadoutDumper const* dumper;
    std::string indent;
    std::string firstIndent;
  };

  friend std::ostream& operator<<(std::ostream& out, StreamAdapter const& dumpAdapter);

public:
  /// Constructor: acquires geometry service providers.
  WireReadoutDumper(geo::GeometryCore const* geom,
                    geo::WireReadoutGeom const* wireGeom,
                    geo::AuxDetGeometryCore const* auxGeom);

  /**
   * @brief Dumps the full geometry information into a stream.
   * @param out the stream to write the information in
   * @param indent string to be prepended to each output line except the first
   * @param firstIndent string to be prepended to the first line only
   *
   * The output is in human-readable format, structured as follows:
   *  * all the information for each cryostat, including:
   *     * all the information for each TPC, including:
   *         * all the information for each readout plane, including:
   *             * all the information for each sensitive element (e.g. wire)
   *     * all the optical detector information
   *  * all information for the "auxiliary" detectors (typically, scintillators
   *    for cosmic ray tagging), including:
   *     * all the information for the "sensitive detectors" in each module.
   *
   * The output starts on the current line of the stream and it ends with a new
   * line.
   *
   */
  void dump(std::ostream& out, std::string const& indent, std::string const& firstIndent) const;

  /**
   * @brief Dumps the full geometry information into a stream.
   * @param out the stream to write the information in
   * @param indent (default: none) string to be prepended to all output lines
   * @see `dump(std::ostream&, std::string const&, std::string const&) const`
   *
   * Dumps the geometry like
   * `dump(std::ostream&, std::string const&, std::string const&) const` does,
   * but using the same indentation string for all lines.
   */
  void dump(std::ostream& out, std::string const& indent = "") const { dump(out, indent, indent); }

  /**
   * @brief Adapter to send the dump to a C++ output stream via insertion.
   * @param indent (default: none) string to be prepended to all output lines
   * @param firstIndent string to be prepended to the first line only
   * @return an opaque adapter for insertion into a C++ output stream
   * @see `dump()`
   *
   * Example: after a `geo::WireReadoutDumper dumper` has been initialized,
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << dumper.toStream("> ");
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will dump the full geometry into `std::cout`.
   *
   */
  StreamAdapter toStream(std::string indent, std::string firstIndent) const;

  /**
   * @brief Adapter to send the dump to a C++ output stream via insertion.
   * @param indent (default: none) string to be prepended to all output lines
   * @return an opaque adapter for insertion into a C++ output stream
   * @see `dump()`
   *
   * Example: after a `geo::WireReadoutDumper dumper` has been initialized,
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << dumper.toStream("> ");
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will dump the full geometry into `std::cout`.
   *
   */
  auto toStream(std::string indent = "") const { return toStream(indent, std::move(indent)); }

private:
  /// Dumps the general detector information.
  void dumpDetectorInfo(std::ostream& out,
                        std::string const& indent,
                        std::string const& firstIndent) const;

  /// Dumps all content of the specified `cryostat`.
  void dumpCryostat(std::ostream& out,
                    geo::CryostatGeo const& cryostat,
                    std::string const& indent,
                    std::string const& firstIndent) const;

  /// Dumps the information from the specified `plane` (including wires).
  void dumpTPCplane(std::ostream& out,
                    geo::PlaneGeo const& plane,
                    std::string const& indent,
                    std::string const& firstIndent) const;

  /// Dumps all auxiliary detector information.
  void dumpAuxiliaryDetectors(std::ostream& out,
                              std::string const& indent,
                              std::string const& firstIndent) const;

}; // geo::WireReadoutDumper

// -----------------------------------------------------------------------------
namespace geo {
  /// Helper for `geo::WireReadoutDumper::toStream()`.
  std::ostream& operator<<(std::ostream& out, WireReadoutDumper::StreamAdapter const& dumpAdapter);
}

// -----------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_WIREREADOUTDUMPER_H
