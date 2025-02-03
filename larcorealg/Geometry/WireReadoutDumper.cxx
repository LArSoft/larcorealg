/**
 * @file    larcorealg/Geometry/WireReadoutDumper.cxx
 * @brief   Algorithm dumping the full detector geometry down to wires.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    January 31, 2025
 * @file    larcorealg/Geometry/WireReadoutDumper.h
 */

// library header
#include "larcorealg/Geometry/WireReadoutDumper.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

// C++ libraries
#include <ostream>

// -----------------------------------------------------------------------------
geo::WireReadoutDumper::WireReadoutDumper(geo::GeometryCore const* geom,
                                          geo::WireReadoutGeom const* wireGeom,
                                          geo::AuxDetGeometryCore const* auxGeom)
  : fGeom{geom}, fWireGeom{wireGeom}, fAuxGeom{auxGeom}
{}

// -----------------------------------------------------------------------------
void geo::WireReadoutDumper::dump(std::ostream& out,
                                  std::string const& indent,
                                  std::string const& firstIndent) const
{

  dumpDetectorInfo(out, indent, firstIndent);

  if (fGeom) {
    for (auto const& cryostat : fGeom->Iterate<geo::CryostatGeo>()) {
      out << "\n";
      dumpCryostat(out, cryostat, indent + "    ", indent + "  ");
    }
  }

  if (fAuxGeom) {
    out << "\n";
    dumpAuxiliaryDetectors(out, indent, indent);
  }

  out << '\n';

} // geo::WireReadoutDumper::dump()

// -----------------------------------------------------------------------------
geo::WireReadoutDumper::StreamAdapter geo::WireReadoutDumper::toStream(
  std::string indent,
  std::string firstIndent) const
{
  return {this, std::move(indent), std::move(firstIndent)};
}

// -----------------------------------------------------------------------------
void geo::WireReadoutDumper::dumpDetectorInfo(std::ostream& out,
                                              std::string const& indent,
                                              std::string const& firstIndent) const
{

  out << firstIndent << "Detector " << (fGeom ? fGeom->DetectorName() : "") << " has "
      << (fGeom ? std::to_string(fGeom->Ncryostats()) : "unknown") << " cryostats and "
      << (fAuxGeom ? std::to_string(fAuxGeom->NAuxDets()) : "unknown") << " auxiliary detectors:";

  if (fGeom) {
    geo::BoxBoundedGeo const& box = fGeom->DetectorEnclosureBox();

    out << "\n"
        << indent << "  Detector enclosure: " << box.Min() << " -- " << box.Max() << " cm => ( "
        << box.SizeX() << " x " << box.SizeY() << " x " << box.SizeZ() << " ) cm^3";
  }

} // geo::WireReadoutDumper::dumpDetectorInfo()

// -----------------------------------------------------------------------------
void geo::WireReadoutDumper::dumpCryostat(std::ostream& out,
                                          geo::CryostatGeo const& cryostat,
                                          std::string const& indent,
                                          std::string const& firstIndent) const
{

  // the indentation game:
  //   `PrintXxxxInfo()` methods don't indent the first line;
  //   our `dumpCryostat()` do as directed, but don't start with a new line

  out << firstIndent;
  cryostat.PrintCryostatInfo(out, indent, cryostat.MaxVerbosity);

  //
  // TPC
  //
  for (geo::TPCGeo const& tpc : cryostat.IterateTPCs()) {

    out << "\n" << indent;
    tpc.PrintTPCInfo(out, indent + "  ", tpc.MaxVerbosity);

    // for (unsigned int const planeNo: tpc.Nplanes()) {
    //   geo::PlaneGeo const& plane = fWireGeom->Plane({ tpc.ID(), planeNo });

    for (auto const& plane : fWireGeom->Iterate<geo::PlaneGeo>(tpc.ID())) {
      out << "\n";
      dumpTPCplane(out, plane, indent + "    ", indent + "  ");

    } // for plane
  }   // for TPC

  //
  // optical detectors
  //
  for (unsigned int const iOpDet : util::counter(cryostat.NOpDet())) {
    geo::OpDetGeo const& opDet = cryostat.OpDet(iOpDet);
    out << "\n" << indent << "  [OpDet #" << iOpDet << "] ";
    opDet.PrintOpDetInfo(out, indent + "    ", opDet.MaxVerbosity);
  } // for

} // geo::WireReadoutDumper::dumpCryostat()

// -----------------------------------------------------------------------------
void geo::WireReadoutDumper::dumpTPCplane(std::ostream& out,
                                          geo::PlaneGeo const& plane,
                                          std::string const& indent,
                                          std::string const& firstIndent) const
{

  const unsigned int nWires = plane.Nwires();

  out << firstIndent;
  plane.PrintPlaneInfo(out, indent, plane.MaxVerbosity);
  SigType_t const sigType = fWireGeom->SignalType(plane.ID());
  out << "\n"
      << indent << "signal type: " << SignalTypeName(sigType) << " (" << static_cast<int>(sigType)
      << ")";

  for (unsigned int const wireNo : util::counter(nWires)) {
    geo::WireID const wireID{plane.ID(), wireNo};
    geo::WireGeo const& wire = plane.Wire(wireID);
    out << "\n" << indent << wireID << " ";
    wire.PrintWireInfo(out, indent, wire.MaxVerbosity);
  } // for wire

} // geo::WireReadoutDumper::dumpTPCplane()

// -----------------------------------------------------------------------------
void geo::WireReadoutDumper::dumpAuxiliaryDetectors(std::ostream& out,
                                                    std::string const& indent,
                                                    std::string const& firstIndent) const
{

  unsigned int const nAuxDets = fAuxGeom->NAuxDets();
  out << firstIndent << "Auxiliary detectors (" << nAuxDets << "):";
  for (unsigned int const iDet : util::counter(nAuxDets)) {
    geo::AuxDetGeo const& auxDet = fAuxGeom->AuxDet(iDet);

    out << '\n' << indent << "[#" << iDet << "] ";
    auxDet.PrintAuxDetInfo(out, indent + "  ", auxDet.MaxVerbosity);

    unsigned int const nSensitive = auxDet.NSensitiveVolume();
    switch (nSensitive) {
    case 0: break;
    case 1: {
      AuxDetSensitiveGeo const& auxDetS = auxDet.SensitiveVolume(0U);
      out << "\n" << indent << "  ";
      auxDetS.PrintAuxDetInfo(out, indent + "    ", auxDetS.MaxVerbosity);
      break;
    }
    default:
      for (unsigned int const iSens : util::counter(nSensitive)) {
        out << "\n" << indent << "  [#" << iSens << "] ";
        AuxDetSensitiveGeo const& auxDetS = auxDet.SensitiveVolume(iSens);
        auxDetS.PrintAuxDetInfo(out, indent + "    ", auxDetS.MaxVerbosity);
      } // for
      break;
    } // if sensitive detectors

  } // for auxiliary detector

} // geo::WireReadoutDumper::dumpAuxiliaryDetectors()

// -----------------------------------------------------------------------------
std::ostream& geo::operator<<(std::ostream& out,
                              geo::WireReadoutDumper::StreamAdapter const& dumpAdapter)
{
  dumpAdapter.dumper->dump(out, dumpAdapter.indent, dumpAdapter.firstIndent);
  return out;
}

// -----------------------------------------------------------------------------
