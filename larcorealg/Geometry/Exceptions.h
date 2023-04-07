/**
 * @file   larcorealg/Geometry/Exceptions.h
 * @brief  Collection of exceptions for Geometry system
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 24, 2016
 *
 * This is currently a header-only library.
 *
 * It offers:
 * - InvalidWireError (for bad wire numbers)
 *
 */

#ifndef LARCOREALG_GEOMETRY_EXCEPTIONS_H
#define LARCOREALG_GEOMETRY_EXCEPTIONS_H

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::PlaneID, ...

// Framework libraries
#include "cetlib_except/exception.h"

// C/C++ standard library
#include <limits> // std::numeric_limits<>

namespace geo {

  /** **************************************************************************
   * @brief Exception thrown on invalid wire number
   *
   * This class is thrown, e.g., by Geometry::NearestWireID().
   * It contains the erroneous wire number, a suggestion of which one might be
   * the right one, and a plane ID.
   *
   * The wire numbers are signed.
   *
   *
   */
  class InvalidWireError : public cet::exception {
  public:
    /// Value used to represent an invalid wire number
    static constexpr int InvalidWireNo = std::numeric_limits<int>::max();

    /// Constructor with the complete information
    InvalidWireError(std::string cat, PlaneID const& planeID, int badWireNo, int betterWireNo)
      : cet::exception(cat), fPlaneID(planeID), fWireNumber(badWireNo), fBetterWireNo(betterWireNo)
    {}

    /// Constructor: no wire suggestions
    InvalidWireError(std::string cat, PlaneID const& planeID, int badWireNo)
      : cet::exception(cat), fPlaneID(planeID), fWireNumber(badWireNo)
    {}

    /// Constructor: no plane information
    InvalidWireError(std::string cat, int badWireNo, int betterWireNo)
      : cet::exception(cat), fWireNumber(badWireNo), fBetterWireNo(betterWireNo)
    {}

    /// Constructor: no plane information and no suggestion
    InvalidWireError(std::string cat, int badWireNo) : cet::exception(cat), fWireNumber(badWireNo)
    {}

    /// @{
    /// @name Access to bad wire

    /// Returns whether we known the bad wire number
    bool hasBadWire() const { return fWireNumber != InvalidWireNo; }

    /// Returns the bad wire number
    int badWire() const { return fWireNumber; }

    /// Returns the bad wire ID
    WireID badWireID() const { return makeWireID(fWireNumber); }

    /// @}

    /// @{
    /// @name Access to suggested wire

    /// Returns whether we known a better wire number
    bool hasSuggestedWire() const { return fBetterWireNo != InvalidWireNo; }

    /// Returns a better wire number
    int suggestedWire() const { return fBetterWireNo; }

    /// Returns a better wire ID
    WireID suggestedWireID() const { return makeWireID(fBetterWireNo); }

    /// @}

    /// @{
    /// @name Plane access

    /// Return whether a plane is recorded with the exception
    bool hasPlane() const { return fPlaneID.isValid; }

    /// Return the plane ID recorded with the exception
    PlaneID const& planeID() const { return fPlaneID; }

    /// @}

  private:
    PlaneID fPlaneID; ///< plane the wire belongs to

    /// the invalid wire number
    int fWireNumber = InvalidWireNo;

    /// a suggestion for a good wire number
    int fBetterWireNo = InvalidWireNo;

    /// Transform a wire number into wire ID
    WireID makeWireID(int wireNo) const { return {fPlaneID, (WireID::WireID_t)wireNo}; }

  }; // class InvalidWireError

} // namespace geo

#endif // LARCOREALG_GEOMETRY_EXCEPTIONS_H
