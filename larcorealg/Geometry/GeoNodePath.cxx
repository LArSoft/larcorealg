/**
 * @file   larcorealg/Geometry/GeoNodePath.cxx
 * @brief  Class representing a path in ROOT geometry.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/Geometry/GeoNodePath.h`
 */

// LArSoft libraries
#include "larcorealg/Geometry/GeoNodePath.h"

// ROOT libraries
#include "TGeoNode.h"

geo::GeoNodePath::operator std::string() const
{
  std::string s = "[";
  auto it = fNodes.cbegin(), end = fNodes.cend();
  if (it != end) {
    s += (*it++)->GetName();
    while (++it != fNodes.cend()) {
      s += '/';
      s += (*it)->GetName();
    }
  } // if
  return s + "]";
}
