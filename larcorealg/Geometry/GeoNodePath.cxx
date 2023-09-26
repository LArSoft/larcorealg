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

// Boost includes
#include "boost/functional/hash.hpp"

// C++ STL libraries
#include <iostream>

namespace {
  std::size_t combine(std::size_t a, std::size_t const b)
  {
    boost::hash_combine(a, b);
    return a;
  }

  std::hash<TGeoNode const*> node_hash;
  geo::GeoNodePath::Entry entry_for(TGeoNode const* node) { return {node, node_hash(node)}; }
  geo::GeoNodePath::Entry entry_for(geo::GeoNodePath::Entry const& parent, TGeoNode const* node)
  {
    return {node, combine(parent.hash_value, node_hash(node))};
  }
}

namespace geo {

  GeoNodePath::GeoNodePath(TGeoNode const* topNode) : fNodes{entry_for(topNode)} {}

  bool GeoNodePath::empty() const { return fNodes.empty(); }

  GeoNodePath::Depth_t GeoNodePath::depth() const { return fNodes.size(); }

  TGeoNode const* GeoNodePath::current() const { return current_entry().node; }

  auto GeoNodePath::current_entry() const -> Entry { return fNodes.back(); }

  auto GeoNodePath::parent_entry() const -> Entry
  {
    if (auto const size = depth(); size > 1ull) { return fNodes[size - 2ull]; }
    return {};
  }

  void GeoNodePath::append(TGeoNode const* node)
  {
    auto entry = empty() ? entry_for(node) : entry_for(current_entry(), node);
    fNodes.push_back(std::move(entry));
  }

  void GeoNodePath::pop() { fNodes.pop_back(); }

  GeoNodePath::operator std::string() const
  {
    std::string s = "[";
    auto it = fNodes.cbegin(), end = fNodes.cend();
    if (it != end) {
      s += it->node->GetName();
      while (++it != fNodes.cend()) {
        s += '/';
        s += it->node->GetName();
      }
    } // if
    return s + ']';
  }
}
