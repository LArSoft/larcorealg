#include "larcorealg/CoreUtils/SearchPathPlusRelative.h"

#include "cetlib/search_path.h"

std::string lar::searchPathPlusRelative(std::string relativePath, std::string fileName)
{
  // Add a final directory separator ("/") to relPath if not already there.
  if (!relativePath.empty() && (relativePath.back() != '/')) relativePath += '/';

  relativePath += std::move(fileName);

  // Search all reasonable locations for the file; we see if by any chance FW_SEARCH_PATH
  // directory is set and try there.
  cet::search_path const sp("FW_SEARCH_PATH");
  return sp.find_file(relativePath);
}
