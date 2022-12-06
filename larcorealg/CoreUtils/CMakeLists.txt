cet_make_library(SOURCE
  makeValueIndex.h
  DebugUtils.cxx
  NumericUtils.cxx
  SortByPointers.h
  LIBRARIES PUBLIC
  range-v3::range-v3
  cetlib_except::cetlib_except
)

cet_make_library(LIBRARY_NAME ParticleFilters INTERFACE
  SOURCE ParticleFilters.h
  LIBRARIES INTERFACE
  ROOT::Geom
  ROOT::Matrix
  ROOT::Physics
)

install_headers()
install_source()
