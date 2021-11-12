# ======================================================================
#  larcorealg main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larcorealg/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(larcorealg VERSION 09.03.01.01 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_boost()
find_ups_root()

find_ups_product( canvas )
find_ups_product( cetlib )
find_ups_product( cetlib_except )
find_ups_product(messagefacility)
find_ups_product( larcoreobj )
find_ups_product(clhep)

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(CetMake)
include(BasicPlugin)

# source
add_subdirectory(larcorealg)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
