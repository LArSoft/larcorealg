# ======================================================================
#  larcorealg main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larcorealg/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(larcorealg VERSION 09.13.01 LANGUAGES CXX)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(canvas_root_io REQUIRED EXPORT)
find_package(canvas REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)
find_package(cetlib REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)

find_package(Boost COMPONENTS headers unit_test_framework REQUIRED EXPORT)
find_package(CLHEP COMPONENTS Geometry Vector REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core GenVector Geom MathCore Matrix Physics REQUIRED EXPORT)

find_package(larcoreobj REQUIRED EXPORT)

include(BuildDictionary)

# source
add_subdirectory(larcorealg)

# tests
add_subdirectory(test)

# ups - table and config files

# packaging utility
cet_cmake_config()
