# ======================================================================
#  larcorealg main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larcorealg/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(larcorealg VERSION 10.00.04 LANGUAGES CXX)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

# Force C++20. LArSoft v10_21_01 / e26 defaults to C++17 (CETPKG_CXX_STANDARD=17,
# buildtool passes -DCMAKE_CXX_STANDARD=17), but the range-v3 removal in
# Geometry/Iterable.h uses std::ranges, which is only available in C++20.
# cetmodules delegates the standard to CMAKE_CXX_STANDARD (cet_set_compiler_flags
# has no CXX_STANDARD option); this normal variable overrides the cache value.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  # -Wno-error=cpp: the cvmfs ROOT (v6_28_12) is built for C++17 and its headers
  # emit a #warning when consumed under a different standard; WERROR would turn
  # that into an error. (Diagnostic only -- C++20 code against a C++17 ROOT is not
  # a supported configuration.)
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-error=cpp
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(canvas_root_io REQUIRED EXPORT)
find_package(canvas REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)
find_package(cetlib REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)

find_package(Boost COMPONENTS unit_test_framework filesystem REQUIRED EXPORT)
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
