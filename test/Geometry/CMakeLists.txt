cet_enable_asserts()
cet_test_env_prepend(FHICL_FILE_PATH ${PROJECT_BINARY_DIR}/fcl ${PROJECT_BINARY_DIR}/job)
cet_test_env(FW_SEARCH_PATH=${PROJECT_BINARY_DIR}/gdml)

cet_test(geo_vectors_utils_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::geo_vectors_utils_TVector
  larcorealg::geo_vectors_utils
)

cet_test(LineClosestPoint_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::LineClosestPoint
  larcoreobj::geo_vectors
)

# test libraries
set(GeometryTestLib_SOURCES
  GeometryTestAlg.cxx
  GeometryGeoIDTestAlg.cxx
  GeometryIteratorTestAlg.cxx
  GeometryIteratorLoopTestAlg.cxx
  ChannelMapStandardTestAlg.cxx
)

list(TRANSFORM GeometryTestLib_SOURCES REPLACE "\\.cxx$" ".h"
  OUTPUT_VARIABLE GeometryTestLib_HEADERS)

cet_make_library(LIBRARY_NAME GeometryTestLib USE_PROJECT_NAME USE_BOOST_UNIT
  SOURCE ${GeometryTestLib_SOURCES}
  LIBRARIES
  PUBLIC
  larcorealg::TestUtils
  larcoreobj::geo_vectors
  PRIVATE
  larcorealg::Exceptions
  larcorealg::Geometry
  larcorealg::geo
  larcorealg::geo_vectors_utils
  larcoreobj::SimpleTypesAndConstants
  messagefacility::MF_MessageLogger
  cetlib::cetlib
  cetlib_except::cetlib_except
  ROOT::Core
)

# geometry test on "standard" geometry

# unit test (use geometry configuration from config file, the rest is ignored)
# this uses BOOST for the test
cet_test(geometry_iterator_test USE_BOOST_UNIT
  SOURCE geometry_iterator_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS -- ./test_geometry.fcl
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcorealg::GeometryTestLib
  larcorealg::geometry_unit_test_base
)

# unit test (use geometry configuration from config file, the rest is ignored)
# this uses BOOST for the test
cet_test(geometry_geoid_test USE_BOOST_UNIT
  SOURCE geometry_geoid_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS -- ./test_geometry.fcl
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcorealg::GeometryTestLib
  larcorealg::geometry_unit_test_base
)

# unit test (use the hard-coded channel mapping for "standard" LArTPCdetector)
cet_test(geometry_test
  SOURCE geometry_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS ./test_geometry.fcl
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcorealg::GeometryTestLib
  larcorealg::geometry_unit_test_base
)

# test of standalone geometry loading (use the hard-coded channel mapping for "standard" LArTPCdetector)
cet_test(geometry_loader_test
  SOURCE geometry_loader_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS ./test_geometry.fcl
  LIBRARIES PRIVATE
  larcorealg::GeometryTestLib
  larcorealg::Geometry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

# test to verify loops on geometry elements by geometry iterators (hard-coded channel mapping for "standard" LArTPCdetector)
cet_test(geometry_iterator_loop_test
  SOURCE geometry_iterator_loop_test.cxx
  DATAFILES test_geometry_iterator_loop.fcl
  TEST_ARGS ./test_geometry_iterator_loop.fcl
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcorealg::GeometryTestLib
  larcorealg::geometry_unit_test_base
  messagefacility::MF_MessageLogger
)


# test for standard channel mapping (BOOST unit test)
cet_test(geometry_standardchannelmapping_test USE_BOOST_UNIT
  SOURCE geometry_standardchannelmapping_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS ./test_geometry.fcl
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcorealg::GeometryTestLib
  larcorealg::geometry_unit_test_base
)


# simple unit test
# this uses BOOST for the test
cet_test(geometry_thirdplaneslope_test USE_BOOST_UNIT
  SOURCE geometry_thirdplaneslope_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS -- ./test_geometry.fcl
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcorealg::GeometryTestLib
  larcorealg::geometry_unit_test_base
  larcoreobj::headers
)

# decomposition engine tests
cet_test(Decomposer_test USE_BOOST_UNIT
  SOURCE Decomposer_test.cxx
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcorealg::geo_vectors_utils_TVector
  ROOT::Physics
  ROOT::GenVector
)

cet_test(SimpleGeo_test USE_BOOST_UNIT
  SOURCE SimpleGeo_test.cxx
)

cet_test(driftvolumes_test SOURCE driftvolumes_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS ./test_geometry.fcl
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcorealg::geometry_unit_test_base
  messagefacility::MF_MessageLogger
  TEST_PROPERTIES ENVIRONMENT FW_SEARCH_PATH=${PROJECT_BINARY_DIR}/gdml
)

cet_test(geometrydatacontainers_test USE_BOOST_UNIT
  SOURCE geometrydatacontainers_test.cxx
  LIBRARIES PRIVATE
  larcorealg::GeometryDataContainers
  larcorealg::CoreUtils
  larcoreobj::SimpleTypesAndConstants
)

cet_test(geoidmapper_test USE_BOOST_UNIT
  SOURCE geoidmapper_test.cxx
  LIBRARIES PRIVATE
  larcorealg::GeometryIDmapper
  larcoreobj::SimpleTypesAndConstants
)

cet_test(readoutdatacontainers_test USE_BOOST_UNIT
  SOURCE readoutdatacontainers_test.cxx
  LIBRARIES PRIVATE
  larcorealg::ReadoutDataContainers
  larcoreobj::SimpleTypesAndConstants
)

cet_test(topology_test USE_BOOST_UNIT
  SOURCE topology_test.cxx
  LIBRARIES PRIVATE
  larcorealg::Geometry
)

file(COPY ${GeometryTestLib_HEADERS}
  DESTINATION "${PROJECT_BINARY_DIR}/larcorealg/test/Geometry")
install_headers(SUBDIRNAME larcorealg)
install_fhicl()
install_source()
