cet_enable_asserts()

cet_test( geo_vectors_utils_test
  LIBRARIES
    ROOT::Core
    ROOT::Physics
    ROOT::GenVector
  USE_BOOST_UNIT
  )

# test libraries
cet_make_library( LIBRARY_NAME GeometryTestLib
                  SOURCE
                    GeometryTestAlg.cxx
                    GeometryGeoIDTestAlg.cxx
                    GeometryIteratorTestAlg.cxx
                    GeometryIteratorLoopTestAlg.cxx
                    ChannelMapStandardTestAlg.cxx
                  LIBRARIES larcorealg_Geometry
                        larcorealg_TestUtils
                        canvas
                        messagefacility::headers
                        fhiclcpp::fhiclcpp fhiclcpp::types
                        cetlib cetlib_except
                        Boost::unit_test_framework
          )


# geometry test on "standard" geometry

# unit test (use geometry configuration from config file, the rest is ignored)
# this uses BOOST for the test
cet_test(geometry_iterator_test
  SOURCES geometry_iterator_test.cxx
  LIBRARIES larcorealg_Geometry
            GeometryTestLib
            messagefacility::headers
            fhiclcpp::fhiclcpp fhiclcpp::types
            cetlib cetlib_except
  USE_BOOST_UNIT
  DATAFILES test_geometry.fcl
  TEST_ARGS -- ./test_geometry.fcl
)

# unit test (use geometry configuration from config file, the rest is ignored)
# this uses BOOST for the test
cet_test(geometry_geoid_test
  SOURCES geometry_geoid_test.cxx
  LIBRARIES larcorealg_Geometry
            GeometryTestLib
            messagefacility::headers
            fhiclcpp::fhiclcpp fhiclcpp::types
            cetlib cetlib_except
  USE_BOOST_UNIT
  DATAFILES test_geometry.fcl
  TEST_ARGS -- ./test_geometry.fcl
)

# unit test (use the hard-coded channel mapping for "standard" LArTPCdetector)
cet_test(geometry_test
  SOURCES geometry_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS ./test_geometry.fcl
  LIBRARIES larcorealg_Geometry
            GeometryTestLib
            messagefacility::headers
            fhiclcpp::fhiclcpp fhiclcpp::types
            cetlib cetlib_except
)

# test of standalone geometry loading (use the hard-coded channel mapping for "standard" LArTPCdetector)
cet_test(geometry_loader_test
  SOURCES geometry_loader_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS ./test_geometry.fcl
  LIBRARIES larcorealg_Geometry
            GeometryTestLib
            messagefacility::headers
            fhiclcpp::fhiclcpp fhiclcpp::types
            cetlib cetlib_except
)

# test to verify loops on geometry elements by geometry iterators (hard-coded channel mapping for "standard" LArTPCdetector)
cet_test(geometry_iterator_loop_test
  SOURCES geometry_iterator_loop_test.cxx
  DATAFILES test_geometry_iterator_loop.fcl
  TEST_ARGS ./test_geometry_iterator_loop.fcl
  LIBRARIES larcorealg_Geometry
            GeometryTestLib
            messagefacility::headers
            fhiclcpp::fhiclcpp fhiclcpp::types
            cetlib cetlib_except
)


# test for standard channel mapping (BOOST unit test)
cet_test(geometry_standardchannelmapping_test
  SOURCES geometry_standardchannelmapping_test.cxx
  LIBRARIES larcorealg_Geometry
            GeometryTestLib
            messagefacility::headers
            fhiclcpp::fhiclcpp fhiclcpp::types
            cetlib cetlib_except
  USE_BOOST_UNIT
  DATAFILES test_geometry.fcl
  TEST_ARGS ./test_geometry.fcl
)


# simple unit test
# this uses BOOST for the test
cet_test(geometry_thirdplaneslope_test
  SOURCES geometry_thirdplaneslope_test.cxx
  LIBRARIES larcorealg_Geometry
            messagefacility::headers
            fhiclcpp::fhiclcpp fhiclcpp::types
            cetlib cetlib_except
  USE_BOOST_UNIT
  DATAFILES test_geometry.fcl
  TEST_ARGS -- ./test_geometry.fcl
)


# decomposition engine tests
cet_test(Decomposer_test
  SOURCES Decomposer_test.cxx
  LIBRARIES
    ROOT::Physics
    ROOT::GenVector
    ROOT::Core # TVersionCheck
  USE_BOOST_UNIT
  )

cet_test(SimpleGeo_test
  SOURCES SimpleGeo_test.cxx
  USE_BOOST_UNIT
)


cet_test(driftvolumes_test
  SOURCES driftvolumes_test.cxx
  DATAFILES test_geometry.fcl
  TEST_ARGS ./test_geometry.fcl
  LIBRARIES larcorealg_Geometry
            messagefacility::headers
            fhiclcpp::fhiclcpp fhiclcpp::types
            cetlib cetlib_except
  USE_BOOST_UNIT
  )

cet_test(geometrydatacontainers_test
  SOURCES geometrydatacontainers_test.cxx
  USE_BOOST_UNIT
  )

cet_test(geoidmapper_test
  SOURCES geoidmapper_test.cxx
  USE_BOOST_UNIT
  )

cet_test(readoutdatacontainers_test
  SOURCES readoutdatacontainers_test.cxx
  USE_BOOST_UNIT
  )

cet_test(topology_test
  SOURCES topology_test.cxx
  LIBRARIES larcorealg_Geometry
  USE_BOOST_UNIT
  )


# set FHICL_FILE_PATH and FW_SEARCH_PATH for the case where building without mrb
# this is not generic, we take advantage of knowing the subdirectories for the fcl and gdml files
set( mrb_build_dir $ENV{MRB_BUILDDIR} )
if( NOT mrb_build_dir )
  #message(STATUS "MRB_BUILDDIR is undefined, set FHICL_FILE_PATH and FW_SEARCH_PATH")
  set_tests_properties( geometry_iterator_test geometry_test geometry_loader_test
                        geometry_iterator_loop_test geometry_standardchannelmapping_test
                        geometry_geoid_test geometry_thirdplaneslope_test
                        PROPERTIES ENVIRONMENT
                        "FHICL_FILE_PATH=.:${CMAKE_BINARY_DIR}/job;FW_SEARCH_PATH=${CMAKE_BINARY_DIR}/gdml" )
endif()


install_headers()
install_fhicl()
install_source()
