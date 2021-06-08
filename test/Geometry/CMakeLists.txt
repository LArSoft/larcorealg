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
                        messagefacility::MF_MessageLogger 
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
