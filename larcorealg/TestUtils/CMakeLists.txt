cet_make_library(LIBRARY_NAME ProviderTestHelpers INTERFACE
  SOURCE
  ProviderList.h
  ProviderTestHelpers.h
  LIBRARIES INTERFACE
  cetlib_except::cetlib_except
)

cet_make_library(LIBRARY_NAME StopWatch INTERFACE
  SOURCE StopWatch.h
)

cet_make_library(LIBRARY_NAME unit_test_base INTERFACE
  SOURCE unit_test_base.h
  LIBRARIES INTERFACE
  larcorealg::ProviderTestHelpers
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  cetlib::cetlib
)

cet_make_library(LIBRARY_NAME boost_unit_test_base INTERFACE
  SOURCE boost_unit_test_base.h
  LIBRARIES INTERFACE
  larcorealg::unit_test_base
  Boost::unit_test_framework
)

cet_make_library(LIBRARY_NAME geometry_unit_test_base INTERFACE
  SOURCE geometry_unit_test_base.h
  LIBRARIES INTERFACE
  larcorealg::unit_test_base
  larcorealg::Geometry
  messagefacility::MF_MessageLogger
  cetlib::cetlib
)

cet_make_library(SOURCE NameSelector.cxx
  LIBRARIES PRIVATE
  canvas::canvas
)

install_headers()
install_source()
