
cet_test(NameSelector_test
  SOURCES NameSelector_test.cxx
  LIBRARIES
    larcorealg_TestUtils
    canvas::canvas
  USE_BOOST_UNIT
)

cet_test(ProviderList_test
  USE_BOOST_UNIT
)

cet_test(StopWatch_test)
