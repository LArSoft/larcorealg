# ======================================================================
#
# Testing
#
# ======================================================================

cet_test(printBacktrace_test
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
cet_test(UncopiableAndUnmovableClass_test USE_BOOST_UNIT)
cet_test(ProviderPack_test USE_BOOST_UNIT)
cet_test(RealComparisons_test USE_BOOST_UNIT)
cet_test(NumericUtils_test USE_BOOST_UNIT)
cet_test(DumpUtils_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  ROOT::Physics
)
cet_test(SortByPointers_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
cet_test(span_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
  Boost::headers
)
cet_test(MetaUtils_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
cet_test(ContainerMeta_test USE_BOOST_UNIT)
cet_test(DebugUtils_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
cet_test(counter_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
cet_test(zip_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
cet_test(enumerate_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
cet_test(values_test USE_BOOST_UNIT)
cet_test(get_elements_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
cet_test(operations_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::CoreUtils
)
