# ======================================================================
#
# Testing
#
# ======================================================================
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
cet_test( printBacktrace_test
  LIBRARIES
    larcorealg_CoreUtils
  )

cet_test(UncopiableAndUnmovableClass_test
  USE_BOOST_UNIT
  )


cet_test( ProviderPack_test USE_BOOST_UNIT )
cet_test( RealComparisons_test USE_BOOST_UNIT )
cet_test( NumericUtils_test USE_BOOST_UNIT )
cet_test( DereferenceIterator_test USE_BOOST_UNIT )
cet_test( DumpUtils_test
  LIBRARIES
    ROOT::Core
    ROOT::Physics
  USE_BOOST_UNIT
  )
cet_test( SortByPointers_test USE_BOOST_UNIT )
cet_test( span_test USE_BOOST_UNIT )
cet_test( MetaUtils_test USE_BOOST_UNIT )
cet_test( ContainerMeta_test USE_BOOST_UNIT )
cet_test( DebugUtils_test USE_BOOST_UNIT )
cet_test( counter_test USE_BOOST_UNIT )
cet_test( zip_test USE_BOOST_UNIT )
cet_test( enumerate_test USE_BOOST_UNIT )
cet_test( values_test USE_BOOST_UNIT )
cet_test( get_elements_test USE_BOOST_UNIT )
cet_test( operations_test USE_BOOST_UNIT )

# test quiet_Math_Functor.h
get_property(IncludeDirectories DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)

cet_test(MathFunctor_workaround HANDBUILT
  TEST_EXEC ./quiet_Math_Functor_test.sh
  TEST_ARGS ${IncludeDirectories}
  DATAFILES quiet_Math_Functor_test.sh
  )
