cet_make_library(LIBRARY_NAME AuxDetGeoObjectSorter INTERFACE
  SOURCE AuxDetGeoObjectSorter.h
  LIBRARIES INTERFACE
  larcoreobj::SimpleTypesAndConstants
  fhiclcpp::fhiclcpp
)

cet_make_library(LIBRARY_NAME Exceptions INTERFACE
  SOURCE Exceptions.h
  LIBRARIES INTERFACE
  larcoreobj::SimpleTypesAndConstants
  cetlib_except::cetlib_except
)

cet_make_library(LIBRARY_NAME GeometryIDmapper INTERFACE
  SOURCE GeometryIDmapper.h
  LIBRARIES INTERFACE
  larcoreobj::SimpleTypesAndConstants
)

cet_make_library(LIBRARY_NAME GeometryDataContainers INTERFACE
  SOURCE GeometryDataContainers.h
  LIBRARIES INTERFACE
  larcorealg::GeometryIDmapper
  larcoreobj::SimpleTypesAndConstants
  Boost::boost
)

cet_make_library(LIBRARY_NAME geo_vectors_utils INTERFACE
  SOURCE geo_vectors_utils.h
  LIBRARIES INTERFACE
  larcoreobj::geo_vectors
  ROOT::GenVector
)

cet_make_library(LIBRARY_NAME geo_vectors_utils_TVector INTERFACE
  SOURCE geo_vectors_utils_TVector.h
  LIBRARIES INTERFACE
  larcorealg::geo_vectors_utils
  ROOT::Physics
)

cet_make_library(LIBRARY_NAME LineClosestPoint INTERFACE
  SOURCE
  LineClosestPoint.h
  LineClosestPoint.tcc
  LIBRARIES INTERFACE
  larcorealg::geo_vectors_utils
  cetlib::cetlib
)

cet_make_library(LIBRARY_NAME Partitions INTERFACE
  SOURCE Partitions.h
  LIBRARIES INTERFACE
  larcorealg::CoreUtils
)

cet_make_library(LIBRARY_NAME ReadoutIDmapper INTERFACE
  SOURCE ReadoutIDmapper.h
  LIBRARIES INTERFACE
  larcorealg::GeometryIDmapper
  larcoreobj::SimpleTypesAndConstants
)

cet_make_library(LIBRARY_NAME ReadoutDataContainers INTERFACE
  SOURCE ReadoutDataContainers.h
  LIBRARIES INTERFACE
  larcorealg::GeometryDataContainers
  larcorealg::ReadoutIDmapper
  larcoreobj::SimpleTypesAndConstants
)

cet_make_library(LIBRARY_NAME StandaloneBasicSetup INTERFACE
  SOURCE StandaloneBasicSetup.h
  LIBRARIES INTERFACE
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
)

cet_make_library(LIBRARY_NAME TransformationMatrix INTERFACE
  SOURCE TransformationMatrix.h
  LIBRARIES INTERFACE
  larcorealg::Geometry
  ROOT::Matrix
)

cet_make_library(LIBRARY_NAME geo INTERFACE
  SOURCE geo.h
  LIBRARIES INTERFACE
  cetlib_except::cetlib_except
)

cet_make_library(LIBRARY_NAME geo_vectors_fhicl INTERFACE
  SOURCE geo_vectors_fhicl.h
  LIBRARIES INTERFACE
  larcoreobj::geo_vectors
  fhiclcpp::types
)

cet_make_library(LIBRARY_NAME geoVectorUtils INTERFACE
  NO_SOURCE
  LIBRARIES INTERFACE
  larcorealg::geo_vectors_fhicl
  larcorealg::geo_vectors_utils_TVector
  larcorealg::geo_vectors_utils
)

cet_make_library(SOURCE
  AuxDetGeo.cxx
  AuxDetGeometryCore.cxx
  AuxDetGeoObjectSorter.cxx
  AuxDetGeoObjectSorterStandard.cxx
  AuxDetReadoutGeom.cxx
  AuxDetSensitiveGeo.cxx
  BoxBoundedGeo.cxx
  CryostatGeo.cxx
  Decomposer.h
  GeometryBuilder.h
  GeometryBuilderStandard.cxx
  GeometryCore.cxx
  GeoNodePath.cxx
  GeoObjectSorter.cxx
  GeoObjectSorterStandard.cxx
  GeoVectorLocalTransformation.cxx
  Intersections.cxx
  LocalTransformation.cxx
  OpDetGeo.cxx
  PlaneGeo.cxx
  ROOTGeometryNavigator.h
  StandaloneGeometrySetup.cxx
  TPCGeo.cxx
  WireGeo.cxx
  WireReadoutGeom.cxx
  WireReadoutGeomBuilderStandard.cxx
  WireReadoutSorterStandard.cxx
  WireReadoutStandardGeom.cxx
  details/extractMaxGeometryElements.h
  details/GeometryIterationPolicy.cxx
  details/ReadoutIterationPolicy.cxx
  details/maybe_default_detector_name.cxx
  LIBRARIES
  PUBLIC
  larcorealg::CoreUtils
  larcorealg::GeometryDataContainers
  larcorealg::GeometryIDmapper
  larcorealg::LineClosestPoint
  larcorealg::Partitions
  larcorealg::ReadoutDataContainers
  larcorealg::geo_vectors_utils
  larcorealg::TransformationMatrix
  larcoreobj::SimpleTypesAndConstants
  larcoreobj::geo_vectors
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  CLHEP::Geometry
  CLHEP::Vector
  range-v3::range-v3
  ROOT::Core
  ROOT::GenVector
  ROOT::Geom
  ROOT::Physics
  PRIVATE
  Boost::filesystem
  larcorealg::Exceptions
  larcorealg::geo_vectors_utils_TVector
  messagefacility::MF_MessageLogger
  cetlib::container_algorithms
  ROOT::MathCore
)

install_headers(SUBDIRS "details")
install_fhicl(SUBDIRS "details")
install_source(SUBDIRS "details")

add_subdirectory(gdml)
