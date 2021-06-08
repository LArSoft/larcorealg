cet_make(
  SUBDIRS
    "details"
  LIBRARIES
          larcoreobj_SimpleTypesAndConstants
          messagefacility::MF_MessageLogger 
          fhiclcpp::fhiclcpp fhiclcpp::types
          cetlib
          cetlib_except
          ROOT::Core
          ROOT::Physics
          ROOT::Geom
          ROOT::GenVector
        )

install_headers(SUBDIRS "details")
install_fhicl(SUBDIRS "details")
install_source(SUBDIRS "details")

add_subdirectory(gdml)
