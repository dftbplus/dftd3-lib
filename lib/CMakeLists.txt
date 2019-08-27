set(sources
  api.f90
  common.f90
  core.f90
  pars.f90
  sizes.f90)

add_library(dftd3 ${sources})

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(dftd3 PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})

target_include_directories(dftd3 PUBLIC ${includedir})

install(TARGETS dftd3
  EXPORT dftd3-targets
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR})

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${includedir}/ DESTINATION ${INSTALL_MOD_DIR})
endif()

