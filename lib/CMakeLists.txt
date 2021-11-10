set(sources
  api.f90
  common.f90
  core.f90
  pars.f90
  sizes.f90)

add_library(dftd3 ${sources})

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(dftd3 PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})

target_include_directories(dftd3 PUBLIC
  $<BUILD_INTERFACE:${includedir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODDIR}>)

install(TARGETS dftd3
  EXPORT ${INSTALL_EXPORT_NAME}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${includedir}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODDIR})
endif()

