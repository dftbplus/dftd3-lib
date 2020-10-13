set(sources
  extras.f90
  main.f90)

add_executable(dftd3_standalone ${sources})

target_link_libraries(dftd3_standalone dftd3)

set_target_properties(dftd3_standalone PROPERTIES OUTPUT_NAME "dftd3")

install(TARGETS dftd3_standalone
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
