cmake_minimum_required(VERSION 3.5)

project(libdftd3 VERSION 0.9.2 LANGUAGES Fortran)

option(BUILD_SHARED_LIBS "Whether the library should be a shared one" FALSE)

option(INSTALL_INCLUDE_FILES "Whether include / module files should be installed" TRUE)

set(INSTALL_MODDIR "dftd3/modfiles" CACHE PATH
  "Installation directory for Fortran module files (within standard include folder)")

add_subdirectory(lib)
if(NOT LIBRARY_ONLY)
  add_subdirectory(prg)
endif()
