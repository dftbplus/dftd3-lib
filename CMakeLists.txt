cmake_minimum_required(VERSION 3.5)

project(libdftd3 VERSION 0.9.2 LANGUAGES Fortran)

set(LIBRARY_ONLY FALSE CACHE BOOL "Whether only library should be compiled")

set(INSTALL_BIN_DIR "bin")
set(INSTALL_MOD_DIR "include")
set(INSTALL_LIB_DIR "lib")

add_subdirectory(lib)
if(NOT LIBRARY_ONLY)
  add_subdirectory(prg)
endif()
