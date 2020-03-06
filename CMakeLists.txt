cmake_minimum_required(VERSION 3.12)
project(atlas_utils LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo." FORCE)
endif()

find_package(eckit REQUIRED)
find_package(atlas REQUIRED)
find_library(NETCDF_LIBRARY netcdf_c++4)

if (NETCDF_LIBRARY-NOTFOUND)
  message(FATAL_ERROR "netcdf not found")
endif()

if(NOT eckit_FOUND)
  message(FATAL_ERROR "Pass path containing eckit-config.cmake using -Deckit_DIR=<path/to/eckit>")
endif()

if(NOT atlas_FOUND)
  message(FATAL_ERROR "Pass path containing atlas-config.cmake using -Datlas_DIR=<path/to/atlas>")
endif()

add_subdirectory(libs)
add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(stencils)
add_subdirectory(scripts)

configure_file(resources/icon_160.nc ${PROJECT_BINARY_DIR}/tests COPYONLY)