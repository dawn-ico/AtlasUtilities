cmake_minimum_required(VERSION 3.12)
project(atlas_utils LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake" CACHE STRING "Modules for CMake" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo." FORCE)
endif()

find_package(eckit REQUIRED)
find_package(atlas REQUIRED)
find_package(NetCDF REQUIRED)
find_package (NetCDFcxx4 REQUIRED)

add_subdirectory(libs)
add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(stencils)
add_subdirectory(scripts)

configure_file(resources/icon_160.nc ${PROJECT_BINARY_DIR}/tests COPYONLY)
install(EXPORT atlasutils DESTINATION cmake NAMESPACE au:: FILE atlas_utils-config.cmake)
