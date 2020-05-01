project(atlasUtilsLibrary)
add_library(atlasUtilsLib STATIC
  AtlasCartesianWrapper.cpp
  AtlasCartesianWrapper.h
  AtlasExtractSubmesh.cpp
  AtlasExtractSubmesh.h
  AtlasFromNetcdf.cpp
  AtlasFromNetcdf.h
  AtlasProjectMesh.cpp
  AtlasProjectMesh.h
  AtlasToNetcdf.cpp
  AtlasToNetcdf.h  
  GenerateRectAtlasMesh.cpp
  GenerateRectAtlasMesh.h
  GenerateRectToylibMesh.cpp
  GenerateRectToylibMesh.h
  ToylibGeomHelper.cpp
  ToylibGeomHelper.h
)
target_link_libraries(atlasUtilsLib PUBLIC toylib atlas eckit)
target_link_libraries(atlasUtilsLib PUBLIC ${NETCDF_CXX4_LIBRARIES})

target_include_directories(atlasUtilsLib PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/utils> $<INSTALL_INTERFACE:include/atlas_utils/utils>)
target_include_directories(atlasUtilsLib PUBLIC ${NETCDF_CXX4_INCLUDES})

set_target_properties(atlasUtilsLib PROPERTIES PUBLIC_HEADER 
     "AtlasCartesianWrapper.h;AtlasExtractSubmesh.h;AtlasFromNetcdf.h;AtlasProjectMesh.h;AtlasToNetcdf.h;GenerateRectAtlasMesh.h;GenerateRectToylibMesh.h;ToylibGeomHelper.h"
)

install(TARGETS atlasUtilsLib EXPORT atlasutils 
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/atlas_utils/utils 
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ 
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)



