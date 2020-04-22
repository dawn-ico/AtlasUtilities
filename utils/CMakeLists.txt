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
  CMakeLists.txt
  GenerateRectAtlasMesh.cpp
  GenerateRectAtlasMesh.h
  GenerateRectToylibMesh.cpp
  GenerateRectToylibMesh.h
)
target_link_libraries(atlasUtilsLib toylib atlas eckit)
target_include_directories(atlasUtilsLib PUBLIC .)