project(atlasIOLibrary)
add_library(atlasIOLib STATIC
  atlasIO.cpp
  atlasIO.h
  toylibIO.cpp
  toylibIO.h
)
target_link_libraries(atlasIOLib atlasUtilsLib atlas eckit)
target_include_directories(atlasIOLib PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/stencils/io>)

set_target_properties(atlasIOLib PROPERTIES PUBLIC_HEADER 
     "atlasIO.h;toylibIO.h"
)

install(TARGETS atlasIOLib EXPORT atlasutils 
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/atlas_utils/stencils/io  
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ 
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)
