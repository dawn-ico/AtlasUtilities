project(toylibLibrary)
add_library(toylib STATIC
                toylib.hpp
                toylib.cpp
)
target_include_directories(toylib PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/libs> $<INSTALL_INTERFACE:include/atlas_utils/libs> $<INSTALL_INTERFACE:include/>)

set_target_properties(toylib PROPERTIES PUBLIC_HEADER "toylib.hpp")

install(TARGETS toylib EXPORT atlasutils 
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/atlas_utils/libs
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ 
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)
