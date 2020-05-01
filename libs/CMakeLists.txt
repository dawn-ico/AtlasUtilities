project(toylibLibrary)
add_library(toylib STATIC
                toylib.hpp
                toylib.cpp
)
target_include_directories(toylib PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/libs>)
set_target_properties(toylib PROPERTIES PUBLIC_HEADER "toylib.hpp")

install(TARGETS toylib EXPORT atlasutils 
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ 
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)