project(atlasIOLibrary)
add_library(atlasIOLib STATIC
  atlasIO.cpp
  atlasIO.h
  toylibIO.cpp
  toylibIO.h
)
target_link_libraries(atlasIOLib atlasUtilsLib atlas eckit)
target_include_directories(atlasIOLib PUBLIC .)