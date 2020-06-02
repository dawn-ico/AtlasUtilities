add_library(atlasUtilsLib STATIC IMPORTED)
find_library(atlas_utils_LIBRARY_PATH atlasUtilsLib HINTS "${CMAKE_CURRENT_LIST_DIR}/../../")
set_target_properties(atlasUtilsLib PROPERTIES IMPORTED_LOCATION "${atlas_utils_LIBRARY_PATH}")

set(atlas_utils_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/../../../include/)