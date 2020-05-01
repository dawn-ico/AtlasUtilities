install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/stencils/interfaces 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/atlas_utils/stencils/
    FILES_MATCHING PATTERN "*.h*")

