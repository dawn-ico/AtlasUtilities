add_subdirectory(io)
add_subdirectory(interfaces)

add_executable(atlasIconLaplaceDriver atlasIconLaplaceDriver.cpp)
target_link_libraries(atlasIconLaplaceDriver atlas eckit atlasUtilsLib atlasIOLib ${NetCDF_C_LIBRARY})

add_executable(atlasIconDiamondLaplacianDriver atlasIconDiamondLaplacianDriver.cpp)
target_link_libraries(atlasIconDiamondLaplacianDriver atlas eckit atlasUtilsLib atlasIOLib ${NetCDF_C_LIBRARY})

add_executable(atlasShallowWater shallowWater.cpp)
target_link_libraries(atlasShallowWater atlas eckit atlasUtilsLib ${NetCDF_C_LIBRARY})

add_executable(mylibIconLaplaceDriver mylibIconLaplaceDriver.cpp)
target_link_libraries(mylibIconLaplaceDriver atlasUtilsLib toylib atlasIOLib ${NetCDF_C_LIBRARY})

install(TARGETS atlasIconLaplaceDriver atlasIconDiamondLaplacianDriver atlasShallowWater mylibIconLaplaceDriver 
       EXPORT atlas_utils DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/)

