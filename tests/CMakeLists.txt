add_executable(TestAtlasFromNetcdf TestAtlasFromNetcdf.cpp)
target_link_libraries(TestAtlasFromNetcdf atlas eckit atlasUtilsLib ${NETCDF_CXX4_LIBRARY}  ${NetCDF_C_LIBRARY})

add_executable(TestAtlasProjectMesh TestAtlasProjectMesh.cpp)
target_link_libraries(TestAtlasProjectMesh atlas eckit atlasUtilsLib ${NETCDF_CXX4_LIBRARY}  ${NetCDF_C_LIBRARY})

add_executable(TestAtlasToNetcdf TestAtlasToNetcdf.cpp)
target_link_libraries(TestAtlasToNetcdf atlas eckit atlasUtilsLib ${NETCDF_CXX4_LIBRARY} ${NetCDF_C_LIBRARY})

add_executable(TestAtlasMeshGen TestAtlasMeshGen.cpp)
target_link_libraries(TestAtlasMeshGen atlas eckit atlasUtilsLib ${NETCDF_CXX4_LIBRARY} ${NetCDF_C_LIBRARY})

add_executable(TestStrMeshGen TestStrMeshGen.cpp)
target_link_libraries(TestStrMeshGen atlas eckit atlasUtilsLib ${NETCDF_CXX4_LIBRARY} ${NetCDF_C_LIBRARY})

