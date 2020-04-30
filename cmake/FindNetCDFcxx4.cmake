# - Find NetCDF
# Find the NetCDF C++ includes and library
#
#  NETCDF_CXX4_INCLUDES    - where to find netcdf.h, etc
#  NETCDF_CXX4_LIBRARIES   - Link these libraries when using NetCDF
#  NETCDF_CXX4_FOUND       - True if NetCDF found including required interfaces (see below)
#
# Normal usage would be:
#  find_package (NetCDFcxx4 REQUIRED)
#  target_link_libraries (uses_f90_interface ${NETCDF_LIBRARIES})
#  target_link_libraries (only_uses_c_interface ${NETCDF_LIBRARIES_C})

find_path(NETCDF_CXX4_INCLUDES netcdf
     HINTS ${netcdfcxx4_DIR}/include)
find_library(NETCDF_CXX4_LIBRARIES netcdf_c++4
     HINTS ${netcdfcxx4_DIR}/lib)
message("LLLLLLL ${NETCDF_CXX4_INCLUDES} in ${netcdfcxx4_DIR} ${NETCDF_CXX4_LIBRARIES}")

set (NETCDF_LIBRARIES "${NetCDF_libs}" CACHE STRING "All NetCDF libraries required for interface level")

# handle the QUIETLY and REQUIRED arguments and set NETCDF_CXX4_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_CXX4_INCLUDES NETCDF_CXX4_LIBRARIES)

