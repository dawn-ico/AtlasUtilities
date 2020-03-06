# Atlas Utilites

Short description of each utility

* `AtlasCartesianWrapper` various helper functions to treat a Atlas mesh as if it was a planar mesh in cartesian coordinates. Can compute stuff like cell centroids, edge midpoint and the like. Some functions quietly assume that the mesh is triangular
* `AtlasExtractSubmesh` as the name suggests a submesh can be extracted from a Atlas mesh by providing a list of cell indices. Depending on which version is called, only the minimal or complete set of neighbor are copied over
* `AtlasFromNetcdf` reads a netcdf file and puts the results into the Atlas data structures. The resulting mesh is compatible with most of atlas, but not with parallelization, so no function spaces and no halos. The netcdf file is expected to follow the DWD naming conventions. Again, either all neighbor lists present in the netcdf are read or only the minimal set. For the latter option Atlas actions can be used to retrieve the complete set of neighbor lists again
* `AtlasToNetcdf` as above, but the other way around.
* `GenerateRectAtlasMesh` a Atlas mesh generator that generates a rectangular mesh of equilateral triangles in a "up, down" topology. Uses `AtlasExtractSubmesh`. Again, no parallelization and no halo regions.
* `GenerateRectMylibMesh` same as above, but for our toy library. Thus, strictly speaking not a Atlas utility. 