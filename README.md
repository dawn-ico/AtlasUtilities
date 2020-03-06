# AtlasUtilities

This repoistory contains some utilities to transform Atlas meshes as well as the first ICON stencil prototype computing a Laplacian.

# Building

Requirements:

* libnetcdc-c++4
* Atlas
* eckit

Please make sure that libnetcdf is installed into the system path. Then:
```
git clone https://github.com/mroethlin/AtlasUtilities.git
cd AtlasUtilities
mkdir build
cd build
cmake .. -Deckit_DIR=<path/to/eckit> -Datlas_DIR=<path/to/atlas>
```
# Testing

Testing is _very_ basic. Simply change into the testing directory and run the bash script therein:

```
cd tests
bash run_tests.sh
```

All tests operate on the `icon_160.nc` input mesh provided. 

The idea is that asserts would trigger if something goes wrong. Furthermore, some output files can be inspected to gain further confidence that the utilities are working correctly:

* `earthV(1|2).msh` are minimal and complete Atlas meshes read from netcdf and can be visualized using gmsh. You should see a globe
* `outWriteNetcdf.nc` "round tripped" netcdf file. I.e. this was read and then written to disk again.
* `outProject.nc` a rectangular subsection of the input mesh which was cut out, then projected onto the plane and regularized into a equilateral, strcture triangle mesh.

In my (very brief) research I couldn't find a way to directly visualize those `*.nc` files. However, they can be passed to the `TestAtlasFromNetcdf` utility, which will convert them to `*.msh` files, which can in turn be visualized using `gmsh`. 

# Stencils

The stencils are located in `stencils`. Two versions are provided, one leveraging Atlas, the other using our toy library. Usage is simple:

```
./(mylib|atlas)IconLaplaceDriver <ny>
```

where `<ny>` is the horizontal resolution. A mesh of resultion `[nx,ny] = [2*ny, ny]` will be generated, and various error norms will be printed. Additionally, the divergence, curl and (normal) vector laplacian fields will be written to disk (`laplICON(mylib|atlas)_div.txt`, `laplICON(mylib|atlas)_rot.txt`, `laplICON(mylib|atlas)_out.txt`). The format is simply:

```
x y value
```
in ascii. 

There is also a python script that succesively increases the resolution to collect convergence data. Usage is:

```
python3 convergence_plot.py ./(mylib|atlas)IconLaplaceDriver
```

Various error norms are collected for divergence, curl and (normal) vector Laplacian. 

Various Octave scripts are provided to assist in visualizing debug and output data. For example, to get a convergence plot:

```
octave:1> convergence_plot_dec('<path/to/conv.csv') 
```
