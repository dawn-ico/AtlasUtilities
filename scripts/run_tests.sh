#!/bin/bash
echo "Testing reading from netcdf...."
./TestAtlasFromNetcdf icon_160.nc
echo ""

echo "Testing writing to netcdf...."
./TestAtlasToNetcdf icon_160.nc outWriteNetcdf.nc 
echo ""

echo "Testing mesh projection....."
./TestAtlasProjectMesh icon_160.nc outProject.nc 
echo ""