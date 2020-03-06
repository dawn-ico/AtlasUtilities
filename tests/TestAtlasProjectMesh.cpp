//===--------------------------------------------------------------------------------*- C++ -*-===//
//                          _
//                         | |
//                       __| | __ ___      ___ ___
//                      / _` |/ _` \ \ /\ / / '_  |
//                     | (_| | (_| |\ V  V /| | | |
//                      \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include <assert.h>
#include <iostream>
#include <optional>
#include <string>

#include <atlas/library/Library.h>
#include <atlas/mesh/HybridElements.h>
#include <atlas/mesh/Mesh.h>
#include <atlas/mesh/Nodes.h>
#include <atlas/mesh/actions/BuildEdges.h>
#include <atlas/output/Gmsh.h>
#include <atlas/util/CoordinateEnums.h>

#include "AtlasFromNetcdf.h"
#include "AtlasProjectMesh.h"
#include "AtlasToNetcdf.h"

int main(int argc, char const* argv[]) {
  if(argc != 3) {
    std::cout << "intended use is\n"
              << argv[0] << " input_file.nc"
              << " output_file.nc" << std::endl;
    return -1;
  }
  std::string inFname(argv[1]);
  std::string outFname(argv[2]);

  auto meshInOpt = AtlasMeshFromNetCDFComplete(inFname);
  assert(meshInOpt.has_value());
  atlas::Mesh meshIn = meshInOpt.value();
  const int startFace = 5;
  const int numFace = 5;
  auto meshProjectedOpt = AtlasProjectMesh(meshIn, startFace, numFace);
  assert(meshProjectedOpt.has_value());
  atlas::Mesh meshProjected = meshProjectedOpt.value();
  AtlasToNetCDF(meshProjected, outFname);

  std::cout << "projection ran succesfully!\n";
}