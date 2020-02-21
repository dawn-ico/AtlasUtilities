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

namespace {
void debugDump(const atlas::Mesh& mesh, const std::string prefix) {
  auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
  const atlas::mesh::HybridElements::Connectivity& node_connectivity =
      mesh.cells().node_connectivity();

  {
    char buf[256];
    sprintf(buf, "%sT.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      int nodeIdx0 = node_connectivity(cellIdx, 0) + 1;
      int nodeIdx1 = node_connectivity(cellIdx, 1) + 1;
      int nodeIdx2 = node_connectivity(cellIdx, 2) + 1;
      fprintf(fp, "%d %d %d\n", nodeIdx0, nodeIdx1, nodeIdx2);
    }
    fclose(fp);
  }

  {
    char buf[256];
    sprintf(buf, "%sP.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      double x = xy(nodeIdx, atlas::LON);
      double y = xy(nodeIdx, atlas::LAT);
      fprintf(fp, "%f %f \n", x, y);
    }
    fclose(fp);
  }
}
} // namespace

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
  debugDump(meshProjected, "PROJ_" + inFname);
  AtlasToNetCDF(meshProjected, outFname);
}