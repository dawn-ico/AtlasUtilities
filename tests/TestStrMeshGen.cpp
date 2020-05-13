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
#include <mesh/Mesh.h>
#include <optional>
#include <string>

#include <atlas/array.h>
#include <atlas/grid.h>
#include <atlas/library/Library.h>
#include <atlas/mesh.h>
#include <atlas/mesh/actions/BuildEdges.h>
#include <atlas/meshgenerator.h>
#include <atlas/output/Gmsh.h>
#include <atlas/util/CoordinateEnums.h>

#include "../utils/GenerateStrIndxAtlasMesh.h"

int main(int argc, char const* argv[]) {
  int nx = 4;
  int ny = 4;
  auto m = AtlasStrIndxMesh(nx, ny);
  if(m.cells().size() != nx * ny) {
    printf("mesh size: %d expected: %d nx: %d ny: %d\n", m.cells().size(), nx * ny, nx, ny);
    assert(false);
  }
  auto& cell2edges = m.cells().edge_connectivity();
  auto& cell2vertices = m.cells().node_connectivity();
  for(int cidx = 0; cidx != m.cells().size(); ++cidx) {
    printf("cell %d : ", cidx);
    for(int color = 0; color < 3; color++) {
      auto edge = cell2edges(cidx, color);
      printf("edge%d = %d ", color, edge);
    }
    for(int color = 0; color < 3; color++) {
      auto vert = cell2vertices(cidx, color);
      printf("node%d = %d ", color, vert);
    }
    printf("\n");
  }
  auto& edge2verts = m.edges().node_connectivity();
  for(int eidx = 0; eidx != m.edges().size(); ++eidx) {
    printf("edge %d : ", eidx);
    for(int color = 0; color < 2; color++) {
      auto vertex = edge2verts(eidx, color);
      printf("vertex%d = %d ", color, vertex);
    }
    printf("\n");
  }

  printf("num vertices: %d\n", m.nodes().size());

  atlas::output::Gmsh gmsh("mesh.msh");
  gmsh.write(m);
}
