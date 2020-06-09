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

#include "../utils/GenerateRectAtlasMesh.h"

int main(int argc, char const* argv[]) {
  for(int nx = 3; nx < 30; nx++) {
    for(int ny = 1; ny < 30; ny++) {
      auto m = AtlasMeshRect(nx, ny);
      if(m.cells().size() != nx * ny) {
        printf("mesh size: %d expected: %d nx: %d ny: %d\n", m.cells().size(), nx * ny, nx, ny);
        assert(false);
      }
    }
  }

  auto m = AtlasMeshRect(4, 5);
  auto xy = atlas::array::make_view<double, 2>(m.nodes().xy());
  for(int vidx = 0; vidx != m.nodes().size(); ++vidx) {
    printf("v: %d   x: %f y: %f\n", vidx, xy(vidx, 0), xy(vidx, 1));
  }
  auto& cell2nodes = m.cells().node_connectivity();
  for(int cidx = 0; cidx != m.cells().size(); ++cidx) {
    printf("c: %d   v0: %d  v1: %d  v2: %d\n", cidx, cell2nodes(cidx, 0), cell2nodes(cidx, 1),
           cell2nodes(cidx, 2));
  }
  atlas::output::Gmsh gmsh("AtlasMeshRect4by4.msh");
  gmsh.write(m);
}
