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
}
