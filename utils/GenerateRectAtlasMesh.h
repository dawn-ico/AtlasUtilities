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

// Generate a equilateral structured triangle mesh using our Atlas with nx = 2*ny
// and [xmin, xmax] = [-180,180], [ymin, max] = [-90, 90]

#pragma once

#include <atlas/mesh.h>

std::tuple<atlas::Grid, atlas::Mesh> AtlasMeshRect(int ny);
std::tuple<atlas::Grid, atlas::Mesh> AtlasMeshSquare(int ny);
void generateCell2CellTable(atlas::Mesh& mesh, bool allocate);