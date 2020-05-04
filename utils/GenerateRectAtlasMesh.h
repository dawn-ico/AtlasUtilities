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

// Generate equilateral structured triangle meshes in atlas

#pragma once

#include <atlas/mesh.h>

atlas::Mesh AtlasMeshRect(int ny);
atlas::Mesh AtlasMeshRect(int nx, int ny);
atlas::Mesh AtlasMeshSquare(int ny);