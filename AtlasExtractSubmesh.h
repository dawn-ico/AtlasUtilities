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

#include <optional>
#include <string>

#include <atlas/mesh/Mesh.h>

// This module offers facilities to extract a submesh from a Atlas mesh
//
// NOTES: Minimal means that only the minimal set of of neighbors lists are adapted (i.e. only
// cellToNode). Complete adapts all neighbor lists

atlas::Mesh AtlasExtractSubMeshMinimal(const atlas::Mesh& mesh, std::pair<int, int> rangeCells);
atlas::Mesh AtlasExtractSubMeshComplete(const atlas::Mesh& mesh, std::pair<int, int> rangeCells);
atlas::Mesh AtlasExtractSubMeshMinimal(const atlas::Mesh& mesh,
                                       const std::vector<int> keptCellIndices);
atlas::Mesh AtlasExtractSubMeshComplete(const atlas::Mesh& mesh,
                                        const std::vector<int> keptCellIndices);