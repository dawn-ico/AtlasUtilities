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

#pragma once

#include <atlas/mesh/Mesh.h>
#include <optional>
#include <string>

// This utility projects all elements on a range of the icosahdral faces onto the plane. This is
// HIGHLY EXPERIMENTAL and some very strong assumptions are made:
//  - all triangles on a ico-face are continous in memory
//  - all triangles on the ico-face in range [startFace, startFace+numFaces] are orientable by their
//    coordinate in cartesian 3 space. This seems not to be the case for the ICON base grids since
//    they underwent spring dynamics optimization. Use generated meshes with spring system optim
//    explicitly turned off for this
//
// The resulting triangles are mapped to [-180,180] x [-90, 90] to perform manufactured solution
// tests. Currently, there is no proper error handling and the method will most likely assert if the
// mesh does not conform to the assumptions above

std::optional<atlas::Mesh> AtlasProjectMesh(const atlas::Mesh& in, int startFace, int numFaces);