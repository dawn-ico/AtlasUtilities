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

#include "mylib.hpp"

// Generate a equilateral structured triangle mesh using our toy library with nx = 2*ny
// and [xmin, xmax] = [-180,180], [ymin, max] = [-90, 90]

mylib::Grid MylibMeshRect(int ny);