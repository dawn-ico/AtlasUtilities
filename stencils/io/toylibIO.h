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

#include "../interfaces/toylib_interface.hpp"

#include "ToylibGeomHelper.h"

#include <string.h>
#include <optional>

void dumpMesh(const toylib::Grid& m, const std::string& fname);
void dumpDualMesh(const toylib::Grid& m, const std::string& fname);

void dumpSparseData(const toylib::Grid& mesh, const toylib::SparseVertexData<double>& sparseData,
                    int level, int edgesPerVertex, const std::string& fname);
void dumpSparseData(const toylib::Grid& mesh, const toylib::SparseFaceData<double>& sparseData,
                    int level, int edgesPerCell, const std::string& fname);

void dumpField(const std::string& fname, const toylib::Grid& mesh,
               const toylib::EdgeData<double>& field, int level,
               std::optional<toylib::edge_color> color = std::nullopt);
void dumpField(const std::string& fname, const toylib::Grid& mesh,
               const toylib::EdgeData<double>& field_x, const toylib::EdgeData<double>& field_y,
               int level, std::optional<toylib::edge_color> color = std::nullopt);
void dumpField(const std::string& fname, const toylib::Grid& mesh,
               const toylib::FaceData<double>& field, int level,
               std::optional<toylib::face_color> color = std::nullopt);
void dumpField(const std::string& fname, const toylib::Grid& mesh,
               const toylib::VertexData<double>& field, int level);
void debugDumpMesh(const toylib::Grid& mesh, const std::string prefix);
