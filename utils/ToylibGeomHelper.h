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

#include "../libs/toylib.hpp"

#include <tuple>

std::tuple<double, double> EdgeMidpoint(const toylib::Edge& e);
std::tuple<double, double> CellCircumcenter(const toylib::Face& c);
std::tuple<double, double> PrimalNormal(const toylib::Edge& e);
std::tuple<double, double> CellMidPoint(const toylib::Face& c);

double EdgeLength(const toylib::Edge& e);
double TriangleArea(const toylib::Vertex& v0, const toylib::Vertex& v1, const toylib::Vertex& v2);
double CellArea(const toylib::Face& c);
double DualCellArea(const toylib::Vertex& center);
double DualEdgeLength(const toylib::Edge& e);
double TangentOrientation(const toylib::Edge& e);

std::vector<toylib::Face> innerCells(const toylib::Grid& m);
std::vector<toylib::Edge> innerEdges(const toylib::Grid& m);
std::vector<toylib::Vertex> innerNodes(const toylib::Grid& m);