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
#include <optional>
#include <string>

#include "mylib.hpp"

namespace {
void debugDump(const mylib::Grid& mesh, const std::string prefix) {
  {
    char buf[256];
    sprintf(buf, "%sT.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(const auto& cellIt : mesh.faces()) {
      fprintf(fp, "%d %d %d\n", cellIt.vertex(0).id() + 1, cellIt.vertex(1).id() + 1,
              cellIt.vertex(2).id() + 1);
    }
    fclose(fp);
  }

  {
    char buf[256];
    sprintf(buf, "%sP.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(const auto& nodeIt : mesh.vertices()) {
      double x = nodeIt.x();
      double y = nodeIt.y();
      fprintf(fp, "%f %f \n", x, y);
    }
    fclose(fp);
  }
}
} // namespace

mylib::Grid makeMylibMeshRect(int ny) {
  const auto& mesh = mylib::Grid(3 * ny, ny, false);

  double newHeight = (ny - 1) * sqrt(3) / 2.;
  double length = newHeight * 2;

  std::vector<int> keep;
  std::tuple<double, double> bblo{0., -std::numeric_limits<double>::max()};
  std::tuple<double, double> bbhi{length + 0.1 * (length / (2 * ny)),
                                  std::numeric_limits<double>::max()};
  auto inBB = [&](double x, double y) {
    return x > std::get<0>(bblo) && y > std::get<1>(bblo) && x < std::get<0>(bbhi) &&
           y < std::get<1>(bbhi);
  };
  for(const auto cellIt : mesh.faces()) {
    bool in = false;
    for(const auto nIt : cellIt.vertices()) {
      if(inBB(nIt->x(), nIt->y())) {
        in = true;
        break;
      }
    }
    if(in) {
      keep.push_back(cellIt.id());
    }
  }

  auto rectMesh = mylib::Grid(mesh, keep);
  double xMin = std::numeric_limits<double>::max();
  double yMin = std::numeric_limits<double>::max();
  double xMax = -std::numeric_limits<double>::max();
  double yMax = -std::numeric_limits<double>::max();
  for(auto nodeIt : rectMesh.vertices()) {
    double x = nodeIt.x();
    double y = nodeIt.y();
    xMin = fmin(x, xMin);
    yMin = fmin(y, yMin);
    xMax = fmax(x, xMax);
    yMax = fmax(y, yMax);
  }
  double lX = xMax - xMin;
  double lY = yMax - yMin;
  // re-center
  rectMesh.shift(-xMin - lX / 2, -yMin - lY / 2);
  // scale (single scale factor to exactly preserve equilateral edge lengths)
  double scale = 180 / lY;
  rectMesh.scale(scale);
  return rectMesh;
}