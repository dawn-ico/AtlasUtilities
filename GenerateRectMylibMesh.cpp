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

#include "GenerateRectMylibMesh.h"

mylib::Grid MylibMeshRect(int ny) {
  bool periodic = false;
  const auto& mesh = mylib::Grid(3 * ny, ny, periodic);

  double newHeight = (ny - 1) * sqrt(3) / 2.;
  double length = newHeight * 2;

  std::tuple<double, double> bblo{0., -std::numeric_limits<double>::max()};
  std::tuple<double, double> bbhi{length + 0.1 * (length / (2 * ny)),
                                  std::numeric_limits<double>::max()};

  auto inBB = [&](double x, double y) {
    return x > std::get<0>(bblo) && y > std::get<1>(bblo) && x < std::get<0>(bbhi) &&
           y < std::get<1>(bbhi);
  };
  auto anyInBB = [&](const mylib::Face& f) {
    for(const auto& v : f.vertices()) {
      if(inBB(v->x(), v->y())) {
        return true;
      }
    }
    return false;
  };

  std::vector<mylib::Face> keptCells;
  std::vector<int> keptIndices;
  std::copy_if(mesh.faces().begin(), mesh.faces().end(), std::back_inserter(keptCells), anyInBB);
  std::transform(keptCells.begin(), keptCells.end(), std::back_inserter(keptIndices),
                 [](const mylib::Face& c) { return c.id(); });
  auto rectMesh = mylib::Grid(mesh, keptIndices);
  double xMin =
      std::min_element(rectMesh.vertices().begin(), rectMesh.vertices().end(),
                       [](const mylib::Vertex& a, const mylib::Vertex& b) { return a.x() < b.x(); })
          ->x();
  double xMax =
      std::max_element(rectMesh.vertices().begin(), rectMesh.vertices().end(),
                       [](const mylib::Vertex& a, const mylib::Vertex& b) { return a.x() < b.x(); })
          ->x();
  double yMin =
      std::min_element(rectMesh.vertices().begin(), rectMesh.vertices().end(),
                       [](const mylib::Vertex& a, const mylib::Vertex& b) { return a.y() < b.y(); })
          ->y();
  double yMax =
      std::max_element(rectMesh.vertices().begin(), rectMesh.vertices().end(),
                       [](const mylib::Vertex& a, const mylib::Vertex& b) { return a.y() < b.y(); })
          ->y();

  double lX = xMax - xMin;
  double lY = yMax - yMin;
  // re-center
  rectMesh.shift(-xMin - lX / 2, -yMin - lY / 2);
  // scale (single scale factor to exactly preserve equilateral edge lengths)
  double scale = 180 / lY;
  rectMesh.scale(scale);
  return rectMesh;
}