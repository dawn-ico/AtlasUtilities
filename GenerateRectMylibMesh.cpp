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

mylib::Grid makeMylibMeshRect(int ny) {
  const auto& mesh = mylib::Grid(3 * ny, ny, false);

  double newHeight = (ny - 1) * sqrt(3) / 2.;
  double length = newHeight * 2;

  std::vector<int> keep;
  std::tuple<double, double> bblo{0., 0.};
  std::tuple<double, double> bbhi{length, std::numeric_limits<double>::max()};
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
  return rectMesh;
}

} // namespace

int main(int argc, char const* argv[]) {
  if(argc != 3) {
    std::cout << "intended use is\n"
              << argv[0] << " yRes"
              << " output_file.nc" << std::endl;
    return -1;
  }

  int ny = atoi(argv[1]);
  std::string outFname(argv[2]);

  auto mesh = makeMylibMeshRect(ny);

  debugDump(mesh, "test");
}