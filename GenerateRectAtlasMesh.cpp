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

#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/meshgenerator.h"
#include <atlas/library/Library.h>
#include <atlas/mesh/HybridElements.h>
#include <atlas/mesh/Mesh.h>
#include <atlas/mesh/Nodes.h>
#include <atlas/mesh/actions/BuildEdges.h>
#include <atlas/output/Gmsh.h>
#include <atlas/util/CoordinateEnums.h>

#include "AtlasExtractSubmesh.h"
#include "AtlasToNetcdf.h"

namespace {
bool TriangleInBB(const atlas::Mesh& mesh, int cellIdx, std::tuple<double, double> bblo,
                  std::tuple<double, double> bbhi) {
  const atlas::mesh::HybridElements::Connectivity& cellToNode = mesh.cells().node_connectivity();
  int node0 = cellToNode(cellIdx, 0);
  int node1 = cellToNode(cellIdx, 1);
  int node2 = cellToNode(cellIdx, 2);

  auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
  auto getXY = [xy](int nodeIdx) -> std::tuple<double, double> {
    return {xy(nodeIdx, atlas::LON), xy(nodeIdx, atlas::LAT)};
  };

  auto [x0, y0] = getXY(node0);
  auto [x1, y1] = getXY(node1);
  auto [x2, y2] = getXY(node2);

  auto inBB = [&](double x, double y) {
    return x > std::get<0>(bblo) && y > std::get<1>(bblo) && x < std::get<0>(bbhi) &&
           y < std::get<1>(bbhi);
  };

  return inBB(x0, y0) || inBB(x1, y1) || inBB(x2, y2);
}

atlas::Mesh makeAtlasMeshRect(int ny) {
  atlas::Grid grid;
  int nx = 3 * ny;

  // Create grid

  // this is adapted from
  // https://github.com/ecmwf/atlas/blob/a0017406f7ae54d306c9585113201af18d86fa40/src/tests/grid/test_grids.cc#L352
  //
  //    here, the grid is simple right triangles with strict up/down orientation. a transform will
  //    be applied later using the AtlasToCartesian wrapper to make the tris equilateral
  {
    using XSpace = atlas::StructuredGrid::XSpace;
    using YSpace = atlas::StructuredGrid::YSpace;

    // grid = atlas::StructuredGrid{XSpace{xspace}, YSpace{yspace}};
    auto x = atlas::grid::LinearSpacing(0, nx, nx, false);
    auto y = atlas::grid::LinearSpacing(0, ny, ny, false);
    grid = atlas::StructuredGrid{x, y};
  }

  auto meshgen = atlas::StructuredMeshGenerator{atlas::util::Config("angle", -1.)};
  auto mesh = meshgen.generate(grid);

  auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    double x = xy(nodeIdx, atlas::LON);
    double y = xy(nodeIdx, atlas::LAT);
    x = x - 0.5 * y;
    y = y * sqrt(3) / 2.;
    xy(nodeIdx, atlas::LON) = x;
    xy(nodeIdx, atlas::LAT) = y;
  }

  double newHeight = (ny - 1) * sqrt(3) / 2.;
  double length = newHeight * 2;

  std::vector<int> keep;
  std::tuple<double, double> lo{0., 0.};
  std::tuple<double, double> hi{length, std::numeric_limits<double>::max()};
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    if(TriangleInBB(mesh, cellIdx, lo, hi)) {
      keep.push_back(cellIdx);
    }
  }

  auto rectMesh = AtlasExtractSubMeshMinimal(mesh, keep);
  return rectMesh;
}

void debugDump(const atlas::Mesh& mesh, const std::string prefix) {
  auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
  const atlas::mesh::HybridElements::Connectivity& node_connectivity =
      mesh.cells().node_connectivity();

  {
    char buf[256];
    sprintf(buf, "%sT.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      int nodeIdx0 = node_connectivity(cellIdx, 0) + 1;
      int nodeIdx1 = node_connectivity(cellIdx, 1) + 1;
      int nodeIdx2 = node_connectivity(cellIdx, 2) + 1;
      fprintf(fp, "%d %d %d\n", nodeIdx0, nodeIdx1, nodeIdx2);
    }
    fclose(fp);
  }

  {
    char buf[256];
    sprintf(buf, "%sP.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      double x = xy(nodeIdx, atlas::LON);
      double y = xy(nodeIdx, atlas::LAT);
      fprintf(fp, "%f %f \n", x, y);
    }
    fclose(fp);
  }
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

  const auto& mesh = makeAtlasMeshRect(ny);

  debugDump(mesh, "test");
}