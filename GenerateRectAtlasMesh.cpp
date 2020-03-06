#include "GenerateRectAtlasMesh.h"

#include "AtlasExtractSubmesh.h"

#include <atlas/array.h>
#include <atlas/grid.h>
#include <atlas/meshgenerator.h>
#include <atlas/util/CoordinateEnums.h>

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

double TriangleArea(const atlas::Mesh& mesh, int cellIdx) {
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

  return fabs(((x0 * (y1 - y2)) + (x1 * (y2 - y0)) + (x2 * (y0 - y1))) * 0.5);
}

void debugDumpMeshRect(const atlas::Mesh& mesh, const std::string prefix) {
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

atlas::Mesh AtlasMeshRect(int ny) {
  atlas::Grid grid;
  int nx = 3 * ny;
  const bool dbgOut = false;

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

  if(dbgOut) {
    debugDumpMeshRect(mesh, "rightMesh");
  }

  auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    double x = xy(nodeIdx, atlas::LON);
    double y = xy(nodeIdx, atlas::LAT);
    x = x - 0.5 * y;
    y = y * sqrt(3) / 2.;
    xy(nodeIdx, atlas::LON) = x;
    xy(nodeIdx, atlas::LAT) = y;
  }

  if(dbgOut) {
    debugDumpMeshRect(mesh, "equiMesh");
  }

  double newHeight = (ny - 1) * sqrt(3) / 2.;
  double length = newHeight * 2;

  std::vector<int> keep;
  std::tuple<double, double> lo{0., -std::numeric_limits<double>::max()};
  std::tuple<double, double> hi{length + length / (nx)*0.1, std::numeric_limits<double>::max()};
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    if(TriangleInBB(mesh, cellIdx, lo, hi) && !TriangleArea(mesh, cellIdx) < 1e-12) {
      // in higher resolutions, atlas constructs some degenerate triangles (area = 0) at the bottom
      // of the mesh. I'm not sure if I'm using Atlas wrong here or if this is a bug in Atlas. For
      // now, use the hack above to filter such triangles
      keep.push_back(cellIdx);
    }
  }

  auto rectMesh = AtlasExtractSubMeshMinimal(mesh, keep);
  if(dbgOut) {
    debugDumpMeshRect(rectMesh, "rectMesh");
  }
  auto xyRect = atlas::array::make_view<double, 2>(rectMesh.nodes().xy());
  double xMin = std::numeric_limits<double>::max();
  double yMin = std::numeric_limits<double>::max();
  double xMax = -std::numeric_limits<double>::max();
  double yMax = -std::numeric_limits<double>::max();
  for(int nodeIdx = 0; nodeIdx < rectMesh.nodes().size(); nodeIdx++) {
    double x = xyRect(nodeIdx, atlas::LON);
    double y = xyRect(nodeIdx, atlas::LAT);
    xMin = fmin(x, xMin);
    yMin = fmin(y, yMin);
    xMax = fmax(x, xMax);
    yMax = fmax(y, yMax);
  }
  double lX = xMax - xMin;
  double lY = yMax - yMin;
  // re-center
  for(int nodeIdx = 0; nodeIdx < rectMesh.nodes().size(); nodeIdx++) {
    double x = xyRect(nodeIdx, atlas::LON);
    double y = xyRect(nodeIdx, atlas::LAT);
    xyRect(nodeIdx, atlas::LON) = x - xMin - lX / 2;
    xyRect(nodeIdx, atlas::LAT) = y - yMin - lY / 2;
  }

  // scale (single scale factor to exactly preserve equilateral edge lengths)
  double scale = 180 / lY;
  for(int nodeIdx = 0; nodeIdx < rectMesh.nodes().size(); nodeIdx++) {
    double x = xyRect(nodeIdx, atlas::LON);
    double y = xyRect(nodeIdx, atlas::LAT);
    xyRect(nodeIdx, atlas::LON) = x * scale;
    xyRect(nodeIdx, atlas::LAT) = y * scale;
  }

  if(dbgOut) {
    debugDumpMeshRect(mesh, "rectMeshScaleMove");
  }

  return rectMesh;
}