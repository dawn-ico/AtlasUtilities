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

#include <cmath>
#include <numeric>
#include <vector>

#include "AtlasCartesianWrapper.h"
#include "AtlasExtractSubmesh.h"
#include "AtlasFromNetcdf.h"
#include "AtlasToNetcdf.h"

#include <atlas/array.h>
#include <atlas/grid.h>
#include <atlas/meshgenerator.h>
#include <atlas/util/CoordinateEnums.h>

#include <netcdf>

namespace {
const double R = 6371; // radius earth

template <typename T>
static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

void debugDump(const atlas::Mesh& mesh, const std::string prefix) {
  auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
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
    auto latToRad = [](double rad) { return rad / 90. * (0.5 * M_PI); };
    auto lonToRad = [](double rad) { return rad / 180. * M_PI; };
    char buf[256];
    sprintf(buf, "%sP.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      double latRad = latToRad(lonlat(nodeIdx, atlas::LAT));
      double lonRad = lonToRad(lonlat(nodeIdx, atlas::LON));

      const double R = 6371; // radius earth
      double x = R * cos(latRad) * cos(lonRad);
      double y = R * cos(latRad) * sin(lonRad);
      double z = R * sin(latRad);
      fprintf(fp, "%f %f %f %f %f\n", x, y, z, latRad, lonRad);
    }
    fclose(fp);
  }

  // visualize with
  // P = load('netcdfMeshP.txt');
  // T = load('netcdfMeshT.txt');
  // trisurf(T(1:10,:),P(:,1),P(:,2),P(:,3))
}

void debugDumpUseXY(const atlas::Mesh& mesh, const std::string prefix) {
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

int GetHorizontalEdge(const atlas::Mesh& m, int cellIdx) {
  const auto& edgeToCell = m.edges().cell_connectivity();
  const auto& cellToEdge = m.cells().edge_connectivity();
  const auto& edgeToNode = m.edges().node_connectivity();

  auto latToRad = [](double rad) { return rad / 90. * (0.5 * M_PI); };
  auto lonlat = atlas::array::make_view<double, 2>(m.nodes().lonlat());

  auto verticalExtent = [&](int edgeIdx) {
    double lat0 = latToRad(lonlat(edgeToNode(edgeIdx, 0), atlas::LAT));
    double lat1 = latToRad(lonlat(edgeToNode(edgeIdx, 1), atlas::LAT));
    return fabs(sin(lat0) - sin(lat1));
  };

  double minVExtent = std::numeric_limits<double>::max();
  int horizontalEdge = -1;
  for(int nbhIdx = 0; nbhIdx < 3; nbhIdx++) {
    int edgeIdx = cellToEdge(cellIdx, nbhIdx);
    double vExtent = verticalExtent(edgeIdx);
    if(vExtent < minVExtent) {
      minVExtent = vExtent;
      horizontalEdge = nbhIdx;
    }
  }

  return horizontalEdge;
}

std::vector<int> NbhV(const atlas::Mesh& m, int cellIdx) {
  auto latToRad = [](double rad) { return rad / 90. * (0.5 * M_PI); };
  auto lonlat = atlas::array::make_view<double, 2>(m.nodes().lonlat());

  const auto& cellToEdge = m.cells().edge_connectivity();
  const auto& edgeToCell = m.edges().cell_connectivity();
  const auto& edgeToNode = m.edges().node_connectivity();
  assert(cellToEdge.cols(cellIdx) == 3);

  auto verticalExtent = [&](int edgeIdx) {
    double lat0 = latToRad(lonlat(edgeToNode(edgeIdx, 0), atlas::LAT));
    double lat1 = latToRad(lonlat(edgeToNode(edgeIdx, 1), atlas::LAT));
    return fabs(sin(lat0) - sin(lat1));
  };

  double minVExtent = std::numeric_limits<double>::max();
  int horizontalEdge = GetHorizontalEdge(m, cellIdx);

  std::vector<int> vNbh;
  for(int nbhIdx = 0; nbhIdx < 3; nbhIdx++) {
    int edgeIdx = cellToEdge(cellIdx, nbhIdx);
    bool boundary = edgeToCell(edgeIdx, 0) == edgeToCell.missing_value() ||
                    edgeToCell(edgeIdx, 1) == edgeToCell.missing_value();
    if(boundary) {
      continue;
    }
    if(nbhIdx == horizontalEdge) {
      continue;
    }
    vNbh.push_back(edgeToCell(edgeIdx, 0) == cellIdx ? edgeToCell(edgeIdx, 1)
                                                     : edgeToCell(edgeIdx, 0));
  }
  assert(vNbh.size() == 1 || vNbh.size() == 2);
  return vNbh;
}

int NbhH(const atlas::Mesh& m, int cellIdx) {
  auto latToRad = [](double rad) { return rad / 90. * (0.5 * M_PI); };
  auto lonlat = atlas::array::make_view<double, 2>(m.nodes().lonlat());

  const auto& cellToEdge = m.cells().edge_connectivity();
  const auto& edgeToCell = m.edges().cell_connectivity();
  const auto& edgeToNode = m.edges().node_connectivity();
  assert(cellToEdge.cols(cellIdx) == 3);

  auto verticalExtent = [&](int edgeIdx) {
    double lat0 = latToRad(lonlat(edgeToNode(edgeIdx, 0), atlas::LAT));
    double lat1 = latToRad(lonlat(edgeToNode(edgeIdx, 1), atlas::LAT));
    return fabs(sin(lat0) - sin(lat1));
  };

  double minVExtent = std::numeric_limits<double>::max();
  int horizontalEdge = GetHorizontalEdge(m, cellIdx);

  for(int nbhIdx = 0; nbhIdx < 3; nbhIdx++) {
    int edgeIdx = cellToEdge(cellIdx, nbhIdx);
    bool boundary = edgeToCell(edgeIdx, 0) == edgeToCell.missing_value() ||
                    edgeToCell(edgeIdx, 1) == edgeToCell.missing_value();
    if(boundary) {
      continue;
    }
    if(nbhIdx == horizontalEdge) {
      return (edgeToCell(edgeIdx, 0) == cellIdx ? edgeToCell(edgeIdx, 1) : edgeToCell(edgeIdx, 0));
    }
  }

  return -1;
}

bool HasBoundaryEdge(const atlas::Mesh& m, int cellIdx) {
  const auto& cellToEdge = m.cells().edge_connectivity();
  const auto& edgeToCell = m.edges().cell_connectivity();
  for(int nbhIdx = 0; nbhIdx < cellToEdge.cols(cellIdx); nbhIdx++) {
    int edgeIdx = cellToEdge(cellIdx, nbhIdx);
    if(edgeToCell(edgeIdx, 0) == edgeToCell.missing_value() ||
       edgeToCell(edgeIdx, 1) == edgeToCell.missing_value()) {
      return true;
    }
  }
  return false;
}

// note: returns true if any one node is in the BB
bool TriangleInBB(const atlas::Mesh& m, int cellIdx,
                  const std::vector<std::tuple<double, double>>& xy,
                  std::tuple<double, double> bblo, std::tuple<double, double> bbhi) {
  const atlas::mesh::HybridElements::Connectivity& cellToNode = m.cells().node_connectivity();
  int node0 = cellToNode(cellIdx, 0);
  int node1 = cellToNode(cellIdx, 1);
  int node2 = cellToNode(cellIdx, 2);

  auto [x0, y0] = xy[node0];
  auto [x1, y1] = xy[node1];
  auto [x2, y2] = xy[node2];

  auto inBB = [&](double x, double y) {
    return x > std::get<0>(bblo) && y > std::get<1>(bblo) && x < std::get<0>(bbhi) &&
           y < std::get<1>(bbhi);
  };

  return inBB(x0, y0) || inBB(x1, y1) || inBB(x2, y2);
}
} // namespace

std::optional<atlas::Mesh> AtlasProjectMesh(const atlas::Mesh& parentMesh, int startFace,
                                            int numFaces) {

  const int icosahedralFaces = 20;
  const int cellPerIcoFace = parentMesh.cells().size() / icosahedralFaces;
  auto subMesh = AtlasExtractSubMeshComplete(
      parentMesh,
      std::pair<int, int>{startFace * cellPerIcoFace, (startFace + numFaces) * cellPerIcoFace});

  const bool dbgOut = false;

  if(dbgOut) {
    debugDump(subMesh, "netcdfMesh");
  }

  auto latToRad = [](double rad) { return rad / 90. * (0.5 * M_PI); };
  auto lonToRad = [](double rad) { return rad / 180. * M_PI; };

  auto lonlat = atlas::array::make_view<double, 2>(subMesh.nodes().lonlat());
  const atlas::mesh::HybridElements::Connectivity& cellToNode = subMesh.cells().node_connectivity();
  const atlas::mesh::HybridElements::Connectivity& cellToEdge = subMesh.cells().edge_connectivity();
  const atlas::mesh::HybridElements::Connectivity& edgeToNode = subMesh.edges().node_connectivity();

  // for ico faces 5-? the cartesian orientation is equal to the topological orientation
  std::vector<int> orientation(subMesh.cells().size());
  for(int cellIdx = 0; cellIdx < subMesh.cells().size(); cellIdx++) {
    int nodeIdx0 = cellToNode(cellIdx, 0);
    int nodeIdx1 = cellToNode(cellIdx, 1);
    int nodeIdx2 = cellToNode(cellIdx, 2);

    double latRad0 = latToRad(lonlat(nodeIdx0, atlas::LAT));
    double latRad1 = latToRad(lonlat(nodeIdx1, atlas::LAT));
    double latRad2 = latToRad(lonlat(nodeIdx2, atlas::LAT));

    std::vector<double> z = {R * sin(latRad0), R * sin(latRad1), R * sin(latRad2)};
    std::sort(z.begin(), z.end());

    if(fabs(z[0] - z[1]) < fabs(z[1] - z[2])) {
      orientation[cellIdx] = -1;
    } else {
      orientation[cellIdx] = +1;
    }
  }

  if(dbgOut) {
    FILE* fp = fopen("orientation.txt", "w+");
    for(int cellIdx = 0; cellIdx < subMesh.cells().size(); cellIdx++) {
      fprintf(fp, "%d\n", orientation[cellIdx]);
    }
    fclose(fp);
  }

  auto toCart = [](double latRad, double lonRad, double R) {
    return std::tuple<double, double, double>{R * cos(latRad) * cos(lonRad),
                                              R * cos(latRad) * sin(lonRad), R * sin(latRad)};
  };

  // construct vertical boundary triangle strip
  std::vector<int> startCandidates;
  startCandidates.push_back(0);
  startCandidates.push_back(NbhV(subMesh, 0)[0]);
  bool lookHor = true;
  while(true) {
    if(lookHor) {
      int hNbh = NbhH(subMesh, startCandidates.back());
      if(hNbh == -1) {
        break;
      }
      startCandidates.push_back(hNbh);
      lookHor = false;
    } else {
      auto vNbh = NbhV(subMesh, startCandidates.back());
      startCandidates.push_back(startCandidates.end()[-2] == vNbh[0] ? vNbh[1] : vNbh[0]);
      lookHor = true;
    }
  }

  // drop triangles with no boundary edge
  std::vector<int> startCells;
  std::copy_if(startCandidates.begin(), startCandidates.end(), std::back_inserter(startCells),
               [&](const int candIdx) { return HasBoundaryEdge(subMesh, candIdx); });
  startCells.pop_back(); // pop last corner

  if(dbgOut) {
    FILE* fp = fopen("startStripe.txt", "w+");
    for(auto cellIdx : startCells) {
      for(int node = 0; node < 3; node++) {
        int nodeIdx = cellToNode(cellIdx, node);
        double lonRad = latToRad(lonlat(nodeIdx, atlas::LON));
        double latRad = latToRad(lonlat(nodeIdx, atlas::LAT));
        auto [x, y, z] = toCart(latRad, lonRad, R);
        fprintf(fp, "%f %f %f\n", x, y, z);
      }
    }
    fclose(fp);
  }

  // using the orientation computed each triangle can now be assigned a unique (I,J) index
  std::unordered_map<int, std::tuple<int, int>> cellIdxToIJ;
  for(int vIdx = 0; vIdx < startCells.size(); vIdx++) {
    std::vector<int> stripeI;
    int cellIdx0 = startCells[vIdx];
    int cellIdx1 = NbhV(subMesh, startCells[vIdx])[0];

    stripeI.push_back(cellIdx0);
    stripeI.push_back(cellIdx1);

    cellIdxToIJ.emplace(cellIdx0, std::tuple<int, int>{vIdx, 0});
    cellIdxToIJ.emplace(cellIdx1, std::tuple<int, int>{vIdx, 1});

    int hIdx = 2;

    while(true) {
      auto vNbh = NbhV(subMesh, stripeI.back());

      if(vNbh.size() == 1 && stripeI.end()[-2] == vNbh[0]) {
        break;
      }

      if(vNbh.size() == 1) {
        int cellIdx = vNbh[0];
        stripeI.push_back(cellIdx);
        cellIdxToIJ.emplace(cellIdx, std::tuple<int, int>{vIdx, hIdx});
        hIdx++;
        continue;
      }
      if(vNbh.size() == 2) {
        int cellIdx = stripeI.end()[-2] == vNbh[0] ? vNbh[1] : vNbh[0];
        stripeI.push_back(cellIdx);
        cellIdxToIJ.emplace(cellIdx, std::tuple<int, int>{vIdx, hIdx});
        hIdx++;
      }
    }

    if(dbgOut) {
      char fname[256];
      sprintf(fname, "stripe_%04d.txt", vIdx);
      FILE* fp = fopen(fname, "w+");
      for(auto cellIdx : stripeI) {
        for(int node = 0; node < 3; node++) {
          int nodeIdx = cellToNode(cellIdx, node);
          double lonRad = latToRad(lonlat(nodeIdx, atlas::LON));
          double latRad = latToRad(lonlat(nodeIdx, atlas::LAT));
          auto [x, y, z] = toCart(latRad, lonRad, R);
          fprintf(fp, "%f %f %f\n", x, y, z);
        }
      }
      fclose(fp);
    }
  }

  int nI = std::get<0>(max_element(cellIdxToIJ.begin(), cellIdxToIJ.end(),
                                   [](const std::pair<int, std::tuple<int, int>>& p1,
                                      const std::pair<int, std::tuple<int, int>>& p2) {
                                     return std::get<0>(p1.second) < std::get<0>(p1.second);
                                   })
                           ->second) +
           1;
  int nJ = std::get<1>(max_element(cellIdxToIJ.begin(), cellIdxToIJ.end(),
                                   [](const std::pair<int, std::tuple<int, int>>& p1,
                                      const std::pair<int, std::tuple<int, int>>& p2) {
                                     return std::get<1>(p1.second) < std::get<1>(p1.second);
                                   })
                           ->second) +
           1;

  if(dbgOut) {
    FILE* fpI = fopen("indexI.txt", "w+");
    FILE* fpJ = fopen("indexJ.txt", "w+");
    for(int cellIdx = 0; cellIdx < subMesh.cells().size(); cellIdx++) {
      int i = std::get<0>(cellIdxToIJ[cellIdx]);
      int j = std::get<1>(cellIdxToIJ[cellIdx]);
      fprintf(fpI, "%d\n", i);
      fprintf(fpJ, "%d\n", j);
    }
    fclose(fpI);
    fclose(fpJ);
  }

  std::vector<std::tuple<double, double>> newXY(subMesh.nodes().size(), {20, 20});
  std::vector<bool> written(subMesh.nodes().size(), false);
  const double l = 1.;
  const double h = 0.5 * sqrt(3);

  auto safeWrite = [&](int idx, std::tuple<double, double> in) {
    if(!written[idx]) {
      written[idx] = true;
      newXY[idx] = in;
    } else {
      auto [xOld, yOld] = newXY[idx];
      assert(fabs(std::get<0>(in) - xOld) < std::numeric_limits<double>::epsilon() * 1e3);
      assert(fabs(std::get<1>(in) - yOld) < std::numeric_limits<double>::epsilon() * 1e3);
      newXY[idx] = in;
    }
    newXY[idx] = in;
  };

  for(int cellIdx = 0; cellIdx < subMesh.cells().size(); cellIdx++) {
    // TODO: should I flip vertical index in order to not flip mesh topologically on its head
    // (lowest corner in space should have lowest index)?
    int vIdx = std::get<0>(cellIdxToIJ[cellIdx]);
    int hIdx = std::get<1>(cellIdxToIJ[cellIdx]) / 2;

    // each triangle takes care of its horizontal edge
    int hEdgeIdxLocal = GetHorizontalEdge(subMesh, cellIdx);
    int hEdgeIdx = cellToEdge(cellIdx, hEdgeIdxLocal);
    assert(hEdgeIdx != -1);
    int nodeIdx0 = edgeToNode(hEdgeIdx, 0);
    int nodeIdx1 = edgeToNode(hEdgeIdx, 1);

    std::vector<int> dbgnodes = {cellToNode(cellIdx, 0), cellToNode(cellIdx, 1),
                                 cellToNode(cellIdx, 2)};
    double lon0 = lonlat(nodeIdx0, atlas::LON);
    double lon1 = lonlat(nodeIdx1, atlas::LON);

    int leftNode, rightNode;
    if(lon0 < lon1) {
      leftNode = nodeIdx0;
      rightNode = nodeIdx1;
    } else {
      leftNode = nodeIdx1;
      rightNode = nodeIdx0;
    }

    if(orientation[cellIdx] == 1) { // up
      double x0 = (hIdx - 1) * l + 0.5 * l + 0.5 * l * vIdx;
      double y0 = (vIdx)*h;
      double x1 = (hIdx + 1) * l - 0.5 * l + 0.5 * l * vIdx;
      double y1 = (vIdx)*h;
      assert(std::find(dbgnodes.begin(), dbgnodes.end(), nodeIdx0) != dbgnodes.end());
      assert(std::find(dbgnodes.begin(), dbgnodes.end(), nodeIdx1) != dbgnodes.end());
      safeWrite(leftNode, {x0, y0});
      safeWrite(rightNode, {x1, y1});
    }

    if(orientation[cellIdx] == -1) { // down
      double x0 = (hIdx)*l + 0.5 * l * vIdx;
      double y0 = (vIdx + 1) * h;
      double x1 = (hIdx + 1) * l + 0.5 * l * vIdx;
      double y1 = (vIdx + 1) * h;
      assert(std::find(dbgnodes.begin(), dbgnodes.end(), nodeIdx0) != dbgnodes.end());
      assert(std::find(dbgnodes.begin(), dbgnodes.end(), nodeIdx1) != dbgnodes.end());
      safeWrite(leftNode, {x0, y0});
      safeWrite(rightNode, {x1, y1});
    }
  }

  //         _________________
  //        /                 /
  //       /                 /
  //      /_________________/
  //        ^_............_^
  //
  //        we are only interested in this part of the mesh

  // this would be the full area indicated above
  // std::tuple<double, double> bbLo(nI / 2 * l, -std::numeric_limits<double>::max());
  // std::tuple<double, double> bbHi(nJ / 2 * l - 0.5 * l, std::numeric_limits<double>::max());

  // however, for now, we want an aspect ratio of 1:2
  double height = h * nI;
  double lowX = nI / 2 * l;
  std::tuple<double, double> bbLo(lowX, -std::numeric_limits<double>::max());
  std::tuple<double, double> bbHi(lowX + 2 * height, std::numeric_limits<double>::max());

  std::vector<int> keep;
  for(int cellIdx = 0; cellIdx < subMesh.cells().size(); cellIdx++) {
    if(TriangleInBB(subMesh, cellIdx, newXY, bbLo, bbHi)) {
      keep.push_back(cellIdx);
    }
  }

  // create yet another atlas submesh only containing the rectangular subsection
  auto xyAtlas = atlas::array::make_view<double, 2>(subMesh.nodes().xy());
  for(int i = 0; i < newXY.size(); i++) {
    xyAtlas(i, atlas::LON) = std::get<0>(newXY[i]);
    xyAtlas(i, atlas::LAT) = std::get<1>(newXY[i]);
  }

  auto rectangularMesh = AtlasExtractSubMeshComplete(subMesh, keep);
  auto xyAtlasCompacted = atlas::array::make_view<double, 2>(rectangularMesh.nodes().xy());
  auto lonlatAtlasCompacted = atlas::array::make_view<double, 2>(rectangularMesh.nodes().lonlat());

  double xMin = std::numeric_limits<double>::max();
  double yMin = std::numeric_limits<double>::max();
  double xMax = -std::numeric_limits<double>::max();
  double yMax = -std::numeric_limits<double>::max();
  for(int nodeIdx = 0; nodeIdx < rectangularMesh.nodes().size(); nodeIdx++) {
    double x = xyAtlasCompacted(nodeIdx, atlas::LON);
    double y = xyAtlasCompacted(nodeIdx, atlas::LAT);
    xMin = fmin(x, xMin);
    yMin = fmin(y, yMin);
    xMax = fmax(x, xMax);
    yMax = fmax(y, yMax);
  }
  double lX = xMax - xMin;
  double lY = yMax - yMin;

  // re-center
  for(int nodeIdx = 0; nodeIdx < rectangularMesh.nodes().size(); nodeIdx++) {
    double x = xyAtlasCompacted(nodeIdx, atlas::LON);
    double y = xyAtlasCompacted(nodeIdx, atlas::LAT);
    xyAtlasCompacted(nodeIdx, atlas::LON) = x - xMin - lX / 2;
    xyAtlasCompacted(nodeIdx, atlas::LAT) = y - yMin - lY / 2;
  }

  // scale (single scale factor to exactly preserve equilateral edge lengths)
  double scale = 180 / lY;
  for(int nodeIdx = 0; nodeIdx < rectangularMesh.nodes().size(); nodeIdx++) {
    double x = xyAtlasCompacted(nodeIdx, atlas::LON);
    double y = xyAtlasCompacted(nodeIdx, atlas::LAT);
    xyAtlasCompacted(nodeIdx, atlas::LON) = x * scale;
    xyAtlasCompacted(nodeIdx, atlas::LAT) = y * scale;
    lonlatAtlasCompacted(nodeIdx, atlas::LON) = xyAtlasCompacted(nodeIdx, atlas::LON);
    lonlatAtlasCompacted(nodeIdx, atlas::LAT) = xyAtlasCompacted(nodeIdx, atlas::LAT);
  }

  if(dbgOut) {
    debugDumpUseXY(rectangularMesh, "rectMesh");
  }

  return rectangularMesh;
}