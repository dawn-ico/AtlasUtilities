//===--------------------------------------------------------------------------------*- C++ -*-===//
//                         _       _
//                        | |     | |
//                    __ _| |_ ___| | __ _ _ __   __ _
//                   / _` | __/ __| |/ _` | '_ \ / _` |
//                  | (_| | || (__| | (_| | | | | (_| |
//                   \__, |\__\___|_|\__,_|_| |_|\__, | - GridTools Clang DSL
//                    __/ |                       __/ |
//                   |___/                       |___/
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include "AtlasCartesianWrapper.h"
#include <atlas/util/CoordinateEnums.h>

namespace {
static double length(const Point& p1, const Point& p2) {
  double dx = std::get<0>(p1) - std::get<0>(p2);
  double dy = std::get<1>(p1) - std::get<1>(p2);
  return sqrt(dx * dx + dy * dy);
}

template <typename T>
static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

double TriangleArea(const Point& v0, const Point& v1, const Point& v2) {
  return fabs((std::get<0>(v0) * (std::get<1>(v1) - std::get<1>(v2)) +
               std::get<0>(v1) * (std::get<1>(v2) - std::get<1>(v0)) +
               std::get<0>(v2) * (std::get<1>(v0) - std::get<1>(v1))) *
              0.5);
}
} // namespace

Point AtlasToCartesian::cellMidpoint(const atlas::Mesh& mesh, int cellIdx) const {
  const atlas::mesh::HybridElements::Connectivity& cellNodeConnectivity =
      mesh.cells().node_connectivity();

  const int missingVal = cellNodeConnectivity.missing_value();
  int numNbh = cellNodeConnectivity.cols(cellIdx);
  double midX = 0;
  double midY = 0;
  for(int nbh = 0; nbh < numNbh; nbh++) {
    int nbhIdx = cellNodeConnectivity(cellIdx, nbh);
    assert(nbhIdx != missingVal);
    auto [nbhX, nbhY] = nodeLocation(nbhIdx);
    midX += nbhX;
    midY += nbhY;
  }
  midX /= numNbh;
  midY /= numNbh;
  return {midX, midY};
}

double AtlasToCartesian::cellArea(const atlas::Mesh& mesh, int cellIdx) const {
  const atlas::mesh::HybridElements::Connectivity& cellNodeConnectivity =
      mesh.cells().node_connectivity();

  const int missingVal = cellNodeConnectivity.missing_value();

  // only valid for tringular cells with all node neighbors set
  int numNbh = cellNodeConnectivity.cols(cellIdx);
  assert(numNbh == 3);
  for(int nbh = 0; nbh < numNbh; nbh++) {
    int nbhIdx = cellNodeConnectivity(cellIdx, nbh);
    assert(nbhIdx != missingVal);
  }

  auto v0 = nodeLocation(cellNodeConnectivity(cellIdx, 0));
  auto v1 = nodeLocation(cellNodeConnectivity(cellIdx, 1));
  auto v2 = nodeLocation(cellNodeConnectivity(cellIdx, 2));

  return TriangleArea(v0, v1, v2);
}

double AtlasToCartesian::dualCellArea(const atlas::Mesh& mesh, int nodeIdx) const {
  const atlas::mesh::Nodes::Connectivity& nodeEdgeConnectivity = mesh.nodes().edge_connectivity();
  const atlas::mesh::HybridElements::Connectivity& edgeCellConnectivity =
      mesh.edges().cell_connectivity();
  double totalArea = 0.;
  const int missingValEdge = nodeEdgeConnectivity.missing_value();
  const int missingValCell = edgeCellConnectivity.missing_value();

  auto center = nodeLocation(nodeIdx);

  int numNbh = nodeEdgeConnectivity.cols(nodeIdx);
  if(numNbh != 6) {
    return 0.;
  }

  for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
    int edgeIdx = nodeEdgeConnectivity(nodeIdx, nbhIdx);
    if(edgeIdx == missingValEdge) {
      return 0.;
    }

    int numNbhCells = edgeCellConnectivity.cols(edgeIdx);
    assert(numNbhCells == 2);

    int cellIdxLo = edgeCellConnectivity(edgeIdx, 0);
    int cellIdxHi = edgeCellConnectivity(edgeIdx, 1);

    if(cellIdxLo == missingValCell || cellIdxHi == missingValCell) {
      return 0.;
    }

    auto pLo = cellCircumcenter(mesh, cellIdxLo);
    auto pHi = cellCircumcenter(mesh, cellIdxHi);

    totalArea += TriangleArea(center, pLo, pHi);
  }
  return totalArea;
}

Point AtlasToCartesian::cellCircumcenter(const atlas::Mesh& mesh, int cellIdx) const {
  const atlas::mesh::HybridElements::Connectivity& cellNodeConnectivity =
      mesh.cells().node_connectivity();

  const int missingVal = cellNodeConnectivity.missing_value();

  // only valid for tringular cells with all node neighbors set
  int numNbh = cellNodeConnectivity.cols(cellIdx);
  assert(numNbh == 3);
  for(int nbh = 0; nbh < numNbh; nbh++) {
    int nbhIdx = cellNodeConnectivity(cellIdx, nbh);
    assert(nbhIdx != missingVal);
  }

  auto [Ax, Ay] = nodeLocation(cellNodeConnectivity(cellIdx, 0));
  auto [Bx, By] = nodeLocation(cellNodeConnectivity(cellIdx, 1));
  auto [Cx, Cy] = nodeLocation(cellNodeConnectivity(cellIdx, 2));

  double D = 2 * (Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By));
  double Ux = 1. / D *
              ((Ax * Ax + Ay * Ay) * (By - Cy) + (Bx * Bx + By * By) * (Cy - Ay) +
               (Cx * Cx + Cy * Cy) * (Ay - By));
  double Uy = 1. / D *
              ((Ax * Ax + Ay * Ay) * (Cx - Bx) + (Bx * Bx + By * By) * (Ax - Cx) +
               (Cx * Cx + Cy * Cy) * (Bx - Ax));
  return {Ux, Uy};
}

double AtlasToCartesian::distanceToCircumcenter(const atlas::Mesh& mesh, int cellIdx,
                                                int edgeIdx) const {
  const auto [x0, y0] = cellCircumcenter(mesh, cellIdx);
  const auto [e1, e2] = cartesianEdge(mesh, edgeIdx);
  // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
  const auto [x1, y1] = e1;
  const auto [x2, y2] = e2;
  return fabs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) /
         sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
}

Orientation AtlasToCartesian::edgeOrientation(const atlas::Mesh& mesh, int edgeIdx) const {

  Point lo;
  Point hi;
  {
    const atlas::mesh::HybridElements::Connectivity& edgeNodeConnectivity =
        mesh.edges().node_connectivity();
    const int missingVal = edgeNodeConnectivity.missing_value();

    int numNbh = edgeNodeConnectivity.cols(edgeIdx);
    assert(numNbh == 2);

    int nbhLo = edgeNodeConnectivity(edgeIdx, 0);
    int nbhHi = edgeNodeConnectivity(edgeIdx, 1);

    assert((nbhLo != missingVal && nbhHi != missingVal));

    lo = nodeToCartUnskewed[nbhLo];
    hi = nodeToCartUnskewed[nbhHi];
  }

  auto [loX, loY] = lo;
  auto [hiX, hiY] = hi;
  double dx = fabs(loX - hiX);
  double dy = fabs(loY - hiY);
  double tol = 1e3 * std::numeric_limits<double>::epsilon();
  assert(!(dx < tol && dy < tol));

  if(dx < tol && dy > tol) {
    return Orientation::Vertical;
  }

  if(dx > tol && dy < tol) {
    return Orientation::Horizontal;
  }

  return Orientation::Diagonal;
}

std::tuple<Point, Point> AtlasToCartesian::cartesianEdge(const atlas::Mesh& mesh,
                                                         int edgeIdx) const {
  const atlas::mesh::HybridElements::Connectivity& edgeNodeConnectivity =
      mesh.edges().node_connectivity();
  const int missingVal = edgeNodeConnectivity.missing_value();

  int numNbh = edgeNodeConnectivity.cols(edgeIdx);
  assert(numNbh == 2);

  int nbhLo = edgeNodeConnectivity(edgeIdx, 0);
  int nbhHi = edgeNodeConnectivity(edgeIdx, 1);

  assert((nbhLo != missingVal && nbhHi != missingVal));

  return {nodeToCart[nbhLo], nodeToCart[nbhHi]};
}

double AtlasToCartesian::edgeLength(const atlas::Mesh& mesh, int edgeIdx) const {
  auto [p1, p2] = cartesianEdge(mesh, edgeIdx);
  return length(p1, p2);
}
double AtlasToCartesian::dualEdgeLength(const atlas::Mesh& mesh, int edgeIdx) const {
  const atlas::mesh::HybridElements::Connectivity& edgeCellConnectivity =
      mesh.edges().cell_connectivity();
  const int missingVal = edgeCellConnectivity.missing_value();

  int numNbh = edgeCellConnectivity.cols(edgeIdx);
  assert(numNbh == 2);

  int nbhLo = edgeCellConnectivity(edgeIdx, 0);
  int nbhHi = edgeCellConnectivity(edgeIdx, 1);

  if(nbhLo == missingVal || nbhHi == missingVal) {
    return 0.;
  }

  Point pLo = cellCircumcenter(mesh, nbhLo);
  Point pHi = cellCircumcenter(mesh, nbhHi);

  return length(pLo, pHi);
}

double AtlasToCartesian::tangentOrientation(const atlas::Mesh& mesh, int edgeIdx) const {
  // ! =1 if vector product of vector from vertex1 to vertex 2 (v2-v1) by vector
  // ! from cell c1 to cell c2 (c2-c1) goes outside the sphere
  // ! =-1 if vector product ...       goes inside  the sphere

  const atlas::mesh::HybridElements::Connectivity& edgeCellConnectivity =
      mesh.edges().cell_connectivity();
  const int missingValCell = edgeCellConnectivity.missing_value();

  const atlas::mesh::HybridElements::Connectivity& edgeNodeConnectivity =
      mesh.edges().node_connectivity();
  const int missingValNode = edgeNodeConnectivity.missing_value();

  int numNbhNode = edgeNodeConnectivity.cols(edgeIdx);
  assert(numNbhNode == 2);

  int nbhLoNode = edgeNodeConnectivity(edgeIdx, 0);
  int nbhHiNode = edgeNodeConnectivity(edgeIdx, 1);

  assert((nbhLoNode != missingValNode && nbhHiNode != missingValNode));

  int numNbhCell = edgeCellConnectivity.cols(edgeIdx);
  assert(numNbhCell == 2);

  int nbhLoCell = edgeCellConnectivity(edgeIdx, 0);
  int nbhHiCell = edgeCellConnectivity(edgeIdx, 1);

  if(nbhLoCell == missingValCell || nbhHiCell == missingValCell) {
    return 1.; // not sure about this on the boundaries. chose 1 arbitrarily
  }

  Point pLoCell = cellCircumcenter(mesh, nbhLoCell);
  Point pHiCell = cellCircumcenter(mesh, nbhHiCell);

  Point pLoNode = nodeLocation(nbhLoNode);
  Point pHiNode = nodeLocation(nbhHiNode);

  double c2c1x = std::get<0>(pHiCell) - std::get<0>(pLoCell);
  double c2c1y = std::get<1>(pHiCell) - std::get<1>(pLoCell);

  double v2v1x = std::get<0>(pHiNode) - std::get<0>(pLoNode);
  double v2v1y = std::get<1>(pHiNode) - std::get<1>(pLoNode);

  return sgn(c2c1x * v2v1y - c2c1y * v2v1x);
}

Vector AtlasToCartesian::primalNormal(const atlas::Mesh& mesh, int edgeIdx) const {
  const atlas::mesh::HybridElements::Connectivity& edgeCellConnectivity =
      mesh.edges().cell_connectivity();
  const int missingValCell = edgeCellConnectivity.missing_value();

  int numNbhCell = edgeCellConnectivity.cols(edgeIdx);
  assert(numNbhCell == 2);

  int nbhLoCell = edgeCellConnectivity(edgeIdx, 0);
  int nbhHiCell = edgeCellConnectivity(edgeIdx, 1);

  if(nbhLoCell == missingValCell || nbhHiCell == missingValCell) {
    auto [p1, p2] = cartesianEdge(mesh, edgeIdx);
    double dx = std::get<0>(p1) - std::get<0>(p2);
    double dy = std::get<1>(p1) - std::get<1>(p2);
    double l = sqrt(dx * dx + dy * dy);
    double ex = dx / l;
    double ey = dy / l;
    return {ey, -ex};
  }

  Point pLoCell = cellCircumcenter(mesh, nbhLoCell);
  Point pHiCell = cellCircumcenter(mesh, nbhHiCell);

  double c2c1x = std::get<0>(pHiCell) - std::get<0>(pLoCell);
  double c2c1y = std::get<1>(pHiCell) - std::get<1>(pLoCell);

  double l = sqrt(c2c1x * c2c1x + c2c1y * c2c1y);
  return {c2c1x / l, c2c1y / l};
}

Point AtlasToCartesian::edgeMidpoint(const atlas::Mesh& mesh, int edgeIdx) const {
  auto [from, to] = cartesianEdge(mesh, edgeIdx);
  auto [fromX, fromY] = from;
  auto [toX, toY] = to;
  return {0.5 * (fromX + toX), 0.5 * (fromY + toY)};
}

std::vector<int> AtlasToCartesian::innerEdges(const atlas::Mesh& mesh) const {
  std::set<int> boundaryEdges;
  const auto& conn = mesh.nodes().edge_connectivity();
  auto anyMissing = [&](int nodeIdx) {
    int missingValue = conn.missing_value();
    for(int nbhIdx = 0; nbhIdx < conn.cols(nodeIdx); nbhIdx++) {
      if(conn(nodeIdx, nbhIdx) == missingValue) {
        return true;
      }
    }
    return false;
  };

  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    if(conn.cols(nodeIdx) != 6 || anyMissing(nodeIdx)) {
      for(int nbhIdx = 0; nbhIdx < conn.cols(nodeIdx); nbhIdx++) {
        boundaryEdges.insert(conn(nodeIdx, nbhIdx));
      }
    }
  }
  std::vector<int> innerEdges;
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    if(!boundaryEdges.count(edgeIdx)) {
      innerEdges.push_back(edgeIdx);
    }
  }

  return innerEdges;
}

std::vector<int> AtlasToCartesian::innerNodes(const atlas::Mesh& mesh) const {
  const auto& conn = mesh.nodes().edge_connectivity();
  auto anyMissing = [&](int nodeIdx) {
    int missingValue = conn.missing_value();
    for(int nbhIdx = 0; nbhIdx < conn.cols(nodeIdx); nbhIdx++) {
      if(conn(nodeIdx, nbhIdx) == missingValue) {
        return true;
      }
    }
    return false;
  };
  std::vector<int> innerNodes;
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    if(anyMissing(nodeIdx) || conn.cols(nodeIdx) != 6) {
      continue;
    }
    innerNodes.push_back(nodeIdx);
  }
  return innerNodes;
}

std::vector<int> AtlasToCartesian::innerCells(const atlas::Mesh& mesh) const {
  const auto& cellToEdge = mesh.cells().edge_connectivity();
  const auto& edgeToCell = mesh.edges().cell_connectivity();
  std::vector<int> innerCells;
  auto isBoundaryEdge = [&edgeToCell](int edgeIdx) {
    if(edgeToCell.cols(edgeIdx) != 2) {
      return true;
    }
    return edgeToCell(edgeIdx, 0) == edgeToCell.missing_value() ||
           edgeToCell(edgeIdx, 1) == edgeToCell.missing_value();
  };
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    int e0 = cellToEdge(cellIdx, 0);
    int e1 = cellToEdge(cellIdx, 1);
    int e2 = cellToEdge(cellIdx, 2);
    if(isBoundaryEdge(e0) || isBoundaryEdge(e1) || isBoundaryEdge(e2)) {
      continue;
    }
    innerCells.push_back(cellIdx);
  }
  return innerCells;
}

AtlasToCartesian::AtlasToCartesian(const atlas::Mesh& mesh, double scale, bool skewTrafo,
                                   bool center)
    : nodeToCart(mesh.nodes().size()), nodeToCartUnskewed(mesh.nodes().size()) {
  double lonLo = std::numeric_limits<double>::max();
  double lonHi = -std::numeric_limits<double>::max();
  double latLo = std::numeric_limits<double>::max();
  double latHi = -std::numeric_limits<double>::max();

  auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());

  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    lonLo = fmin(xy(nodeIdx, 0), lonLo);
    lonHi = fmax(xy(nodeIdx, 0), lonHi);
    latLo = fmin(xy(nodeIdx, 1), latLo);
    latHi = fmax(xy(nodeIdx, 1), latHi);
  }

  double uniscale = scale / (lonHi - lonLo);
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    double lon = xy(nodeIdx, 0);
    double lat = xy(nodeIdx, 1);
    double cartx = (lon - lonLo) * uniscale;
    double carty = (lat - latLo) * uniscale;

    if(center) {
      cartx -= 0.5 * scale;
      carty -= 0.5 * scale;
    }

    nodeToCartUnskewed[nodeIdx] = {cartx, carty};

    if(skewTrafo) {
      cartx = cartx - 0.5 * carty;
      carty = carty * sqrt(3) / 2.;
    }

    nodeToCart[nodeIdx] = {cartx, carty};
  }
}

AtlasToCartesian::AtlasToCartesian(const atlas::Mesh& mesh, bool scale)
    : nodeToCart(mesh.nodes().size()), nodeToCartUnskewed(mesh.nodes().size()) {
  auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
  for(int cellIdx = 0; cellIdx < mesh.nodes().size(); cellIdx++) {
    double cartx = scale ? xy(cellIdx, atlas::LON) / 180 * M_PI : xy(cellIdx, atlas::LON);
    double carty = scale ? xy(cellIdx, atlas::LAT) / 90 * 0.5 * M_PI : xy(cellIdx, atlas::LAT);
    nodeToCartUnskewed[cellIdx] = {cartx, carty};
    nodeToCart[cellIdx] = {cartx, carty};
  }
}