#include "ToylibGeomHelper.h"

namespace {
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
} // namespace

double EdgeLength(const toylib::Edge& e) {
  double x0 = e.vertex(0).x();
  double y0 = e.vertex(0).y();
  double x1 = e.vertex(1).x();
  double y1 = e.vertex(1).y();
  double dx = x1 - x0;
  double dy = y1 - y0;
  return sqrt(dx * dx + dy * dy);
}

std::tuple<double, double> EdgeMidpoint(const toylib::Edge& e) {
  double x0 = e.vertex(0).x();
  double y0 = e.vertex(0).y();
  double x1 = e.vertex(1).x();
  double y1 = e.vertex(1).y();
  return {0.5 * (x0 + x1), 0.5 * (y0 + y1)};
}

std::tuple<double, double> CellCircumcenter(const toylib::Face& c) {
  double Ax = c.vertex(0).x();
  double Ay = c.vertex(0).y();
  double Bx = c.vertex(1).x();
  double By = c.vertex(1).y();
  double Cx = c.vertex(2).x();
  double Cy = c.vertex(2).y();

  double D = 2 * (Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By));
  double Ux = 1. / D *
              ((Ax * Ax + Ay * Ay) * (By - Cy) + (Bx * Bx + By * By) * (Cy - Ay) +
               (Cx * Cx + Cy * Cy) * (Ay - By));
  double Uy = 1. / D *
              ((Ax * Ax + Ay * Ay) * (Cx - Bx) + (Bx * Bx + By * By) * (Ax - Cx) +
               (Cx * Cx + Cy * Cy) * (Bx - Ax));
  return {Ux, Uy};
}

std::tuple<double, double> PrimalNormal(const toylib::Edge& e) {
  if(e.faces().size() != 2) {
    return {0., 0.};
  }

  auto [x0, y0] = CellCircumcenter(e.face(0));
  auto [x1, y1] = CellCircumcenter(e.face(1));
  double dx = x1 - x0;
  double dy = y1 - y0;
  double l = sqrt(dx * dx + dy * dy);
  return {dx / l, dy / l};
}

std::tuple<double, double> CellMidPoint(const toylib::Face& c) {
  auto v0 = c.vertex(0);
  auto v1 = c.vertex(1);
  auto v2 = c.vertex(2);
  return {1. / 3. * (v0.x() + v1.x() + v2.x()), 1. / 3. * (v0.y() + v1.y() + v2.y())};
}

double TriangleArea(const toylib::Vertex& v0, const toylib::Vertex& v1, const toylib::Vertex& v2) {
  return fabs(
      (v0.x() * (v1.y() - v2.y()) + v1.x() * (v2.y() - v0.y()) + v2.x() * (v0.y() - v1.y())) * 0.5);
}

double CellArea(const toylib::Face& c) {
  auto v0 = c.vertex(0);
  auto v1 = c.vertex(1);
  auto v2 = c.vertex(2);
  return TriangleArea(v0, v1, v2);
}

double DualCellArea(const toylib::Vertex& center) {
  double totalArea = 0.;
  for(const auto& e : center.edges()) {
    if(e->faces().size() != 2) {
      return 0.;
    }
    auto [leftx, lefty] = CellCircumcenter(e->face(0));
    auto [rightx, righty] = CellCircumcenter(e->face(1));
    toylib::Vertex left(leftx, lefty, -1);
    toylib::Vertex right(rightx, righty, -1);
    totalArea += TriangleArea(center, left, right);
  }
  return totalArea;
}

double DualEdgeLength(const toylib::Edge& e) {
  if(e.faces().size() == 1) { // dual edge length is zero on boundaries!
    return 0.;
  }
  auto c0 = e.face(0);
  auto c1 = e.face(1);
  auto [x0, y0] = CellCircumcenter(c0);
  auto [x1, y1] = CellCircumcenter(c1);
  double dx = x1 - x0;
  double dy = y1 - y0;
  return sqrt(dx * dx + dy * dy);
}

double TangentOrientation(const toylib::Edge& e) {
  if(e.faces().size() == 1) { // not sure about this on the boundaries. chose 1 arbitrarily
    return 1.;
  }

  auto c0 = e.face(0);
  auto c1 = e.face(1);
  auto [x0, y0] = CellCircumcenter(c0);
  auto [x1, y1] = CellCircumcenter(c1);
  double c2c1x = x1 - x0;
  double c2c1y = y1 - y0;

  auto v0 = e.vertex(0);
  auto v1 = e.vertex(1);
  double v2v1x = v1.x() - v0.x();
  double v2v1y = v1.y() - v0.y();

  return sgn(c2c1x * v2v1y - c2c1y * v2v1x);
}

std::vector<toylib::Face> innerCells(const toylib::Grid& m) {
  std::vector<toylib::Face> innerCells;
  for(const auto f : m.faces()) {
    bool hasBoundaryEdge = false;
    for(const auto e : f.edges()) {
      hasBoundaryEdge |= (e->faces().size() != 2);
    }
    if(hasBoundaryEdge) {
      continue;
    }
    innerCells.push_back(f);
  }
  return innerCells;
}
std::vector<toylib::Edge> innerEdges(const toylib::Grid& m) {
  std::vector<toylib::Edge> innerEdges;
  for(const auto e : m.edges()) {
    if(e.get().faces().size() != 2) {
      continue;
    }
    innerEdges.push_back(e);
  }
  return innerEdges;
}
std::vector<toylib::Vertex> innerNodes(const toylib::Grid& m) {
  std::vector<toylib::Vertex> innerVertices;
  for(const auto v : m.vertices()) {
    if(v.edges().size() != 6) {
      continue;
    }
    innerVertices.push_back(v);
  }
  return innerVertices;
}
