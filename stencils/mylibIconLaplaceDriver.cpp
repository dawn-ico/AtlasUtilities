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

// This file implements the nabla2vec function from mo_math_laplace.f90 (see also
// mo_math_divrot.f90) in ICON, using out toylibarry, dawn, and the toy library interface therein.
// Some notes:
//
//  - MOST CODE COMMENTS HAVE BEEN OMITTED, SEE atlasIconLaplaceDriver.cpp for EXPLANATIONS WHAT
//    HAPPENS IN THE CODE
//
//  - names have been kept close to the FORTRAN code, but the "_Location" suffixes have been removed
//    because of the strong typing in C++ and inconsistent application in the FORTRAN source
//
//  - the main function in this file either accepts a vertical resolution nx, or reads a netcdf mesh
//    from disk. in the latter case, the netcdf file needs to contain a structures equilateral
//    triangle mesh
//
//  - this version makes no attempt to compute anything meaningful at the boundaries. values at the
//    boundaries are skipped in outputs, meaningless default values are assigned to various
//    geometrical factors etc.

#include <assert.h>
#include <fstream>
#include <optional>

#include "interfaces/mylib_interface.hpp"
#include "mylib.hpp"

#include "generated_iconLaplace.hpp"

#include "GenerateRectMylibMesh.h"

namespace {
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

//===------------------------------------------------------------------------------------------===//
// geometric helper functions
//===------------------------------------------------------------------------------------------===//
std::tuple<double, double> EdgeMidpoint(const mylib::Edge& e);
std::tuple<double, double> CellCircumcenter(const mylib::Face& c);
std::tuple<double, double> PrimalNormal(const mylib::Edge& e);
std::tuple<double, double> CellMidPoint(const mylib::Face& c);

double EdgeLength(const mylib::Edge& e);
double TriangleArea(const mylib::Vertex& v0, const mylib::Vertex& v1, const mylib::Vertex& v2);
double CellArea(const mylib::Face& c);
double DualCellArea(const mylib::Vertex& center);
double DualEdgeLength(const mylib::Edge& e);
double TangentOrientation(const mylib::Edge& e);

std::vector<mylib::Face> innerCells(const mylib::Grid& m);
std::vector<mylib::Edge> innerEdges(const mylib::Grid& m);
std::vector<mylib::Vertex> innerNodes(const mylib::Grid& m);

//===------------------------------------------------------------------------------------------===//
// output (debugging)
//===------------------------------------------------------------------------------------------===//
void dumpMesh(const mylib::Grid& m, const std::string& fname);
void dumpDualMesh(const mylib::Grid& m, const std::string& fname);

void dumpSparseData(const mylib::Grid& mesh, const mylib::SparseVertexData<double>& sparseData,
                    int level, int edgesPerVertex, const std::string& fname);
void dumpSparseData(const mylib::Grid& mesh, const mylib::SparseFaceData<double>& sparseData,
                    int level, int edgesPerCell, const std::string& fname);

void dumpField(const std::string& fname, const mylib::Grid& mesh,
               const mylib::EdgeData<double>& field, int level,
               std::optional<mylib::edge_color> color = std::nullopt);
void dumpField(const std::string& fname, const mylib::Grid& mesh,
               const mylib::EdgeData<double>& field_x, const mylib::EdgeData<double>& field_y,
               int level, std::optional<mylib::edge_color> color = std::nullopt);
void dumpField(const std::string& fname, const mylib::Grid& mesh,
               const mylib::FaceData<double>& field, int level,
               std::optional<mylib::face_color> color = std::nullopt);
void dumpField(const std::string& fname, const mylib::Grid& mesh,
               const mylib::VertexData<double>& field, int level);
void debugDumpMesh(const mylib::Grid& mesh, const std::string prefix);

//===------------------------------------------------------------------------------------------===//
// error reporting
//===------------------------------------------------------------------------------------------===//
template <typename DataT, typename ElemT>
std::tuple<double, double, double> MeasureError(DataT ref, DataT sol,
                                                std::vector<ElemT> innerElements, int level);

} // namespace

int main(int argc, char const* argv[]) {
  if(argc != 2) {
    std::cout << "intended use is\n" << argv[0] << " ny" << std::endl;
    return -1;
  }
  int w = atoi(argv[1]);

  int k_size = 1;
  const int level = 0;
  double lDomain = M_PI;

  const bool dbg_out = false;

  mylib::Grid mesh = MylibMeshRect(w);

  mesh.scale(M_PI / 180.); // rescale to radians

  if(dbg_out) {
    debugDumpMesh(mesh, "laplICONmylib_Mesh");
  }

  const int edgesPerVertex = 6;
  const int edgesPerCell = 3;

  //===------------------------------------------------------------------------------------------===//
  // input field (field we want to take the laplacian of)
  //===------------------------------------------------------------------------------------------===//
  mylib::EdgeData<double> vec(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // control fields (containing analytical solutions)
  //===------------------------------------------------------------------------------------------===//
  mylib::FaceData<double> divVecSol(mesh, k_size);
  mylib::VertexData<double> rotVecSol(mesh, k_size);
  mylib::EdgeData<double> lapVecSol(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // output field (field containing the computed laplacian)
  //===------------------------------------------------------------------------------------------===//
  mylib::EdgeData<double> nabla2_vec(mesh, k_size);
  // term 1 and term 2 of nabla for debugging
  mylib::EdgeData<double> nabla2t1_vec(mesh, k_size);
  mylib::EdgeData<double> nabla2t2_vec(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // intermediary fields (curl/rot and div of vec_e)
  //===------------------------------------------------------------------------------------------===//
  mylib::VertexData<double> rot_vec(mesh, k_size);
  mylib::FaceData<double> div_vec(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // sparse dimensions for computing intermediary fields
  //===------------------------------------------------------------------------------------------===//
  mylib::SparseVertexData<double> geofac_rot(mesh, k_size, edgesPerVertex);
  mylib::SparseVertexData<double> edge_orientation_vertex(mesh, k_size, edgesPerVertex);

  mylib::SparseFaceData<double> geofac_div(mesh, k_size, edgesPerCell);
  mylib::SparseFaceData<double> edge_orientation_cell(mesh, k_size, edgesPerCell);

  //===------------------------------------------------------------------------------------------===//
  // fields containing geometric information
  //===------------------------------------------------------------------------------------------===//
  mylib::EdgeData<double> tangent_orientation(mesh, k_size);
  mylib::EdgeData<double> primal_edge_length(mesh, k_size);
  mylib::EdgeData<double> dual_edge_length(mesh, k_size);
  mylib::EdgeData<double> dual_normal_x(mesh, k_size);
  mylib::EdgeData<double> dual_normal_y(mesh, k_size);
  mylib::EdgeData<double> primal_normal_x(mesh, k_size);
  mylib::EdgeData<double> primal_normal_y(mesh, k_size);

  mylib::FaceData<double> cell_area(mesh, k_size);
  mylib::VertexData<double> dual_cell_area(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  for(auto const& e : mesh.edges()) {
    primal_edge_length(e, level) = EdgeLength(e);
    dual_edge_length(e, level) = DualEdgeLength(e);
    tangent_orientation(e, level) = TangentOrientation(e);
    auto [xm, ym] = EdgeMidpoint(e);
    auto [nx, ny] = PrimalNormal(e);
    primal_normal_x(e, level) = nx;
    primal_normal_y(e, level) = ny;
    // The primal normal, dual normal
    // forms a left-handed coordinate system
    dual_normal_x(e, level) = ny;
    dual_normal_y(e, level) = -nx;
  }
  if(dbg_out) {
    dumpField("laplICONmylib_tangentOrientation.txt", mesh, tangent_orientation, level);
    dumpField("laplICONmylib_EdgeLength.txt", mesh, primal_edge_length, level);
    dumpField("laplICONmylib_dualEdgeLength.txt", mesh, dual_edge_length, level);
    dumpField("laplICONmylib_nrm.txt", mesh, primal_normal_x, primal_normal_y, level);
    dumpField("laplICONmylib_dnrm.txt", mesh, dual_normal_x, dual_normal_y, level);
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on cells
  //===------------------------------------------------------------------------------------------===//
  for(const auto& c : mesh.faces()) {
    cell_area(c, level) = CellArea(c);
  }
  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on vertices
  //===------------------------------------------------------------------------------------------===//
  for(const auto& v : mesh.vertices()) {
    dual_cell_area(v, level) = DualCellArea(v);
  }

  if(dbg_out) {
    dumpField("laplICONmylib_areaCell.txt", mesh, cell_area, level);
    dumpField("laplICONmylib_areaCellDual.txt", mesh, dual_cell_area, level);
  }

  //===------------------------------------------------------------------------------------------===//
  // analytical solutions
  //===------------------------------------------------------------------------------------------===//
  auto sphericalHarmonic = [](double x, double y) -> std::tuple<double, double> {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return {c1 * cos(2 * x) * cos(y) * cos(y) * sin(y), c2 * cos(x) * cos(y) * sin(y)};
  };
  auto analyticalDivergence = [](double x, double y) {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return -2 * c1 * sin(2 * x) * cos(y) * cos(y) * sin(y) +
           c2 * cos(x) * (cos(y) * cos(y) - sin(y) * sin(y));
  };
  auto analyticalCurl = [](double x, double y) {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    double dudy = c1 * cos(2 * x) * cos(y) * (cos(y) * cos(y) - 2 * sin(y) * sin(y));
    double dvdx = -c2 * cos(y) * sin(x) * sin(y);
    return dvdx - dudy;
  };
  auto analyticalLaplacian = [](double x, double y) -> std::tuple<double, double> {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return {-4 * c1 * cos(2 * x) * cos(y) * cos(y) * sin(y), -4 * c2 * cos(x) * sin(y) * cos(y)};
  };
  for(const auto& e : mesh.edges()) {
    auto [px, py] = EdgeMidpoint(e);
    auto [u, v] = sphericalHarmonic(px, py);
    auto [lu, lv] = analyticalLaplacian(px, py);
    vec(e.get(), level) = primal_normal_x(e, level) * u + primal_normal_y(e, level) * v;
    lapVecSol(e.get(), level) = primal_normal_x(e, level) * lu + primal_normal_y(e, level) * lv;
  }
  for(const auto& v : mesh.vertices()) {
    rotVecSol(v, level) = analyticalCurl(v.x(), v.y());
  }
  for(const auto& c : mesh.faces()) {
    auto [cx, cy] = CellCircumcenter(c);
    divVecSol(c, level) = analyticalDivergence(cx, cy);
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical factors (sparse dimensions)
  //===------------------------------------------------------------------------------------------===//
  auto dot = [](const mylib::Vertex& v1, const mylib::Vertex& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y();
  };

  for(const auto& v : mesh.vertices()) {
    int m_sparse = 0;
    if(v.edges().size() != 6) {
      continue;
    }
    for(const auto& e : v.edges()) {
      mylib::Vertex testVec =
          mylib::Vertex(v.vertex(m_sparse).x() - v.x(), v.vertex(m_sparse).y() - v.y(), -1);
      mylib::Vertex dual = mylib::Vertex(dual_normal_x(*e, level), dual_normal_y(*e, level), -1);
      edge_orientation_vertex(v, m_sparse, level) = sgn(dot(testVec, dual));
      m_sparse++;
    }
  }
  for(const auto& c : mesh.faces()) {
    int m_sparse = 0;
    auto [xm, ym] = CellCircumcenter(c);
    for(const auto& e : c.edges()) {
      mylib::Vertex vOutside(e->vertex(0).x() - xm, e->vertex(0).y() - ym, -1);
      edge_orientation_cell(c, m_sparse, level) =
          sgn(dot(mylib::Vertex(e->vertex(0).x() - xm, e->vertex(0).y() - ym, -1),
                  mylib::Vertex(primal_normal_x(*e, level), primal_normal_y(*e, level), -1)));
      m_sparse++;
    }
  }

  // init sparse quantities for div and rot
  for(const auto& v : mesh.vertices()) {
    int m_sparse = 0;
    for(const auto& e : v.edges()) {
      geofac_rot(v, m_sparse, level) = dual_edge_length(*e, level) *
                                       edge_orientation_vertex(v, m_sparse, level) /
                                       dual_cell_area(v, level);
      m_sparse++;
    }
  }

  for(const auto& c : mesh.faces()) {
    int m_sparse = 0;
    for(const auto& e : c.edges()) {
      geofac_div(c, m_sparse, level) = primal_edge_length(*e, level) *
                                       edge_orientation_cell(c, m_sparse, level) /
                                       cell_area(c, level);
      m_sparse++;
    }
  }

  if(dbg_out) {
    dumpSparseData(mesh, geofac_rot, level, edgesPerVertex,
                   std::string("laplICONmylib_geofacRot.txt"));
    dumpSparseData(mesh, geofac_div, level, edgesPerCell,
                   std::string("laplICONmylib_geofacDiv.txt"));
  }

  //===------------------------------------------------------------------------------------------===//
  // stencil call
  //===------------------------------------------------------------------------------------------===//
  dawn_generated::cxxnaiveico::icon<mylibInterface::mylibTag>(
      mesh, k_size, vec, div_vec, rot_vec, nabla2t1_vec, nabla2t2_vec, nabla2_vec,
      primal_edge_length, dual_edge_length, tangent_orientation, geofac_rot, geofac_div)
      .run();

  if(dbg_out) {
    dumpField("laplICONmylib_nabla2t1.txt", mesh, nabla2t1_vec, level);
    dumpField("laplICONmylib_nabla2t2.txt", mesh, nabla2t2_vec, level);
  }

  //===------------------------------------------------------------------------------------------===//
  // report errors
  //===------------------------------------------------------------------------------------------===//

  {
    auto [Linf, L1, L2] = MeasureError(div_vec, divVecSol, innerCells(mesh), level);
    printf("[div] dx: %e L_inf: %e L_1: %e L_2: %e\n", 180. / w, Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureError(rot_vec, rotVecSol, innerNodes(mesh), level);
    printf("[rot] dx: %e L_inf: %e L_1: %e L_2: %e\n", 180. / w, Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureError(nabla2_vec, lapVecSol, innerEdges(mesh), level);
    printf("[lap] dx: %e L_inf: %e L_1: %e L_2: %e\n", 180. / w, Linf, L1, L2);
  }

  printf("-----\n");

  //===------------------------------------------------------------------------------------------===//
  // dumping a hopefully nice colorful divergence, curl and Laplacian
  //===------------------------------------------------------------------------------------------===//
  dumpField("laplICONmylib_div.txt", mesh, div_vec, level);
  dumpField("laplICONmylib_rot.txt", mesh, rot_vec, level);
  dumpField("laplICONmylib_out.txt", mesh, nabla2_vec, level);
}

namespace {

double EdgeLength(const mylib::Edge& e) {
  double x0 = e.vertex(0).x();
  double y0 = e.vertex(0).y();
  double x1 = e.vertex(1).x();
  double y1 = e.vertex(1).y();
  double dx = x1 - x0;
  double dy = y1 - y0;
  return sqrt(dx * dx + dy * dy);
}

std::tuple<double, double> EdgeMidpoint(const mylib::Edge& e) {
  double x0 = e.vertex(0).x();
  double y0 = e.vertex(0).y();
  double x1 = e.vertex(1).x();
  double y1 = e.vertex(1).y();
  return {0.5 * (x0 + x1), 0.5 * (y0 + y1)};
}

std::tuple<double, double> CellCircumcenter(const mylib::Face& c) {
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

std::tuple<double, double> PrimalNormal(const mylib::Edge& e) {
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

std::tuple<double, double> CellMidPoint(const mylib::Face& c) {
  auto v0 = c.vertex(0);
  auto v1 = c.vertex(1);
  auto v2 = c.vertex(2);
  return {1. / 3. * (v0.x() + v1.x() + v2.x()), 1. / 3. * (v0.y() + v1.y() + v2.y())};
}

double TriangleArea(const mylib::Vertex& v0, const mylib::Vertex& v1, const mylib::Vertex& v2) {
  return fabs(
      (v0.x() * (v1.y() - v2.y()) + v1.x() * (v2.y() - v0.y()) + v2.x() * (v0.y() - v1.y())) * 0.5);
}

double CellArea(const mylib::Face& c) {
  auto v0 = c.vertex(0);
  auto v1 = c.vertex(1);
  auto v2 = c.vertex(2);
  return TriangleArea(v0, v1, v2);
}

double DualCellArea(const mylib::Vertex& center) {
  double totalArea = 0.;
  for(const auto& e : center.edges()) {
    if(e->faces().size() != 2) {
      return 0.;
    }
    auto [leftx, lefty] = CellCircumcenter(e->face(0));
    auto [rightx, righty] = CellCircumcenter(e->face(1));
    mylib::Vertex left(leftx, lefty, -1);
    mylib::Vertex right(rightx, righty, -1);
    totalArea += TriangleArea(center, left, right);
  }
  return totalArea;
}

double DualEdgeLength(const mylib::Edge& e) {
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

double TangentOrientation(const mylib::Edge& e) {
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

std::vector<mylib::Face> innerCells(const mylib::Grid& m) {
  std::vector<mylib::Face> innerCells;
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
std::vector<mylib::Edge> innerEdges(const mylib::Grid& m) {
  std::vector<mylib::Edge> innerEdges;
  for(const auto e : m.edges()) {
    if(e.get().faces().size() != 2) {
      continue;
    }
    innerEdges.push_back(e);
  }
  return innerEdges;
}
std::vector<mylib::Vertex> innerNodes(const mylib::Grid& m) {
  std::vector<mylib::Vertex> innerVertices;
  for(const auto v : m.vertices()) {
    if(v.edges().size() != 6) {
      continue;
    }
    innerVertices.push_back(v);
  }
  return innerVertices;
}

void debugDumpMesh(const mylib::Grid& mesh, const std::string prefix) {
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
      fprintf(fp, "%f %f \n", nodeIt.x(), nodeIt.y());
    }
    fclose(fp);
  }
}

void dumpMesh(const mylib::Grid& m, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(const auto& e : m.edges()) {
    fprintf(fp, "%f %f %f %f\n", e.get().vertex(0).x(), e.get().vertex(0).y(),
            e.get().vertex(1).x(), e.get().vertex(1).y());
  }
  fclose(fp);
}

void dumpDualMesh(const mylib::Grid& m, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(const auto& e : m.edges()) {
    if(e.get().faces().size() != 2) {
      continue;
    }

    // This is WRONG!, leads to a dual mesh which is not orthogonal to
    // primal mesh
    // auto [xm1, ym1] = CellMidPoint(e.get().face(0));
    // auto [xm2, ym2] = CellMidPoint(e.get().face(1));

    auto [xm1, ym1] = CellCircumcenter(e.get().face(0));
    auto [xm2, ym2] = CellCircumcenter(e.get().face(1));
    fprintf(fp, "%f %f %f %f\n", xm1, ym1, xm2, ym2);
  }
  fclose(fp);
}

void dumpSparseData(const mylib::Grid& mesh, const mylib::SparseVertexData<double>& sparseData,
                    int level, int edgesPerVertex, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(const auto& v : mesh.vertices()) {
    int sparse_idx = 0;
    for(const auto& e : v.edges()) {
      double val = sparseData(v, sparse_idx, level);
      auto [emx, emy] = EdgeMidpoint(*e);
      double dx = emx - v.x();
      double dy = emy - v.y();
      fprintf(fp, "%f %f %f\n", v.x() + 0.5 * dx, v.y() + 0.5 * dy, val);
      sparse_idx++;
    }
  }
}

void dumpSparseData(const mylib::Grid& mesh, const mylib::SparseFaceData<double>& sparseData,
                    int level, int edgesPerCell, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(const auto& c : mesh.faces()) {
    int sparse_idx = 0;
    auto [cmx, cmy] = CellCircumcenter(c);
    for(const auto& e : c.edges()) {
      double val = sparseData(c, sparse_idx, level);
      auto [emx, emy] = EdgeMidpoint(*e);
      double dx = emx - cmx;
      double dy = emy - cmy;
      fprintf(fp, "%f %f %f\n", cmx + 0.5 * dx, cmy + 0.5 * dy, val);
      sparse_idx++;
    }
  }
}

void dumpField(const std::string& fname, const mylib::Grid& mesh,
               const mylib::EdgeData<double>& field, int level,
               std::optional<mylib::edge_color> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(auto& e : mesh.edges()) {
    if(color.has_value() && e.get().color() != color.value()) {
      continue;
    }
    auto [x, y] = EdgeMidpoint(e);
    fprintf(fp, "%f %f %f\n", x, y, std::isfinite(field(e, level)) ? field(e, level) : 0.);
  }
  fclose(fp);
}

void dumpField(const std::string& fname, const mylib::Grid& mesh,
               const mylib::EdgeData<double>& field_x, const mylib::EdgeData<double>& field_y,
               int level, std::optional<mylib::edge_color> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(auto& e : mesh.edges()) {
    if(color.has_value() && e.get().color() != color.value()) {
      continue;
    }
    auto [x, y] = EdgeMidpoint(e);
    fprintf(fp, "%f %f %f %f\n", x, y, field_x(e, level), field_y(e, level));
  }
  fclose(fp);
}

void dumpField(const std::string& fname, const mylib::Grid& mesh,
               const mylib::FaceData<double>& field, int level,
               std::optional<mylib::face_color> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(auto& c : mesh.faces()) {
    if(color.has_value() && c.color() != color.value()) {
      continue;
    }
    auto [x, y] = CellCircumcenter(c);
    fprintf(fp, "%f %f %f\n", x, y, field(c, level));
  }
  fclose(fp);
}

void dumpField(const std::string& fname, const mylib::Grid& mesh,
               const mylib::VertexData<double>& field, int level) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(auto& v : mesh.vertices()) {
    fprintf(fp, "%f %f %f\n", v.x(), v.y(), field(v, level));
  }
  fclose(fp);
}

template <typename DataT, typename ElemT>
std::tuple<double, double, double> MeasureError(DataT ref, DataT sol,
                                                std::vector<ElemT> innerElements, int level) {
  double Linf = 0.;
  double L1 = 0.;
  double L2 = 0.;
  for(const auto& elem : innerElements) {
    double dif = ref(elem, level) - sol(elem, level);
    if(!std::isfinite(dif)) {
      continue;
    }
    Linf = fmax(fabs(dif), Linf);
    L1 += fabs(dif);
    L2 += dif * dif;
  }
  L1 /= innerElements.size();
  L2 = sqrt(L2) / sqrt(innerElements.size());
  return {Linf, L1, L2};
}
} // namespace