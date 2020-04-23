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

#include "interfaces/toylib_interface.hpp"
#include "toylib.hpp"

#include "generated_iconLaplace.hpp"

#include "GenerateRectToylibMesh.h"

#include "io/toylibIO.h"

namespace {
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

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

  toylib::Grid mesh = toylibMeshRect(w);

  mesh.scale(M_PI / 180.); // rescale to radians

  if(dbg_out) {
    debugDumpMesh(mesh, "laplICONtoylib_Mesh");
  }

  const int edgesPerVertex = 6;
  const int edgesPerCell = 3;

  //===------------------------------------------------------------------------------------------===//
  // input field (field we want to take the laplacian of)
  //===------------------------------------------------------------------------------------------===//
  toylib::EdgeData<double> vec(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // control fields (containing analytical solutions)
  //===------------------------------------------------------------------------------------------===//
  toylib::FaceData<double> divVecSol(mesh, k_size);
  toylib::VertexData<double> rotVecSol(mesh, k_size);
  toylib::EdgeData<double> lapVecSol(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // output field (field containing the computed laplacian)
  //===------------------------------------------------------------------------------------------===//
  toylib::EdgeData<double> nabla2_vec(mesh, k_size);
  // term 1 and term 2 of nabla for debugging
  toylib::EdgeData<double> nabla2t1_vec(mesh, k_size);
  toylib::EdgeData<double> nabla2t2_vec(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // intermediary fields (curl/rot and div of vec_e)
  //===------------------------------------------------------------------------------------------===//
  toylib::VertexData<double> rot_vec(mesh, k_size);
  toylib::FaceData<double> div_vec(mesh, k_size);

  //===------------------------------------------------------------------------------------------===//
  // sparse dimensions for computing intermediary fields
  //===------------------------------------------------------------------------------------------===//
  toylib::SparseVertexData<double> geofac_rot(mesh, k_size, edgesPerVertex);
  toylib::SparseVertexData<double> edge_orientation_vertex(mesh, k_size, edgesPerVertex);

  toylib::SparseFaceData<double> geofac_div(mesh, k_size, edgesPerCell);
  toylib::SparseFaceData<double> edge_orientation_cell(mesh, k_size, edgesPerCell);

  //===------------------------------------------------------------------------------------------===//
  // fields containing geometric information
  //===------------------------------------------------------------------------------------------===//
  toylib::EdgeData<double> tangent_orientation(mesh, k_size);
  toylib::EdgeData<double> primal_edge_length(mesh, k_size);
  toylib::EdgeData<double> dual_edge_length(mesh, k_size);
  toylib::EdgeData<double> dual_normal_x(mesh, k_size);
  toylib::EdgeData<double> dual_normal_y(mesh, k_size);
  toylib::EdgeData<double> primal_normal_x(mesh, k_size);
  toylib::EdgeData<double> primal_normal_y(mesh, k_size);

  toylib::FaceData<double> cell_area(mesh, k_size);
  toylib::VertexData<double> dual_cell_area(mesh, k_size);

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
    dumpField("laplICONtoylib_tangentOrientation.txt", mesh, tangent_orientation, level);
    dumpField("laplICONtoylib_EdgeLength.txt", mesh, primal_edge_length, level);
    dumpField("laplICONtoylib_dualEdgeLength.txt", mesh, dual_edge_length, level);
    dumpField("laplICONtoylib_nrm.txt", mesh, primal_normal_x, primal_normal_y, level);
    dumpField("laplICONtoylib_dnrm.txt", mesh, dual_normal_x, dual_normal_y, level);
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
    dumpField("laplICONtoylib_areaCell.txt", mesh, cell_area, level);
    dumpField("laplICONtoylib_areaCellDual.txt", mesh, dual_cell_area, level);
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
  auto dot = [](const toylib::Vertex& v1, const toylib::Vertex& v2) {
    return v1.x() * v2.x() + v1.y() * v2.y();
  };

  for(const auto& v : mesh.vertices()) {
    int m_sparse = 0;
    if(v.edges().size() != 6) {
      continue;
    }
    for(const auto& e : v.edges()) {
      toylib::Vertex testVec =
          toylib::Vertex(v.vertex(m_sparse).x() - v.x(), v.vertex(m_sparse).y() - v.y(), -1);
      toylib::Vertex dual = toylib::Vertex(dual_normal_x(*e, level), dual_normal_y(*e, level), -1);
      edge_orientation_vertex(v, m_sparse, level) = sgn(dot(testVec, dual));
      m_sparse++;
    }
  }
  for(const auto& c : mesh.faces()) {
    int m_sparse = 0;
    auto [xm, ym] = CellCircumcenter(c);
    for(const auto& e : c.edges()) {
      toylib::Vertex vOutside(e->vertex(0).x() - xm, e->vertex(0).y() - ym, -1);
      edge_orientation_cell(c, m_sparse, level) =
          sgn(dot(toylib::Vertex(e->vertex(0).x() - xm, e->vertex(0).y() - ym, -1),
                  toylib::Vertex(primal_normal_x(*e, level), primal_normal_y(*e, level), -1)));
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
                   std::string("laplICONtoylib_geofacRot.txt"));
    dumpSparseData(mesh, geofac_div, level, edgesPerCell,
                   std::string("laplICONtoylib_geofacDiv.txt"));
  }

  //===------------------------------------------------------------------------------------------===//
  // stencil call
  //===------------------------------------------------------------------------------------------===//
  dawn_generated::cxxnaiveico::ICON_laplacian_stencil<toylibInterface::toylibTag>(
      mesh, k_size, vec, div_vec, rot_vec, nabla2t1_vec, nabla2t2_vec, nabla2_vec,
      primal_edge_length, dual_edge_length, tangent_orientation, geofac_rot, geofac_div)
      .run();

  if(dbg_out) {
    dumpField("laplICONtoylib_nabla2t1.txt", mesh, nabla2t1_vec, level);
    dumpField("laplICONtoylib_nabla2t2.txt", mesh, nabla2t2_vec, level);
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
  dumpField("laplICONtoylib_div.txt", mesh, div_vec, level);
  dumpField("laplICONtoylib_rot.txt", mesh, rot_vec, level);
  dumpField("laplICONtoylib_out.txt", mesh, nabla2_vec, level);
}

namespace {

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