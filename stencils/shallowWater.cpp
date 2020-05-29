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

//===------------------------------------------------------------------------------------------===//
//
//  UNTESTED, INCOMPLETE
//
//===------------------------------------------------------------------------------------------===//

// Things to checks
//  [X] why are the normals off?
//      => mesh was scaled with two slightly different factors along x and y
//  [X] introduce the stabilization term
//      => doesn't seem to help
//  - what about the boundaries?
//  - is the method really consistent now?

// Shallow water equation solver as described in "A simple and efficient unstructured finite volume
// scheme for solving the shallow water equations in overland flow applications" by Cea and Bladé
// Follows notation in the paper as closely as possilbe

#include <cmath>
#include <cstdio>
#include <fenv.h>
#include <optional>
#include <set>
#include <vector>

// atlas functions
#include <atlas/array.h>
#include <atlas/grid.h>
#include <atlas/mesh.h>
#include <atlas/mesh/actions/BuildEdges.h>
#include <atlas/meshgenerator.h>
#include <atlas/util/CoordinateEnums.h>

// atlas interface for dawn generated code
#include "interfaces/atlas_interface.hpp"

// icon stencil
#include "generated_iconLaplace.hpp"

// atlas utilities
#include "../utils/AtlasCartesianWrapper.h"
#include "../utils/GenerateRectAtlasMesh.h"

template <typename T>
static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

void dumpMesh4Triplot(const atlas::Mesh& mesh, const std::string prefix,
                      const atlasInterface::Field<double>& field,
                      std::optional<AtlasToCartesian> wrapper);

void dumpNodeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level);
void dumpCellField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level);
void dumpCellFieldOnNodes(const std::string& fname, const atlas::Mesh& mesh,
                          AtlasToCartesian wrapper, atlasInterface::Field<double>& field,
                          int level);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level,
                   std::optional<Orientation> color = std::nullopt);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level, std::vector<int> edgeList,
                   std::optional<Orientation> color = std::nullopt);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field_x, atlasInterface::Field<double>& field_y,
                   int level, std::optional<Orientation> color = std::nullopt);

void gradient(const atlas::Mesh& mesh, const atlasInterface::Field<double>& f,
              atlasInterface::Field<double>& nx, atlasInterface::Field<double>& ny,
              atlasInterface::SparseDimension<double>& edge_orientation_cell,
              const atlasInterface::Field<double>& L, const atlasInterface::Field<double>& A,
              atlasInterface::Field<double>& f_x, atlasInterface::Field<double>& f_y, int level) {
  const auto& conn = mesh.cells().edge_connectivity();
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    double lhs_x = 0.;
    double lhs_y = 0.;
    for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
      int edgeIdx = conn(cellIdx, nbhIdx);
      lhs_x += f(edgeIdx, level) * nx(edgeIdx, level) *
               edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
      lhs_y += f(edgeIdx, level) * ny(edgeIdx, level) *
               edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
    }
    f_x(cellIdx, level) = lhs_x / A(cellIdx, level);
    f_y(cellIdx, level) = lhs_y / A(cellIdx, level);

    // if(fabs(lhs_x / A(cellIdx, level)) > 1e-6) {
    //   printf("!\n");
    // }
  }
}

void divergence(const atlas::Mesh& mesh, const atlasInterface::Field<double>& u,
                const atlasInterface::Field<double>& v, atlasInterface::Field<double>& nx,
                atlasInterface::Field<double>& ny,
                atlasInterface::SparseDimension<double>& edge_orientation_cell,
                const atlasInterface::Field<double>& L, const atlasInterface::Field<double>& A,
                atlasInterface::Field<double>& div, int level) {
  const auto& conn = mesh.cells().edge_connectivity();
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    double lhs = 0.;
    for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
      int edgeIdx = conn(cellIdx, nbhIdx);
      lhs += (u(edgeIdx, level) * nx(edgeIdx, level) + v(edgeIdx, level) * ny(edgeIdx, level)) *
             edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
    }
    div(cellIdx, level) = lhs / A(cellIdx, level);
  }
}

void laplacianDiamond(const atlas::Mesh& mesh, const atlasInterface::Field<double>& fV,
                      const atlasInterface::Field<double>& fE,
                      const atlasInterface::Field<double>& L,
                      const atlasInterface::Field<double>& vertVertL,
                      atlasInterface::Field<double>& lapl, int level) {

  static std::vector<int>* diamondNBh = 0;
  if(diamondNBh == 0) {
    diamondNBh = new std::vector<int>(mesh->edges().size() * 4);
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      auto diamondNbhE = atlasInterface::getNeighbors(
          atlasInterface::atlasTag{}, mesh,
          {dawn::LocationType::Edges, dawn::LocationType::Cells, dawn::LocationType::Vertices},
          edgeIdx);

      if(diamondNbhE.size() == 4) {
        (*diamondNBh)[4 * edgeIdx + 0] = diamondNbhE[0];
        (*diamondNBh)[4 * edgeIdx + 1] = diamondNbhE[1];
        (*diamondNBh)[4 * edgeIdx + 2] = diamondNbhE[2];
        (*diamondNBh)[4 * edgeIdx + 3] = diamondNbhE[3];
      } else {
        (*diamondNBh)[4 * edgeIdx + 0] = -1;
        (*diamondNBh)[4 * edgeIdx + 1] = -1;
        (*diamondNBh)[4 * edgeIdx + 2] = -1;
        (*diamondNBh)[4 * edgeIdx + 3] = -1;
      }
    }
  }

  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    std::vector<int> diamondNbhE = {(*diamondNBh)[4 * edgeIdx + 0], (*diamondNBh)[4 * edgeIdx + 1],
                                    (*diamondNBh)[4 * edgeIdx + 2], (*diamondNBh)[4 * edgeIdx + 3]};

    if(std::find(diamondNbhE.begin(), diamondNbhE.end(), -1) != diamondNbhE.end()) {
      lapl(edgeIdx, level) = 0.;
      continue;
    }

    double hx = 0.5 * L(mesh, edgeIdx);
    double hy = 0.5 * vertVertL(mesh, edgeIdx);

    double fxx =
        (fV(diamondNbhE[0], level) + fV(diamondNbhE[1], level) - 2 * (fE(edgeIdx, level))) /
        (hx * hx);
    double fyy =
        (fV(diamondNbhE[2], level) + fV(diamondNbhE[3], level) - 2 * (fE(edgeIdx, level))) /
        (hy * hy);

    lapl(edgeIdx, level) = fxx + fyy;
  }
}

void lerp(const atlas::Mesh& mesh, const atlasInterface::Field<double>& src,
          atlasInterface::Field<double>& dest, atlasInterface::Field<double>& alpha, int level) {

  const auto& conn = mesh.edges().cell_connectivity();
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    double lhs = 0.;
    double weights[2] = {1 - alpha(edgeIdx, level), alpha(edgeIdx, level)};
    for(int nbhIdx = 0; nbhIdx < conn.cols(edgeIdx); nbhIdx++) {
      int cellIdx = conn(edgeIdx, nbhIdx);
      if(cellIdx == conn.missing_value()) {
        assert(weights[nbhIdx] == 0.);
        continue;
      }
      lhs += src(cellIdx, level) * weights[nbhIdx];
    }
    dest(edgeIdx, level) = lhs;
  }
}

void lerpe2v(const atlas::Mesh& mesh, const atlasInterface::Field<double>& src,
             atlasInterface::Field<double>& dest, int level) {

  const auto& conn = mesh.nodes().edge_connectivity();
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    double lhs = 0.;
    for(int nbhIdx = 0; nbhIdx < conn.cols(nodeIdx); nbhIdx++) {
      int edgeIdx = conn(nodeIdx, nbhIdx);
      if(edgeIdx == conn.missing_value()) {
        continue;
      }
      lhs += src(edgeIdx, level);
    }
    dest(nodeIdx, level) = lhs / conn.cols(nodeIdx);
  }
}

int main(int argc, char const* argv[]) {

  // enable floating point exception
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  if(argc != 2) {
    std::cout << "intended use is\n" << argv[0] << " ny" << std::endl;
    return -1;
  }
  int w = atoi(argv[1]);
  // reference level of fluid, make sure to chose this large enough, otherwise initial
  // splash may induce negative fluid height and crash the sim
  const double refHeight = 2.;

  // constants
  const double CFLconst = 0.05;
  const double Grav = -9.81;

  int k_size = 1;
  const int level = 0;
  double lDomain = 10;

  bool useFriction = true;
  double fricCoeff = 0.3;
  double densityWater = 1e3;

  bool useDamping = true;
  double dampingCoeff = 0.5;

  // dump a whole bunch of debug output (meant to be visualized using Octave, but gnuplot and the
  // like will certainly work too)
  const bool dbg_out = false;
  const bool readMeshFromDisk = false;

  auto mesh = AtlasMeshRect(w * sqrt(3), w);
  atlas::mesh::actions::build_edges(mesh, atlas::util::Config("pole_edges", false));
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
  atlas::mesh::actions::build_element_to_edge_connectivity(mesh);

  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    const auto& nodeToEdge = mesh.nodes().edge_connectivity();
    const auto& edgeToCell = mesh.edges().cell_connectivity();
    auto& nodeToCell = mesh.nodes().cell_connectivity();

    std::set<int> nbh;
    for(int nbhEdgeIdx = 0; nbhEdgeIdx < nodeToEdge.cols(nodeIdx); nbhEdgeIdx++) {
      int edgeIdx = nodeToEdge(nodeIdx, nbhEdgeIdx);
      if(edgeIdx == nodeToEdge.missing_value()) {
        continue;
      }
      for(int nbhCellIdx = 0; nbhCellIdx < edgeToCell.cols(edgeIdx); nbhCellIdx++) {
        int cellIdx = edgeToCell(edgeIdx, nbhCellIdx);
        if(cellIdx == edgeToCell.missing_value()) {
          continue;
        }
        nbh.insert(cellIdx);
      }
    }

    assert(nbh.size() <= 6);
    std::vector<int> initData(nbh.size(), nodeToCell.missing_value());
    nodeToCell.add(1, nbh.size(), initData.data());
    int copyIter = 0;
    for(const int n : nbh) {
      nodeToCell.set(nodeIdx, copyIter++, n);
    }
  }

  printf("mesh stats: #cells %d #nodes %d #edges %d\n", mesh->cells().size(), mesh->nodes().size(),
         mesh->edges().size());

  // wrapper with various atlas helper functions
  AtlasToCartesian wrapper(mesh, lDomain, false, true);

  const int edgesPerVertex = 6;
  const int edgesPerCell = 3;

  //===------------------------------------------------------------------------------------------===//
  // helper lambdas to readily construct atlas fields and views on one line
  //===------------------------------------------------------------------------------------------===//
  auto MakeAtlasField = [&](const std::string& name,
                            int size) -> std::tuple<atlas::Field, atlasInterface::Field<double>> {
    atlas::Field field_F{name, atlas::array::DataType::real64(),
                         atlas::array::make_shape(mesh.edges().size(), k_size)};
    return {field_F, atlas::array::make_view<double, 2>(field_F)};
  };

  auto MakeAtlasSparseField =
      [&](const std::string& name, int size,
          int sparseSize) -> std::tuple<atlas::Field, atlasInterface::SparseDimension<double>> {
    atlas::Field field_F{name, atlas::array::DataType::real64(),
                         atlas::array::make_shape(mesh.edges().size(), k_size, sparseSize)};
    return {field_F, atlas::array::make_view<double, 3>(field_F)};
  };

  auto [uC_F, uC] = MakeAtlasField("u", mesh.cells().size());
  auto [vC_F, vC] = MakeAtlasField("v", mesh.cells().size());
  auto [hC_F, hC] = MakeAtlasField("h", mesh.cells().size());

  auto [uE_F, uE] = MakeAtlasField("u", mesh.edges().size());
  auto [vE_F, vE] = MakeAtlasField("v", mesh.edges().size());
  auto [hE_F, hE] = MakeAtlasField("h", mesh.edges().size());

  auto [laplhE_F, laplhE] = MakeAtlasField("laplh", mesh.edges().size());

  auto [hV_F, hV] = MakeAtlasField("h", mesh.nodes().size());

  auto [h_xC_F, h_xC] = MakeAtlasField("h_x", mesh.cells().size());
  auto [h_yC_F, h_yC] = MakeAtlasField("h_y", mesh.cells().size());
  auto [divuvC_F, divuvC] = MakeAtlasField("divuv", mesh.cells().size());

  auto [u_tC_F, u_tC] = MakeAtlasField("u_t", mesh.cells().size());
  auto [v_tC_F, v_tC] = MakeAtlasField("v_t", mesh.cells().size());
  auto [h_tC_F, h_tC] = MakeAtlasField("h_t", mesh.cells().size());

  auto [uinitC_F, uinitC] = MakeAtlasField("u", mesh.cells().size());
  auto [vinitC_F, vinitC] = MakeAtlasField("v", mesh.cells().size());
  auto [hinitC_F, hinitC] = MakeAtlasField("h", mesh.cells().size());

  // CFL per cell
  auto [cfl_F, cfl] = MakeAtlasField("CFL", mesh.cells().size());

  // Geometrical factors on edges
  auto [lambda_F, lambda] = MakeAtlasField("lambda", mesh.edges().size());  // normal velocity
  auto [L_F, L] = MakeAtlasField("L", mesh.edges().size());                 // edge length
  auto [vertVertL_F, vertVertL] = MakeAtlasField("L", mesh.edges().size()); // edge length
  auto [nx_F, nx] = MakeAtlasField("nx", mesh.edges().size());              // normals
  auto [ny_F, ny] = MakeAtlasField("ny", mesh.edges().size());
  auto [alpha_F, alpha] = MakeAtlasField("alpha", mesh.edges().size());

  // Geometrical factors on cells
  auto [A_F, A] = MakeAtlasField("A", mesh.cells().size());
  auto [edge_orientation_cell_F, edge_orientation_cell] =
      MakeAtlasSparseField("edge_orientation_cell", mesh.cells().size(), edgesPerCell);

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    L(edgeIdx, level) = wrapper.edgeLength(mesh, edgeIdx);
    auto [nxe, nye] = wrapper.primalNormal(mesh, edgeIdx);
    nx(edgeIdx, level) = nxe;
    ny(edgeIdx, level) = nye;

    auto diamondNbh = atlasInterface::getNeighbors(
        atlasInterface::atlasTag{}, mesh,
        {dawn::LocationType::Edges, dawn::LocationType::Cells, dawn::LocationType::Vertices},
        edgeIdx);

    auto [x1, y1] = wrapper.nodeLocation(diamondNbh[2]);
    auto [x2, y2] = wrapper.nodeLocation(diamondNbh[3]);

    double dx = x1 - x2;
    double dy = y1 - y2;

    vertVertL(edgeIdx, level) = sqrt(dx * dx + dy * dy);
  }

  {
    const auto& conn = mesh.edges().cell_connectivity();
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      int cellIdx1 = conn(edgeIdx, 0);
      int cellIdx2 = conn(edgeIdx, 1);
      double d1 = (cellIdx1 >= 0) ? wrapper.distanceToCircumcenter(mesh, cellIdx1, edgeIdx) : 0.;
      double d2 = (cellIdx2 >= 0) ? wrapper.distanceToCircumcenter(mesh, cellIdx2, edgeIdx) : 0.;
      alpha(edgeIdx, level) = d2 / (d1 + d2);
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on cells
  //===------------------------------------------------------------------------------------------===//
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    A(cellIdx, level) = wrapper.cellArea(mesh, cellIdx);
  }

  auto dot = [](const Vector& v1, const Vector& v2) {
    return std::get<0>(v1) * std::get<0>(v2) + std::get<1>(v1) * std::get<1>(v2);
  };
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    const atlas::mesh::HybridElements::Connectivity& cellEdgeConnectivity =
        mesh.cells().edge_connectivity();
    auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);

    const int missingVal = cellEdgeConnectivity.missing_value();
    int numNbh = cellEdgeConnectivity.cols(cellIdx);
    assert(numNbh == edgesPerCell);

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = cellEdgeConnectivity(cellIdx, nbhIdx);
      auto [emX, emY] = wrapper.edgeMidpoint(mesh, edgeIdx);
      Vector toOutsdie{emX - xm, emY - ym};
      Vector primal = {nx(edgeIdx, level), ny(edgeIdx, level)};
      edge_orientation_cell(cellIdx, nbhIdx, level) = sgn(dot(toOutsdie, primal));
    }
    // explanation: the vector cellMidpoint -> edgeMidpoint is guaranteed to point outside. The
    // dot product checks if the edge normal has the same orientation. edgeMidpoint is arbitrary,
    // any point on e would work just as well
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize height and other fields
  //===------------------------------------------------------------------------------------------===//
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
    xm -= 0;
    ym -= 0;
    double v = sqrt(xm * xm + ym * ym);
    hC(cellIdx, level) = exp(-5 * v * v) + refHeight;
    // hC(cellIdx, level) = refHeight;
    // h(cellIdx, level) = sin(xm) * sin(ym) + refHeight;
  }

  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    uC(cellIdx, level) = 0.;
    vC(cellIdx, level) = 0.;
  }

  //===------------------------------------------------------------------------------------------===//
  // collect boundary edgse
  //===------------------------------------------------------------------------------------------===//
  std::set<int> boundaryEdgesSet;
  {
    const auto& conn = mesh.edges().cell_connectivity();
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      if(conn.cols(edgeIdx) < 2) {
        boundaryEdgesSet.insert(edgeIdx);
        continue;
      }

      int cellIdx0 = conn(edgeIdx, 0);
      int cellIdx1 = conn(edgeIdx, 1);
      if(cellIdx0 == conn.missing_value() || cellIdx1 == conn.missing_value()) {
        boundaryEdgesSet.insert(edgeIdx);
      }
    }
  }
  std::vector<int> boundaryEdges(boundaryEdgesSet.begin(), boundaryEdgesSet.end());
  std::vector<double> boundary_Nx, boundary_Ny;

  std::vector<int> boundaryCells;
  {
    const auto& conn = mesh.cells().edge_connectivity();
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      int edgeIdx0 = conn(cellIdx, 0);
      int edgeIdx1 = conn(cellIdx, 1);
      int edgeIdx2 = conn(cellIdx, 2);

      if(boundaryEdgesSet.count(edgeIdx0) || boundaryEdgesSet.count(edgeIdx1) ||
         boundaryEdgesSet.count(edgeIdx2)) {
        boundaryCells.push_back(cellIdx);
      }
    }
  }

  {
    FILE* fp = fopen("boundaryCells.txt", "w+");
    for(auto it : boundaryCells) {
      auto [x, y] = wrapper.cellCircumcenter(mesh, it);
      fprintf(fp, "%f %f\n", x, y);
    }
    fclose(fp);
  }

  dumpEdgeField("nrm", mesh, wrapper, nx, ny, level);
  dumpEdgeField("L", mesh, wrapper, L, level);

  // dumpMesh4Triplot(mesh, "init", h, std::nullopt);
  dumpMesh4Triplot(mesh, "init", hC, wrapper);

  double t = 0.;
  double dt = 0;
  double t_final = 16.;
  int step = 0;

  // writing this intentionally close to generated code
  while(t < t_final) {
    // make some splashes
    if(step > 0 && step % 1000 == 0) {
      double splash_x = (rand() / ((double)RAND_MAX) - 0.5) * 2;
      double splash_y = (rand() / ((double)RAND_MAX) - 0.5) * 2;
      printf("splashing at %f %f\n", splash_x, splash_y);
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
        xm += splash_x;
        ym += splash_y;
        double v = sqrt(xm * xm + ym * ym);
        hC(cellIdx, level) += exp(-5 * v * v);
      }
    }

    // dumpCellField("uC", mesh, wrapper, uC, level);
    // dumpCellField("vC", mesh, wrapper, vC, level);
    // dumpCellField("hC", mesh, wrapper, hC, level);

    // init
    for(int cellIdx = 0; cellIdx < mesh->cells().size(); cellIdx++) {
      hinitC(cellIdx, level) = hC(cellIdx, level);
      uinitC(cellIdx, level) = uC(cellIdx, level);
      vinitC(cellIdx, level) = vC(cellIdx, level);
    }

    // predict
    for(int cellIdx = 0; cellIdx < mesh->cells().size(); cellIdx++) {
      hC(cellIdx, level) = hinitC(cellIdx, level) + 0.5 * dt * h_tC(cellIdx, level);
      uC(cellIdx, level) = uinitC(cellIdx, level) + 0.5 * dt * u_tC(cellIdx, level);
      vC(cellIdx, level) = vinitC(cellIdx, level) + 0.5 * dt * v_tC(cellIdx, level);
    }

    lerp(mesh, hC, hE, alpha, level);
    lerp(mesh, uC, uE, alpha, level);
    lerp(mesh, vC, vE, alpha, level);

    // dumpEdgeField("uE", mesh, wrapper, uE, level);
    // dumpEdgeField("vE", mesh, wrapper, vE, level);
    // dumpEdgeField("hE", mesh, wrapper, hE, level);

    for(auto bndIdx : boundaryEdges) {
      uE(bndIdx, level) = 0.;
      vE(bndIdx, level) = 0.;
    }

    // stabilizing(hi pass filtering)
    if(useDamping) {
      lerpe2v(mesh, hE, hV, level);
      laplacianDiamond(mesh, hV, hE, L, vertVertL, laplhE, level);
      for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
        hE(edgeIdx, level) += dampingCoeff * dt * laplhE(edgeIdx, level);
      }
    }

    gradient(mesh, hE, nx, ny, edge_orientation_cell, L, A, h_xC, h_yC, level);
    divergence(mesh, uE, vE, nx, ny, edge_orientation_cell, L, A, divuvC, level);

    for(auto bndIdx : boundaryCells) {
      h_xC(bndIdx, level) = 0.;
      h_yC(bndIdx, level) = 0.;
    }

    // dumpCellField("h_xC", mesh, wrapper, h_xC, level);
    // dumpCellField("h_yC", mesh, wrapper, h_yC, level);
    // dumpCellField("divuvC", mesh, wrapper, divuvC, level);

    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      double absVelRel =
          sqrt(uC(cellIdx, level) * uC(cellIdx, level) + vC(cellIdx, level) * vC(cellIdx, level));
      double fricx = 0.;
      double fricy = 0.;
      if(absVelRel > 1e-8 && useFriction) {
        double fricx = fricCoeff * hC(cellIdx, level) * A(cellIdx, level) * densityWater *
                       fabs(Grav) * uC(cellIdx, level) / absVelRel;
        double fricy = fricCoeff * hC(cellIdx, level) * A(cellIdx, level) * densityWater *
                       fabs(Grav) * uC(cellIdx, level) / absVelRel;
      }
      u_tC(cellIdx, level) = -Grav * h_xC(cellIdx, level) + fricx;
      v_tC(cellIdx, level) = -Grav * h_yC(cellIdx, level) + fricy;
      h_tC(cellIdx, level) = hC(cellIdx, level) * (divuvC(cellIdx, level));
    }

    // dumpCellField("u_tC", mesh, wrapper, u_tC, level);
    // dumpCellField("v_tC", mesh, wrapper, v_tC, level);
    // dumpCellField("h_tC", mesh, wrapper, h_tC, level);

    // correct
    for(int cellIdx = 0; cellIdx < mesh->cells().size(); cellIdx++) {
      hC(cellIdx, level) = hinitC(cellIdx, level) + dt * h_tC(cellIdx, level);
      uC(cellIdx, level) = uinitC(cellIdx, level) + dt * u_tC(cellIdx, level);
      vC(cellIdx, level) = vinitC(cellIdx, level) + dt * v_tC(cellIdx, level);
    }

    // adapt CLF
    // this would probably be in the driver code anyway
    {
      const auto& conn = mesh.cells().edge_connectivity();
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double l0 = L(conn(cellIdx, 0), level);
        double l1 = L(conn(cellIdx, 1), level);
        double l2 = L(conn(cellIdx, 2), level);
        double hi = hC(cellIdx, level);
        double Ux = uC(cellIdx, level) / hi;
        double Uy = vC(cellIdx, level) / hi;
        double U = sqrt(Ux * Ux + Uy * Uy);
        cfl(cellIdx, level) = CFLconst * std::min({l0, l1, l2}) / (U + sqrt(fabs(Grav) * hi));
      }
      double mindt = std::numeric_limits<double>::max();
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        mindt = fmin(cfl(cellIdx, level), mindt);
      }
      dt = mindt;
    }

    t += dt;

    if(step % 20 == 0) {
      char buf[256];
      sprintf(buf, "out/stepH_%04d.txt", step);
      dumpCellFieldOnNodes(buf, mesh, wrapper, hC, level);
    }
    std::cout << "time " << t << " timestep " << step++ << " dt " << dt << "\n";
  }

  dumpMesh4Triplot(mesh, "final", hC, wrapper);
}

void dumpMesh4Triplot(const atlas::Mesh& mesh, const std::string prefix,
                      const atlasInterface::Field<double>& field,
                      std::optional<AtlasToCartesian> wrapper) {
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
      if(wrapper == std::nullopt) {
        double x = xy(nodeIdx, atlas::LON);
        double y = xy(nodeIdx, atlas::LAT);
        fprintf(fp, "%f %f \n", x, y);
      } else {
        auto [x, y] = wrapper.value().nodeLocation(nodeIdx);
        fprintf(fp, "%f %f \n", x, y);
      }
    }
    fclose(fp);
  }

  {
    char buf[256];
    sprintf(buf, "%sC.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      fprintf(fp, "%f\n", field(cellIdx, 0));
    }
    fclose(fp);
  }
}

void dumpNodeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    auto [xm, ym] = wrapper.nodeLocation(nodeIdx);
    fprintf(fp, "%f %f %f\n", xm, ym, field(nodeIdx, level));
  }
  fclose(fp);
}

void dumpCellFieldOnNodes(const std::string& fname, const atlas::Mesh& mesh,
                          AtlasToCartesian wrapper, atlasInterface::Field<double>& field,
                          int level) {
  FILE* fp = fopen(fname.c_str(), "w+");
  const auto& conn = mesh.nodes().cell_connectivity();
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    double h = 0.;
    for(int nbhIdx = 0; nbhIdx < conn.cols(nodeIdx); nbhIdx++) {
      int cIdx = conn(nodeIdx, nbhIdx);
      h += field(cIdx, 0);
    }
    h /= conn.cols(nodeIdx);
    fprintf(fp, "%f\n", h);
  }
  fclose(fp);
}

void dumpCellField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
    fprintf(fp, "%f %f %e\n", xm, ym, field(cellIdx, level));
  }
  fclose(fp);
}

void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level,
                   std::optional<Orientation> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    if(color.has_value() && wrapper.edgeOrientation(mesh, edgeIdx) != color.value()) {
      continue;
    }
    auto [xm, ym] = wrapper.edgeMidpoint(mesh, edgeIdx);
    fprintf(fp, "%f %f %f\n", xm, ym,
            std::isfinite(field(edgeIdx, level)) ? field(edgeIdx, level) : 0.);
  }
  fclose(fp);
}

void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level, std::vector<int> edgeList,
                   std::optional<Orientation> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(int edgeIdx : edgeList) {
    if(color.has_value() && wrapper.edgeOrientation(mesh, edgeIdx) != color.value()) {
      continue;
    }
    auto [xm, ym] = wrapper.edgeMidpoint(mesh, edgeIdx);
    fprintf(fp, "%f %f %f\n", xm, ym,
            std::isfinite(field(edgeIdx, level)) ? field(edgeIdx, level) : 0.);
  }
  fclose(fp);
}

void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field_x, atlasInterface::Field<double>& field_y,
                   int level, std::optional<Orientation> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    if(color.has_value() && wrapper.edgeOrientation(mesh, edgeIdx) != color.value()) {
      continue;
    }
    auto [xm, ym] = wrapper.edgeMidpoint(mesh, edgeIdx);
    fprintf(fp, "%f %f %f %f\n", xm, ym, field_x(edgeIdx, level), field_y(edgeIdx, level));
  }
  fclose(fp);
}

//===------------------------------------------------------------------------------------------===//
// snippet to nudge boundary triangles to get domain which is strictly bound by a square / rectangle
//===------------------------------------------------------------------------------------------===//

// {
//   auto& xy = wrapper.xy();
//   const auto& conn = mesh->edges().node_connectivity();
//   for(auto boundaryEdge : boundaryEdges) {
//     auto [p1, p2] = wrapper.cartesianEdge(mesh, boundaryEdge);
//     if(fabs(std::get<1>(p1) - std::get<1>(p2)) > 1e-8) {
//       if(std::get<0>(p1) < 0) {
//         double x1 = std::get<0>(p1);
//         double x2 = std::get<0>(p2);
//         double xmin = fmin(x1, x2);
//         int i1 = conn(boundaryEdge, 0);
//         int i2 = conn(boundaryEdge, 1);
//         xy[conn(boundaryEdge, 0)] = {xmin, std::get<1>(xy[conn(boundaryEdge, 0)])};
//         xy[conn(boundaryEdge, 1)] = {xmin, std::get<1>(xy[conn(boundaryEdge, 1)])};
//       } else {
//         double x1 = std::get<0>(p1);
//         double x2 = std::get<0>(p2);
//         double xmax = fmax(x1, x2);
//         int i1 = conn(boundaryEdge, 0);
//         int i2 = conn(boundaryEdge, 1);
//         xy[conn(boundaryEdge, 0)] = {xmax, std::get<1>(xy[conn(boundaryEdge, 0)])};
//         xy[conn(boundaryEdge, 1)] = {xmax, std::get<1>(xy[conn(boundaryEdge, 1)])};
//       }
//     }
//   }
// }

//===------------------------------------------------------------------------------------------===//
// snippet to compute boundary normals (guaranteed to point outside)
//===------------------------------------------------------------------------------------------===//

// {
//   const auto& conn = mesh->edges().cell_connectivity();
//   for(auto boundaryEdge : boundaryEdges) {
//     auto [p1, p2] = wrapper.cartesianEdge(mesh, boundaryEdge);
//     double ny = -(std::get<0>(p1) - std::get<0>(p2));
//     double nx = std::get<1>(p1) - std::get<1>(p2);
//     nx /= sqrt(nx * nx + ny * ny);
//     ny /= sqrt(nx * nx + ny * ny);
//     int cellIdx = conn(boundaryEdge, 0) == conn.missing_value() ? conn(boundaryEdge, 1)
//                                                                 : conn(boundaryEdge, 0);
//     auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
//     auto [emX, emY] = wrapper.edgeMidpoint(mesh, boundaryEdge);
//     Vector toOutsdie{emX - xm, emY - ym};
//     Vector primal = {nx, ny};
//     int orientation = sgn(dot(toOutsdie, primal));
//     orientation == -1 ? boundary_Nx.push_back(-nx) : boundary_Nx.push_back(nx);
//     orientation == -1 ? boundary_Ny.push_back(-ny) : boundary_Ny.push_back(ny);
//   }
// }