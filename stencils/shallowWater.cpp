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

//===-----------------------------------------------------------------------------

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

  // use high frequency damping. original damping by Cea and Blade is heavily dissipative, hence the
  // damping can be modulated by a coefficient in this implementation
  const bool use_corrector = true;
  const double DampingCoeff = 0.01;

  // optional bed friction, manning coefficient of 0.01 is roughly equal to flow of water over
  // concrete
  const bool use_friction = true;
  const double ManningCoeff = 0.01;

  int k_size = 1;
  const int level = 0;
  double lDomain = 10;

  // dump a whole bunch of debug output (meant to be visualized using Octave, but gnuplot and the
  // like will certainly work too)
  const bool dbg_out = false;
  const bool readMeshFromDisk = false;

  auto mesh = AtlasMeshSquare(w);
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

  // Edge Fluxes
  auto [Q_F, Q] = MakeAtlasField("Q", mesh.edges().size());    // mass
  auto [Fx_F, Fx] = MakeAtlasField("Fx", mesh.edges().size()); // momentum
  auto [Fy_F, Fy] = MakeAtlasField("Fy", mesh.edges().size());

  // Edge Velocities (to be interpolated from cell circumcenters)
  auto [Ux_F, Ux] = MakeAtlasField("Ux", mesh.edges().size());
  auto [Uy_F, Uy] = MakeAtlasField("Uy", mesh.edges().size());

  // Height on edges (to be interpolated from cell circumcenters)
  auto [hs_F, hs] = MakeAtlasField("hs", mesh.edges().size());

  // Cell Centered Values
  auto [h_F, h] = MakeAtlasField("h", mesh.cells().size());    // fluid height
  auto [qx_F, qx] = MakeAtlasField("qx", mesh.cells().size()); // discharge
  auto [qy_F, qy] = MakeAtlasField("qy", mesh.cells().size());
  auto [Sx_F, Sx] = MakeAtlasField("Sx", mesh.cells().size()); // free surface gradient
  auto [Sy_F, Sy] = MakeAtlasField("Sy", mesh.cells().size());

  // Time Derivative of Cell Centered Values
  auto [dhdt_F, dhdt] = MakeAtlasField("h", mesh.cells().size());    // fluid height
  auto [dqxdt_F, dqxdt] = MakeAtlasField("qx", mesh.cells().size()); // discharge
  auto [dqydt_F, dqydt] = MakeAtlasField("qy", mesh.cells().size());

  // CFL per cell
  auto [cfl_F, cfl] = MakeAtlasField("CFL", mesh.cells().size());

  // upwinded edge values for fluid height, discharge
  auto [hU_F, hU] = MakeAtlasField("h", mesh.cells().size());
  auto [qUx_F, qUx] = MakeAtlasField("qx", mesh.cells().size());
  auto [qUy_F, qUy] = MakeAtlasField("qy", mesh.cells().size());

  // Geometrical factors on edges
  auto [lambda_F, lambda] = MakeAtlasField("lambda", mesh.edges().size()); // normal velocity
  auto [L_F, L] = MakeAtlasField("L", mesh.edges().size());                // edge length
  auto [nx_F, nx] = MakeAtlasField("nx", mesh.edges().size());             // normals
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
    h(cellIdx, level) = exp(-5 * v * v) + refHeight;
    // h(cellIdx, level) = refHeight;
    // h(cellIdx, level) = sin(xm) * sin(ym) + refHeight;
  }

  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    qx(cellIdx, level) = 0.;
    qy(cellIdx, level) = 0.;
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
  dumpMesh4Triplot(mesh, "init", h, wrapper);

  double t = 0.;
  double dt = 0;
  double t_final = 16.;
  int step = 0;

  // writing this intentionally close to generated code
  while(t < t_final) {

    // make some splashes
    if(step > 0 && step % 1000 == 0) {
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
        xm -= 0;
        ym -= 0;
        double v = sqrt(xm * xm + ym * ym);
        h(cellIdx, level) += exp(-5 * v * v);
      }
    }

    // convert cell centered discharge to velocity and lerp to edges
    {
      const auto& conn = mesh.edges().cell_connectivity();
      for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
        double lhs = 0.;
        double weights[2] = {1 - alpha(edgeIdx, level),
                             alpha(edgeIdx, level)}; // currently not supported in dawn
        for(int nbhIdx = 0; nbhIdx < conn.cols(edgeIdx); nbhIdx++) {
          int cellIdx = conn(edgeIdx, nbhIdx);
          if(cellIdx == conn.missing_value()) {
            assert(weights[nbhIdx] == 0.);
            continue;
          }
          lhs += qx(cellIdx, level) / h(cellIdx, level) * weights[nbhIdx];
        }
        Ux(edgeIdx, level) = lhs;
      }
    }
    {
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
          lhs += qy(cellIdx, level) / h(cellIdx, level) * weights[nbhIdx];
        }
        Uy(edgeIdx, level) = lhs;
      }
    }
    {
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
          lhs += h(cellIdx, level) * weights[nbhIdx];
        }
        hs(edgeIdx, level) = lhs;
      }
    }

    // dumpEdgeField("hs", mesh, wrapper, hs, level);

    // normal edge velocity
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      lambda(edgeIdx, level) =
          nx(edgeIdx, level) * Ux(edgeIdx, level) + ny(edgeIdx, level) * Uy(edgeIdx, level);
    }

    // upwinding for edge values
    //  this pattern is currently unsupported
    {
      const auto& conn = mesh.edges().cell_connectivity();
      for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
        int lo = conn(edgeIdx, 0);
        int hi = conn(edgeIdx, 1);
        hU(edgeIdx, level) = (lambda(edgeIdx, level) < 0) ? h(hi, level) : h(lo, level);
        qUx(edgeIdx, level) = (lambda(edgeIdx, level) < 0) ? qx(hi, level) : qx(lo, level);
        qUy(edgeIdx, level) = (lambda(edgeIdx, level) < 0) ? qy(hi, level) : qy(lo, level);
      }
    }

    // update edge fluxes
    {
      const auto& conn = mesh.edges().cell_connectivity();
      for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
        int cLo = conn(edgeIdx, 0);
        int cHi = conn(edgeIdx, 1);
        bool innerCell = cLo != conn.missing_value() && cHi != conn.missing_value();
        Q(edgeIdx, level) = lambda(edgeIdx, level) * (hU(edgeIdx, level)) * L(edgeIdx, level);
        if(use_corrector && innerCell) {
          double hj = h(cHi, level);
          double hi = h(cLo, level);
          double deltaij = hi - hj;
          Q(edgeIdx, level) -= DampingCoeff * 0.5 * deltaij *
                               sqrt(fabs(Grav) * hU(edgeIdx, level)) * L(edgeIdx, level);
        }
      }
    }
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      Fx(edgeIdx, level) = lambda(edgeIdx, level) * qUx(edgeIdx, level) * L(edgeIdx, level);
    }
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      Fy(edgeIdx, level) = lambda(edgeIdx, level) * qUy(edgeIdx, level) * L(edgeIdx, level);
    }

    // boundary conditions (zero flux)
    // currently not supported in dawn
    for(auto it : boundaryEdges) {
      Q(it, level) = 0;
      Fx(it, level) = 0;
      Fy(it, level) = 0;
      // hs(it, level) = refHeight;
    }

    // dumpEdgeField("L", mesh, wrapper, L, level);
    // return 0;

    // evolve cell values
    {
      const auto& conn = mesh.cells().edge_connectivity();
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double lhs = 0.;
        for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
          int edgeIdx = conn(cellIdx, nbhIdx);
          lhs += Q(edgeIdx, level) * edge_orientation_cell(cellIdx, nbhIdx, level);
        }
        dhdt(cellIdx, level) = lhs;
      }
    }
    {
      const auto& conn = mesh.cells().edge_connectivity();
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double lhs = 0.;
        for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
          int edgeIdx = conn(cellIdx, nbhIdx);
          lhs += Fx(edgeIdx, level) * edge_orientation_cell(cellIdx, nbhIdx, level);
        }
        dqxdt(cellIdx, level) = lhs / A(cellIdx, level);
        if(use_friction) {
          double lenq = sqrt(qx(cellIdx, level) * qx(cellIdx, level) +
                             qy(cellIdx, level) * qy(cellIdx, level));
          dqxdt(cellIdx, level) -= Grav * ManningCoeff * ManningCoeff /
                                   pow(h(cellIdx, level), 10. / 3.) * lenq * qx(cellIdx, level);
        }
      }
    }
    {
      const auto& conn = mesh.cells().edge_connectivity();
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double lhs = 0.;
        for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
          int edgeIdx = conn(cellIdx, nbhIdx);
          lhs += Fy(edgeIdx, level) * edge_orientation_cell(cellIdx, nbhIdx, level);
        }
        dqydt(cellIdx, level) = lhs / A(cellIdx, level);
        if(use_friction) {
          double lenq = sqrt(qx(cellIdx, level) * qx(cellIdx, level) +
                             qy(cellIdx, level) * qy(cellIdx, level));
          dqydt(cellIdx, level) -= Grav * ManningCoeff * ManningCoeff /
                                   pow(h(cellIdx, level), 10. / 3.) * lenq * qy(cellIdx, level);
        }
      }
    }
    {
      const auto& conn = mesh.cells().edge_connectivity();
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double lhs = 0.;
        for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
          int edgeIdx = conn(cellIdx, nbhIdx);
          lhs -= hs(edgeIdx, level) * nx(edgeIdx, level) *
                 edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
        }
        Sx(cellIdx, level) = lhs / A(cellIdx, level);
      }
    }
    {
      const auto& conn = mesh.cells().edge_connectivity();
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double lhs = 0.;
        for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
          int edgeIdx = conn(cellIdx, nbhIdx);
          lhs -= hs(edgeIdx, level) * ny(edgeIdx, level) *
                 edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
        }
        Sy(cellIdx, level) = lhs / A(cellIdx, level);
      }
    }
    for(auto it : boundaryCells) {
      Sx(it, level) = 0.;
      Sy(it, level) = 0.;
    }
    // dumpEdgeField("hs", mesh, wrapper, hs, level);
    // dumpCellField("Sx", mesh, wrapper, Sx, level);
    // dumpCellField("Sy", mesh, wrapper, Sy, level);

    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      dhdt(cellIdx, level) = dhdt(cellIdx, level) / A(cellIdx, level) * dt;
      dqxdt(cellIdx, level) =
          (dqxdt(cellIdx, level) - Grav * (h(cellIdx, level)) * Sx(cellIdx, level)) * dt;
      dqydt(cellIdx, level) =
          (dqydt(cellIdx, level) - Grav * (h(cellIdx, level)) * Sy(cellIdx, level)) * dt;
    }
    for(auto it : boundaryCells) {
      dhdt(it, level) = 0.;
      dqxdt(it, level) = 0.;
      dqydt(it, level) = 0.;
    }
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      h(cellIdx, level) = h(cellIdx, level) + dhdt(cellIdx, level);
      qx(cellIdx, level) = qx(cellIdx, level) - dqxdt(cellIdx, level);
      qy(cellIdx, level) = qy(cellIdx, level) - dqydt(cellIdx, level);
    }

    // dumpCellField("h", mesh, wrapper, h, level);
    // dumpCellField("dhdt", mesh, wrapper, dhdt, level);
    // dumpCellField("dqxdt", mesh, wrapper, dqxdt, level);
    // dumpCellField("dqydt", mesh, wrapper, dqydt, level);
    // if(step > 1) {
    //   return 0;
    // }

    // adapt CLF
    // this would probably be in the driver code anyway
    {
      const auto& conn = mesh.cells().edge_connectivity();
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double l0 = L(conn(cellIdx, 0), level);
        double l1 = L(conn(cellIdx, 1), level);
        double l2 = L(conn(cellIdx, 2), level);
        double hi = h(cellIdx, level);
        double Ux = qx(cellIdx, level) / hi;
        double Uy = qy(cellIdx, level) / hi;
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
      // sprintf(buf, "out/step_%04d.txt", step);

      sprintf(buf, "out/stepH_%04d.txt", step);
      dumpCellField(buf, mesh, wrapper, h, level);
      // dumpCellFieldOnNodes(buf, mesh, wrapper, h, level);
    }
    std::cout << "time " << t << " timestep " << step++ << " dt " << dt << "\n";
  }

  dumpMesh4Triplot(mesh, "final", h, wrapper);
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