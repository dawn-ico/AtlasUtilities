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

// Tests grad / div / curl on the surface of a sphere using 2 dimensional FV primitives

#include <cmath>
#include <cstdio>
#include <fenv.h>
#include <glm/fwd.hpp>
#include <mesh/actions/BuildCellCentres.h>
#include <optional>
#include <set>
#include <vector>

#include <glm/glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>

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
#include "../utils/AtlasFromNetcdf.h"

template <typename T>
static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

const double RADIUS_SPHERE = 1;
// const double RADIUS_SPHERE = 6378e3;

void dumpMesh4Triplot(const atlas::Mesh& mesh, const std::string prefix,
                      const atlasInterface::Field<double>& field,
                      const std::vector<glm::dvec3>& xyz);

void dumpMesh4Triplot(const atlas::Mesh& mesh, const std::string prefix,
                      const std::vector<glm::dvec3>& xyz);

void dumpNodeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level);
void dumpCellField(const std::string& fname, const atlas::Mesh& mesh,
                   atlasInterface::Field<double>& field, int level);
void dumpCellFieldOnNodes(const std::string& fname, const atlas::Mesh& mesh,
                          AtlasToCartesian wrapper, atlasInterface::Field<double>& field,
                          int level);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh,
                   const std::vector<glm::dvec3>& xyz, atlasInterface::Field<double>& field,
                   int level);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level, std::vector<int> edgeList,
                   std::optional<Orientation> color = std::nullopt);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field_x, atlasInterface::Field<double>& field_y,
                   int level, std::optional<Orientation> color = std::nullopt);

//===-----------------------------------------------------------------------------

std::tuple<glm::dvec3, glm::dvec3> cartEdge(const atlas::Mesh& mesh,
                                            const std::vector<glm::dvec3>& xyz, size_t edgeIdx) {
  const auto& conn = mesh.edges().node_connectivity();
  return {xyz[conn(edgeIdx, 0)], xyz[conn(edgeIdx, 1)]};
}

double edgeLength(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz, size_t edgeIdx) {
  auto [p1, p2] = cartEdge(mesh, xyz, edgeIdx);
  return glm::length(p1 - p2);
}

glm::dvec3 edgeMidpoint(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz,
                        size_t edgeIdx) {
  auto [p1, p2] = cartEdge(mesh, xyz, edgeIdx);
  return 0.5 * (p1 + p2);
}

glm::dvec3 cellCircumcenter(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz,
                            int cellIdx) {
  const auto& cellNodeConnectivity = mesh.cells().node_connectivity();
  const int missingVal = cellNodeConnectivity.missing_value();

  // only valid for tringular cells with all node neighbors set
  int numNbh = cellNodeConnectivity.cols(cellIdx);
  assert(numNbh == 3);
  for(int nbh = 0; nbh < numNbh; nbh++) {
    int nbhIdx = cellNodeConnectivity(cellIdx, nbh);
    assert(nbhIdx != missingVal);
  }

  glm::dvec3 a = xyz[(cellNodeConnectivity(cellIdx, 0))];
  glm::dvec3 b = xyz[(cellNodeConnectivity(cellIdx, 1))];
  glm::dvec3 c = xyz[(cellNodeConnectivity(cellIdx, 2))];

  // https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
  glm::dvec3 ac = c - a;
  glm::dvec3 ab = b - a;
  glm::dvec3 abXac = glm::cross(ab, ac);

  // this is the vector from a TO the circumsphere center
  glm::dvec3 toCircumsphereCenter =
      (glm::cross(abXac, ab) * glm::length2(ac) + glm::cross(ac, abXac) * glm::length2(ab)) /
      (2.f * glm::length2(abXac));

  return a + toCircumsphereCenter;
}

double dualEdgeLength(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz, size_t edgeIdx) {
  const auto& conn = mesh.edges().cell_connectivity();
  auto p1 = cellCircumcenter(mesh, xyz, conn(edgeIdx, 0));
  auto p2 = cellCircumcenter(mesh, xyz, conn(edgeIdx, 1));
  return glm::length(p1 - p2);
}

glm::dvec3 primalNormal(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz,
                        size_t edgeIdx) {
  const auto& conn = mesh.edges().cell_connectivity();
  glm::dvec3 c0 = cellCircumcenter(mesh, xyz, conn(edgeIdx, 0));
  glm::dvec3 c1 = cellCircumcenter(mesh, xyz, conn(edgeIdx, 1));
  return glm::normalize(c1 - c0);
}

double distanceToCircumcenter(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz,
                              size_t cellIdx, size_t edgeIdx) {
  glm::dvec3 x0 = cellCircumcenter(mesh, xyz, cellIdx);
  auto [x1, x2] = cartEdge(mesh, xyz, edgeIdx);
  // https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  return glm::length(glm::cross(x0 - x1, x0 - x2)) / glm::length(x2 - x1);
}

double triangleArea(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) {
  return 0.5 * glm::length(glm::cross(v1 - v0, v2 - v0));
}

double cellArea(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz, size_t cellIdx) {
  const auto& cellNodeConnectivity = mesh.cells().node_connectivity();
  const int missingVal = cellNodeConnectivity.missing_value();

  // only valid for triangular cells with all node neighbors set
  int numNbh = cellNodeConnectivity.cols(cellIdx);
  assert(numNbh == 3);
  for(int nbh = 0; nbh < numNbh; nbh++) {
    int nbhIdx = cellNodeConnectivity(cellIdx, nbh);
    assert(nbhIdx != missingVal);
  }

  glm::dvec3 v0 = xyz[(cellNodeConnectivity(cellIdx, 0))];
  glm::dvec3 v1 = xyz[(cellNodeConnectivity(cellIdx, 1))];
  glm::dvec3 v2 = xyz[(cellNodeConnectivity(cellIdx, 2))];

  return triangleArea(v0, v1, v2);
}

double dualCellArea(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz, int nodeIdx) {
  const atlas::mesh::Nodes::Connectivity& nodeEdgeConnectivity = mesh.nodes().edge_connectivity();
  const atlas::mesh::HybridElements::Connectivity& edgeCellConnectivity =
      mesh.edges().cell_connectivity();
  double totalArea = 0.;
  const int missingValEdge = nodeEdgeConnectivity.missing_value();
  const int missingValCell = edgeCellConnectivity.missing_value();

  auto center = xyz[nodeIdx];

  int numNbh = nodeEdgeConnectivity.cols(nodeIdx);
  assert(numNbh == 6 || numNbh == 5);

  for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
    int edgeIdx = nodeEdgeConnectivity(nodeIdx, nbhIdx);
    assert(edgeIdx != missingValEdge);

    int numNbhCells = edgeCellConnectivity.cols(edgeIdx);
    assert(numNbhCells == 2);

    int cellIdxLo = edgeCellConnectivity(edgeIdx, 0);
    int cellIdxHi = edgeCellConnectivity(edgeIdx, 1);

    assert(cellIdxLo != missingValCell && cellIdxHi != missingValCell);

    auto pLo = cellCircumcenter(mesh, xyz, cellIdxLo);
    auto pHi = cellCircumcenter(mesh, xyz, cellIdxHi);

    totalArea += triangleArea(center, pLo, pHi);
  }
  return totalArea;
}

glm::dvec2 myprojLTPC_ENU(const glm::dvec3& pCart, const glm::dvec3& pRefCart) {
  double lambda = atan2(pRefCart.y, pRefCart.x);
  double phi = asin((pRefCart.z) / (RADIUS_SPHERE));

  // glm uses column major ordering
  glm::dmat3x3 R(-sin(lambda), -sin(phi) * cos(lambda), cos(phi) * cos(lambda), cos(lambda),
                 -sin(phi) * sin(lambda), cos(phi) * sin(lambda), 0., cos(phi), sin(phi));

  glm::dvec3 ret = R * (pRefCart - pCart); // this is in east, north, up
  // assert(fabs(ret.z) < 1e-2);              // expecting this to be small
  return {ret.x, ret.y}; // do not need up
}

std::tuple<double, double, double> MeasureErrors(int numEl,
                                                 const atlasInterface::Field<double>& ref,
                                                 const atlasInterface::Field<double>& sol,
                                                 int level);

//===-----------------------------------------------------------------------------

int main(int argc, char const* argv[]) {
  // enable floating point exception
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  if(argc != 2) {
    std::cout << "intended use is\n" << argv[0] << " <mesh>.nc" << std::endl;
    return -1;
  }

  int k_size = 1;
  const int level = 0;

  // dump a whole bunch of debug output (meant to be visualized using Octave, but gnuplot and the
  // like will certainly work too)
  const bool dbg_out = true;

  atlas::Mesh mesh = AtlasMeshFromNetCDFMinimal(argv[1]).value();
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

  // netcdf mesh has lon/lat set in degrees, we also want cartesian coordinates here
  std::vector<glm::dvec3> xyz(mesh.nodes().size());
  {
    auto latToRad = [](double rad) { return rad / 90. * (0.5 * M_PI); };
    auto lonToRad = [](double rad) { return rad / 180. * M_PI; };
    auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
    auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
    for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      double lon = lonToRad(lonlat(nodeIdx, atlas::LON));
      double lat = latToRad(lonlat(nodeIdx, atlas::LAT));

      double x = RADIUS_SPHERE * cos(lat) * cos(lon);
      double y = RADIUS_SPHERE * cos(lat) * sin(lon);
      double z = RADIUS_SPHERE * sin(lat);

      double latRet = asin(z / RADIUS_SPHERE);
      double lonRet = atan2(y, x);

      assert(fabs(latRet - lat) < 1e-6);
      assert(fabs(lonRet - lon) < 1e-6);

      xyz[nodeIdx] = {x, y, z};
      lonlat(nodeIdx, atlas::LON) = lon;
      lonlat(nodeIdx, atlas::LAT) = lat;
    }
  }

  if(dbg_out) {
    dumpMesh4Triplot(mesh, "shallow", xyz);
  }

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

  auto [u_F, u] = MakeAtlasField("u", mesh.edges().size());
  auto [v_F, v] = MakeAtlasField("v", mesh.edges().size());
  auto [uvn_F, uvn] = MakeAtlasField("v", mesh.edges().size());

  auto [div_F, div] = MakeAtlasField("div_uvn", mesh.cells().size());
  auto [grad_u_x_F, grad_u_x] = MakeAtlasField("grad_u_x", mesh.cells().size());
  auto [grad_u_y_F, grad_u_y] = MakeAtlasField("grad_u_y", mesh.cells().size());
  auto [grad_v_x_F, grad_v_x] = MakeAtlasField("grad_v_x", mesh.cells().size());
  auto [grad_v_y_F, grad_v_y] = MakeAtlasField("grad_v_y", mesh.cells().size());
  auto [curl_F, curl] = MakeAtlasField("curl", mesh.nodes().size());

  auto [divSol_F, divSol] = MakeAtlasField("divSol_uvn", mesh.cells().size());
  auto [gradSol_u_x_F, gradSol_u_x] = MakeAtlasField("gradSol_u_x", mesh.cells().size());
  auto [gradSol_u_y_F, gradSol_u_y] = MakeAtlasField("gradSol_u_y", mesh.cells().size());
  auto [gradSol_v_x_F, gradSol_v_x] = MakeAtlasField("gradSol_v_x", mesh.cells().size());
  auto [gradSol_v_y_F, gradSol_v_y] = MakeAtlasField("gradSol_v_y", mesh.cells().size());
  auto [curlSol_F, curlSol] = MakeAtlasField("curlSol", mesh.nodes().size());

  // Geometrical factors on edges
  auto [lambda_F, lambda] = MakeAtlasField("lambda", mesh.edges().size()); // normal velocity
  auto [L_F, L] = MakeAtlasField("L", mesh.edges().size());                // edge length
  auto [dualL_F, dualL] = MakeAtlasField("dualL", mesh.edges().size());    // edge length
  auto [nx_F, nx] = MakeAtlasField("nx", mesh.edges().size());             // normals
  auto [ny_F, ny] = MakeAtlasField("ny", mesh.edges().size());
  auto [nz_F, nz] = MakeAtlasField("nz", mesh.edges().size());

  // Geometrical factors on cells
  auto [A_F, A] = MakeAtlasField("A", mesh.cells().size());
  auto [edge_orientation_cell_F, edge_orientation_cell] =
      MakeAtlasSparseField("edge_orientation_cell", mesh.cells().size(), edgesPerCell);
  auto [nxLCce_F, nxLCce] = MakeAtlasSparseField("nxLCce", mesh.cells().size(), edgesPerCell);
  auto [nyLCce_F, nyLCce] = MakeAtlasSparseField("nyLCce", mesh.cells().size(), edgesPerCell);

  // Geometrical factors on nodes
  auto [dualA_F, dualA] = MakeAtlasField("dualA", mesh.nodes().size());
  auto [edge_orientation_vertex_F, edge_orientation_vertex] =
      MakeAtlasSparseField("edge_orientation_vertex", mesh.nodes().size(), edgesPerVertex);

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on cells
  //===------------------------------------------------------------------------------------------===//
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    A(cellIdx, level) = cellArea(mesh, xyz, cellIdx);
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on nodes
  //===------------------------------------------------------------------------------------------===//
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    dualA(nodeIdx, level) = dualCellArea(mesh, xyz, nodeIdx);
  }

  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    const auto& nodeEdgeConnectivity = mesh.nodes().edge_connectivity();
    const auto& edgeNodeConnectivity = mesh.edges().node_connectivity();

    const int missingVal = nodeEdgeConnectivity.missing_value();
    int numNbh = nodeEdgeConnectivity.cols(nodeIdx);

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = nodeEdgeConnectivity(nodeIdx, nbhIdx);

      int n0 = edgeNodeConnectivity(edgeIdx, 0);
      int n1 = edgeNodeConnectivity(edgeIdx, 1);

      int centerIdx = (n0 == nodeIdx) ? n0 : n1;
      int farIdx = (n0 == nodeIdx) ? n1 : n0;

      glm::dvec3 edge = xyz[farIdx] - xyz[centerIdx];
      glm::dvec3 norm = primalNormal(mesh, xyz, edgeIdx);
      glm::dvec2 normENU = myprojLTPC_ENU(xyz[nodeIdx] + norm, xyz[nodeIdx]);
      glm::dvec2 edgeENU = myprojLTPC_ENU(xyz[nodeIdx] + edge, xyz[nodeIdx]);
      glm::dvec2 dualNormENU{normENU.y, -normENU.x};

      edge_orientation_vertex(nodeIdx, nbhIdx, level) = sgn(glm::dot(edgeENU, dualNormENU));
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    L(edgeIdx, level) = edgeLength(mesh, xyz, edgeIdx);
    dualL(edgeIdx, level) = dualEdgeLength(mesh, xyz, edgeIdx);
  }

  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    const atlas::mesh::HybridElements::Connectivity& cellEdgeConnectivity =
        mesh.cells().edge_connectivity();

    glm::dvec3 c = cellCircumcenter(mesh, xyz, cellIdx);

    int numNbh = cellEdgeConnectivity.cols(cellIdx);
    assert(numNbh == edgesPerCell);

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = cellEdgeConnectivity(cellIdx, nbhIdx);
      glm::dvec3 toOutside = edgeMidpoint(mesh, xyz, edgeIdx) - c;
      glm::dvec3 n = primalNormal(mesh, xyz, edgeIdx);
      glm::dvec2 nENU = myprojLTPC_ENU(c + n, c);
      glm::dvec2 toOutsideENU = myprojLTPC_ENU(c + toOutside, c);
      nxLCce(cellIdx, nbhIdx, level) = nENU.x;
      nyLCce(cellIdx, nbhIdx, level) = nENU.y;
      edge_orientation_cell(cellIdx, nbhIdx, level) = -sgn(glm::dot(toOutsideENU, nENU));
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize test function & analytical solutions
  //===------------------------------------------------------------------------------------------===//
  auto sphericalHarmonic = [](double lambda, double phi) -> std::tuple<double, double> {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return {c1 * cos(2 * lambda) * cos(phi) * cos(phi) * sin(phi),
            c2 * cos(lambda) * cos(phi) * sin(phi)};
  };
  // NOTE: del operator on curvilinar derivatives is different from del operator in cartesian
  // coordinates, see
  // https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates#Del_formula However,
  // ALSO NOTE that phi in geoscience is phase shifted by PI/2!!!!1
  auto analyticalGradientU = [](double lambda, double phi) -> std::tuple<double, double> {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    return {(1 / sin(phi + M_PI_2)) * c1 * -2 * sin(2 * lambda) * cos(phi) * cos(phi) * sin(phi),
            c1 * cos(2 * lambda) * cos(phi) * (cos(phi) * cos(phi) - 2 * sin(phi) * sin(phi))};
  };
  auto analyticalGradientV = [](double lambda, double phi) -> std::tuple<double, double> {
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return {(1 / sin(phi + M_PI_2)) * c2 * sin(lambda) * sin(phi) * (-cos(phi)),
            c2 * cos(lambda) * (cos(phi) * cos(phi) - sin(phi) * sin(phi))};
  };
  auto analyticalDivergence = [](double lambda, double phi) {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));

    double u = c1 * cos(2 * lambda) * cos(phi) * cos(phi) * sin(phi);
    double v = c2 * cos(lambda) * cos(phi) * sin(phi);

    double u_lambda = c1 * -2 * sin(2 * lambda) * cos(phi) * cos(phi) * sin(phi);
    double v_phi = c2 * cos(lambda) * (cos(phi) * cos(phi) - sin(phi) * sin(phi));

    return 1 / sin(phi + M_PI_2) * u_lambda + v_phi + cos(phi + M_PI_2) / sin(phi + M_PI_2) * v;
  };
  auto analyticalCurl = [](double lambda, double phi) {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));

    double u = c1 * cos(2 * lambda) * cos(phi) * cos(phi) * sin(phi);
    double v = c2 * cos(lambda) * cos(phi) * sin(phi);

    double dudlambda =
        c1 * cos(2 * lambda) * cos(phi) * (cos(phi) * cos(phi) - 2 * sin(phi) * sin(phi));
    double dvdphi = -c2 * cos(phi) * sin(lambda) * sin(phi);
    return dudlambda + cos(phi + M_PI_2) / sin(phi + M_PI_2) * u - dvdphi / sin(phi + M_PI_2);
  };

  // test fun
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    glm::dvec3 m = edgeMidpoint(mesh, xyz, edgeIdx);
    glm::dvec3 n = primalNormal(mesh, xyz, edgeIdx);
    glm::dvec2 nENU = myprojLTPC_ENU(m + n, m);
    double lambdaCyl = atan2(m.y, m.x);
    double phiCyl = asin((m.z) / (RADIUS_SPHERE));
    auto [uIdx, vIdx] = sphericalHarmonic(lambdaCyl, phiCyl);
    u(edgeIdx, level) = uIdx;
    v(edgeIdx, level) = vIdx;
    uvn(edgeIdx, level) = uIdx * nENU.x + vIdx * nENU.y;
  }
  // solutions
  for(int cellIdx = 0; cellIdx < mesh->cells().size(); cellIdx++) {
    glm::dvec3 m = cellCircumcenter(mesh, xyz, cellIdx);
    double lambdaCyl = atan2(m.y, m.x);
    double phiCyl = asin((m.z) / (RADIUS_SPHERE));
    divSol(cellIdx, level) = analyticalDivergence(lambdaCyl, phiCyl);
    auto [gradUx, gradUy] = analyticalGradientU(lambdaCyl, phiCyl);
    auto [gradVx, gradVy] = analyticalGradientV(lambdaCyl, phiCyl);
    gradSol_u_x(cellIdx, level) = gradUx;
    gradSol_u_y(cellIdx, level) = gradUy;
    gradSol_v_x(cellIdx, level) = gradVx;
    gradSol_v_y(cellIdx, level) = gradVy;
  }
  for(int nodeIdx = 0; nodeIdx < mesh->nodes().size(); nodeIdx++) {
    auto m = xyz[nodeIdx];
    double phiCyl = asin((m.z) / (RADIUS_SPHERE));
    double lambdaCyl = atan2(m.y, m.x);
    curlSol(nodeIdx, level) = analyticalCurl(lambdaCyl, phiCyl);
  }

  //===------------------------------------------------------------------------------------------===//
  // gradients
  //===------------------------------------------------------------------------------------------===//
  {
    const auto& conn = mesh.cells().edge_connectivity();
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      double lhs_x = 0.;
      double lhs_y = 0.;
      for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
        int edgeIdx = conn(cellIdx, nbhIdx);
        lhs_x += u(edgeIdx, level) * nxLCce(cellIdx, nbhIdx, level) *
                 edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
        lhs_y += u(edgeIdx, level) * nyLCce(cellIdx, nbhIdx, level) *
                 edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
      }
      grad_u_x(cellIdx, level) = lhs_x / A(cellIdx, level);
      grad_u_y(cellIdx, level) = lhs_y / A(cellIdx, level);
    }
  }
  {
    const auto& conn = mesh.cells().edge_connectivity();
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      double lhs_x = 0.;
      double lhs_y = 0.;
      for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
        int edgeIdx = conn(cellIdx, nbhIdx);
        lhs_x += v(edgeIdx, level) * nxLCce(cellIdx, nbhIdx, level) *
                 edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
        lhs_y += v(edgeIdx, level) * nyLCce(cellIdx, nbhIdx, level) *
                 edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
      }
      grad_v_x(cellIdx, level) = lhs_x / A(cellIdx, level);
      grad_v_y(cellIdx, level) = lhs_y / A(cellIdx, level);
    }
  }
  //===------------------------------------------------------------------------------------------===//
  // divergence
  //===------------------------------------------------------------------------------------------===//
  {
    const auto& conn = mesh.cells().edge_connectivity();
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      double lhs = 0.;
      for(int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
        int edgeIdx = conn(cellIdx, nbhIdx);
        lhs +=
            uvn(edgeIdx, level) * edge_orientation_cell(cellIdx, nbhIdx, level) * L(edgeIdx, level);
      }
      div(cellIdx, level) = lhs / A(cellIdx, level);
    }
  }
  //===------------------------------------------------------------------------------------------===//
  // curl
  //===------------------------------------------------------------------------------------------===//
  {
    const auto& conn = mesh->nodes().edge_connectivity();
    for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      double lhs = 0.;
      for(int nbhIdx = 0; nbhIdx < conn.cols(nodeIdx); nbhIdx++) {
        int edgeIdx = conn(nodeIdx, nbhIdx);
        lhs += uvn(edgeIdx, level) * edge_orientation_vertex(nodeIdx, nbhIdx, level) *
               dualL(edgeIdx, level);
      }
      curl(nodeIdx, level) = lhs / dualA(nodeIdx, level);
    }
  }

  if(dbg_out) {
    {
      FILE* fp = fopen("gradU.txt", "w+");
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        fprintf(fp, "%f %f %f %f\n", grad_u_x(cellIdx, level), grad_u_y(cellIdx, level),
                gradSol_u_x(cellIdx, level), gradSol_u_y(cellIdx, level));
      }
      fclose(fp);
    }
    {
      FILE* fp = fopen("gradV.txt", "w+");
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        fprintf(fp, "%f %f %f %f\n", grad_v_x(cellIdx, level), grad_v_y(cellIdx, level),
                gradSol_v_x(cellIdx, level), gradSol_v_y(cellIdx, level));
      }
      fclose(fp);
    }
    {
      FILE* fp = fopen("div.txt", "w+");
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        fprintf(fp, "%f %f\n", div(cellIdx, level), divSol(cellIdx, level));
      }
      fclose(fp);
    }
    {
      const auto& conn = mesh->cells().node_connectivity();
      FILE* fp = fopen("curl.txt", "w+");
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double curlIdx = 0.;
        double curlSolIdx = 0.;
        for(int nbhIdx = 0; nbhIdx < 3; nbhIdx++) {
          curlIdx += curl(conn(cellIdx, nbhIdx), level);
          curlSolIdx += curlSol(conn(cellIdx, nbhIdx), level);
        }
        curlIdx /= 3.;
        curlSolIdx /= 3.;
        fprintf(fp, "%f %f\n", curlIdx, curlSolIdx);
      }
      fclose(fp);
    }
    {
      const auto& conn = mesh->cells().edge_connectivity();
      const auto& connNodes = mesh->cells().node_connectivity();
      FILE* fp = fopen("in.csv", "w+");
      fprintf(fp, "lamba, phi, f, f_lambda, f_phi, f_lambda_ref, f_phi_ref, div, div_ref, curl, "
                  "curl_ref\n");
      for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        double uIdx = 0.;
        double vIdx = 0.;
        double curlIdx = 0.;
        double curlSolIdx = 0.;
        for(int nbhIdx = 0; nbhIdx < 3; nbhIdx++) {
          uIdx += u(conn(cellIdx, nbhIdx), level);
          vIdx += v(conn(cellIdx, nbhIdx), level);
        }
        for(int nbhIdx = 0; nbhIdx < 3; nbhIdx++) {
          curlIdx += curl(connNodes(cellIdx, nbhIdx), level);
          curlSolIdx += curlSol(connNodes(cellIdx, nbhIdx), level);
        }
        uIdx /= 3.;
        vIdx /= 3.;
        curlIdx /= 3.;
        curlSolIdx /= 3.;
        auto m = cellCircumcenter(mesh, xyz, cellIdx);
        double lambdaCyl = atan2(m.y, m.x);
        double phiCyl = asin((m.z) / (RADIUS_SPHERE));
        fprintf(fp, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", lambdaCyl, phiCyl, uIdx,
                grad_u_x(cellIdx, level), grad_u_y(cellIdx, level), gradSol_u_x(cellIdx, level),
                gradSol_u_y(cellIdx, level), div(cellIdx, level), divSol(cellIdx, level), curlIdx,
                curlSolIdx);
      }
      fclose(fp);
    }
  }

  {
    auto [Linf, L1, L2] = MeasureErrors(mesh->cells().size(), gradSol_u_x, grad_u_x, level);
    printf("gradUx: %e %e %e\n", Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureErrors(mesh->cells().size(), gradSol_u_y, grad_u_y, level);
    printf("gradUy: %e %e %e\n", Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureErrors(mesh->cells().size(), gradSol_v_x, grad_v_x, level);
    printf("gradVx: %e %e %e\n", Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureErrors(mesh->cells().size(), gradSol_v_y, grad_v_y, level);
    printf("gradVy: %e %e %e\n", Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureErrors(mesh->cells().size(), divSol, div, level);
    printf("div: %e %e %e\n", Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureErrors(mesh->nodes().size(), curlSol, curl, level);
    printf("curl: %e %e %e\n", Linf, L1, L2);
  }
}

std::tuple<double, double, double> MeasureErrors(int numEl,
                                                 const atlasInterface::Field<double>& ref,
                                                 const atlasInterface::Field<double>& sol,
                                                 int level) {
  double Linf = 0.;
  double L1 = 0.;
  double L2 = 0.;
  for(int idx = 0; idx < numEl; idx++) {
    double dif = ref(idx, level) - sol(idx, level);
    Linf = fmax(fabs(dif), Linf);
    L1 += fabs(dif);
    L2 += dif * dif;
  }
  L1 /= numEl;
  L2 = sqrt(L2) / sqrt(numEl);
  return {Linf, L1, L2};
}

void dumpMesh4Triplot(const atlas::Mesh& mesh, const std::string prefix,
                      const atlasInterface::Field<double>& field,
                      const std::vector<glm::dvec3>& xyz) {
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
    for(const auto& it : xyz) {
      fprintf(fp, "%f %f %f\n", it.x, it.y, it.z);
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

void dumpMesh4Triplot(const atlas::Mesh& mesh, const std::string prefix,
                      const std::vector<glm::dvec3>& xyz) {
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
    for(const auto& it : xyz) {
      fprintf(fp, "%f %f %f\n", it.x, it.y, it.z);
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

void dumpCellField(const std::string& fname, const atlas::Mesh& mesh,
                   atlasInterface::Field<double>& field, int level) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    fprintf(fp, "%f\n", field(cellIdx, level));
  }
  fclose(fp);
}

void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh,
                   const std::vector<glm::dvec3>& xyz, atlasInterface::Field<double>& field,
                   int level) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    glm::dvec3 pos = edgeMidpoint(mesh, xyz, edgeIdx);
    fprintf(fp, "%f %f %f %f\n", pos.x, pos.y, pos.z, field(edgeIdx, level));
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
