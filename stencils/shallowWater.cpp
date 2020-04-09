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
//    current state: normals at the poles are weird (not orthogonal to mesh projected to tangent
//    plane). cell midpoints with such normals are printed to console. debug from there!
//
//===------------------------------------------------------------------------------------------===//

// Shallow water equation solver as described in "A simple and efficient unstructured finite volume
// scheme for solving the shallow water equations in overland flow applications" by Cea and Bladé
// Follows notation in the paper as closely as possilbe

#include <cmath>
#include <cstdio>
#include <fenv.h>
#include <glm/fwd.hpp>
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

glm::dvec2 edgeMidpoint2(const atlas::Mesh& mesh, size_t edgeIdx) {
  auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
  const auto& conn = mesh.edges().node_connectivity();
  glm::dvec2 p1{lonlat(conn(edgeIdx, 0), atlas::LON), lonlat(conn(edgeIdx, 0), atlas::LAT)};
  glm::dvec2 p2{lonlat(conn(edgeIdx, 1), atlas::LON), lonlat(conn(edgeIdx, 1), atlas::LAT)};
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

glm::dvec2 cellCircumcenter2(const atlas::Mesh& mesh, int cellIdx) {
  const auto& cellNodeConnectivity = mesh.cells().node_connectivity();
  const int missingVal = cellNodeConnectivity.missing_value();

  // only valid for tringular cells with all node neighbors set
  int numNbh = cellNodeConnectivity.cols(cellIdx);
  assert(numNbh == 3);
  for(int nbh = 0; nbh < numNbh; nbh++) {
    int nbhIdx = cellNodeConnectivity(cellIdx, nbh);
    assert(nbhIdx != missingVal);
  }

  auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());

  auto [Ax, Ay] = std::tuple<double, double>{lonlat(cellNodeConnectivity(cellIdx, 0), atlas::LON),
                                             lonlat(cellNodeConnectivity(cellIdx, 0), atlas::LAT)};
  auto [Bx, By] = std::tuple<double, double>{lonlat(cellNodeConnectivity(cellIdx, 1), atlas::LON),
                                             lonlat(cellNodeConnectivity(cellIdx, 1), atlas::LAT)};
  auto [Cx, Cy] = std::tuple<double, double>{lonlat(cellNodeConnectivity(cellIdx, 2), atlas::LON),
                                             lonlat(cellNodeConnectivity(cellIdx, 2), atlas::LAT)};

  double D = 2 * (Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By));
  double Ux = 1. / D *
              ((Ax * Ax + Ay * Ay) * (By - Cy) + (Bx * Bx + By * By) * (Cy - Ay) +
               (Cx * Cx + Cy * Cy) * (Ay - By));
  double Uy = 1. / D *
              ((Ax * Ax + Ay * Ay) * (Cx - Bx) + (Bx * Bx + By * By) * (Ax - Cx) +
               (Cx * Cx + Cy * Cy) * (Bx - Ax));
  return {Ux, Uy};
}

glm::dvec3 primalNormal(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz,
                        size_t edgeIdx) {
  const auto& conn = mesh.edges().cell_connectivity();
  glm::dvec3 c0 = cellCircumcenter(mesh, xyz, conn(edgeIdx, 0));
  glm::dvec3 c1 = cellCircumcenter(mesh, xyz, conn(edgeIdx, 1));
  return glm::normalize(c1 - c0);
}

glm::dvec2 primalNormal2(const atlas::Mesh& mesh, size_t edgeIdx) {
  const auto& conn = mesh.edges().cell_connectivity();
  glm::dvec2 c0 = cellCircumcenter2(mesh, conn(edgeIdx, 0));
  glm::dvec2 c1 = cellCircumcenter2(mesh, conn(edgeIdx, 1));
  return glm::normalize(c1 - c0);
}

double distanceToCircumcenter(const atlas::Mesh& mesh, const std::vector<glm::dvec3>& xyz,
                              size_t cellIdx, size_t edgeIdx) {
  glm::dvec3 x0 = cellCircumcenter(mesh, xyz, cellIdx);
  auto [x1, x2] = cartEdge(mesh, xyz, edgeIdx);
  // https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  return glm::length(glm::cross(x0 - x1, x0 - x2)) / glm::length(x2 - x1);
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

  return 0.5 * glm::length(glm::cross(v1 - v0, v2 - v0));
}

// lifted from the ICON source. supposed to
//   !! Convert cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$
//   !! into zonal @f$p\_gu@f$ and meridional @f$g\_v@f$ velocity components
//
// this is the trafo here
// https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
// makes sense, east = x, north = y, up = dropped
//
// glm::dvec2 cvec2gvec(const glm::dvec2& lonlat, const glm::dvec3& cartvec) {
//   double p_long = lonlat.x;
//   double p_lat = lonlat.y;
//   return {-sin(p_long) * cartvec.x + cos(p_long) * cartvec.y,
//           -sin(p_lat) * cos(p_long) * cartvec.x - sin(p_lat) * sin(p_long) * cartvec.y +
//               cos(p_lat) * cartvec.z};
// }

glm::dvec2 cvec2gvec(glm::dvec2 lonlat, glm::dvec3 cartvec) {
  double p_long = lonlat.x;
  double p_lat = lonlat.y;

  double p_cu = cartvec.x;
  double p_cv = cartvec.y;
  double p_cw = cartvec.z;

  double z_sln = sin(p_long);
  double z_cln = cos(p_long);
  double z_slt = sin(p_lat);
  double z_clt = cos(p_lat);

  double p_gu = z_cln * p_cv - z_sln * p_cu;
  double p_gv = z_cln * p_cu + z_sln * p_cv;
  p_gv = z_slt * p_gv;
  p_gv = z_clt * p_cw - p_gv;
  return {p_gu, p_gv};
}

// seems to be equivalent to above
// glm::dvec2 cc2tv(const glm::dvec2& position, const glm::dvec3& xx) {
//   double z_sinlo = sin(position.x);
//   double z_coslo = cos(position.x);
//   double z_sinlacoslo = sin(position.y) * z_coslo;
//   double z_sinlasinlo = sin(position.y) * z_sinlo;
//   return {-z_sinlo * xx.x + z_coslo * xx.y,
//           -z_sinlacoslo * xx.x - z_sinlasinlo * xx.y + cos(position.y) * xx.z};
// }

// cartesian to lonlat conversion manually "transpiled" from FORTRAN using
//  ELEMENTAL FUNCTION cc2gc_sphere(p_x) result (p_pos)
// mo_math_utilities.f90
glm::dvec2 cartToLonLatIcon(glm::dvec3 cart) {
  double z_x = cart.x;
  double z_y = cart.y;
  double z_z = cart.z;

  double z_r = z_x * z_x + z_y * z_y;
  z_r = sqrt(z_r);

  glm::dvec2 ret;

  if(fabs(z_r) < std::numeric_limits<double>::epsilon()) {
    if(z_z > 0.) {
      ret.y = M_PI_2;
    } else {
      ret.y = -1. * M_PI_2;
    }
    ret.x = 0.;
  } else {
    ret.y = atan2(z_z, z_r);
    if(fabs(z_x) < std::numeric_limits<double>::epsilon()) {
      if(z_y >= 0) {
        ret.x = M_PI_2;
      } else {
        ret.x = -M_PI_2;
      }
    } else {
      ret.x = atan2(z_y, z_x);
    }
  }
  return ret;
}

// https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates
// this leads to problems in triangles that cross the poles since the orientation of north flips
glm::dvec2 myprojLTPC_NED(const glm::dvec3& pCart, const glm::dvec3& pRefCart) {
  // double phi = asin((pRefCart.z) / (RADIUS_SPHERE));
  // double lambda = atan2(pRefCart.y, pRefCart.x);
  glm::dvec2 refLL = cartToLonLatIcon(pRefCart);
  double lambda = refLL.x;
  double phi = refLL.y;

  // glm uses column major ordering
  glm::dmat3x3 R(-sin(phi) * cos(lambda), -sin(phi) * sin(lambda), cos(phi), -sin(lambda),
                 cos(lambda), 0, -cos(phi) * cos(lambda), -cos(phi) * sin(lambda), -sin(phi));

  glm::dvec3 ret = glm::transpose(R) * (pCart - pRefCart); // this is in north, east, down
  // assert(fabs(ret.z) < 1e-3);                              // expecting this to be small
  return {ret.y, ret.x}; // would rather use x = east, y = north
}

glm::dvec2 myprojLTPC_ENU(const glm::dvec3& pCart, const glm::dvec3& pRefCart) {
  // double phi = asin((pRefCart.z) / (RADIUS_SPHERE));
  // double lambda = atan2(pRefCart.y, pRefCart.x);
  glm::dvec2 refLL = cartToLonLatIcon(pRefCart);
  double lambda = refLL.x;
  double phi = refLL.y;

  // glm uses column major ordering
  glm::dmat3x3 R(-sin(lambda), -sin(phi) * cos(lambda), cos(phi) * cos(lambda), cos(lambda),
                 -sin(phi) * sin(lambda), cos(phi) * sin(lambda), 0., cos(phi), sin(phi));

  glm::dvec3 ret = R * (pRefCart - pCart); // this is in east, north, up
  // assert(fabs(ret.z) < 1e-3);                              // expecting this to be small
  return {ret.x, ret.y}; // do not need up
}

//===-----------------------------------------------------------------------------

int main(int argc, char const* argv[]) {
  // enable floating point exception
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  if(argc != 2) {
    std::cout << "intended use is\n" << argv[0] << " <mesh>.nc" << std::endl;
    return -1;
  }

  // reference level of fluid, make sure to chose this large enough, otherwise initial
  // splash may induce negative fluid height and crash the sim
  const double refHeight = 2.;

  // constants
  const double CFLconst = 0.05;
  const double Grav = -9.81;

  // use high frequency damping. original damping by Cea and Blade is heavily dissipative, hence the
  // damping can be modulated by a coefficient in this implementation
  const bool use_corrector = true;
  const double DampingCoeff = 0.02;

  // optional bed friction, manning coefficient of 0.01 is roughly equal to flow of water over
  // concrete
  const bool use_friction = true;
  const double ManningCoeff = 0.01;

  int k_size = 1;
  const int level = 0;

  // dump a whole bunch of debug output (meant to be visualized using Octave, but gnuplot and the
  // like will certainly work too)
  const bool dbg_out = false;
  const bool readMeshFromDisk = false;

  atlas::Mesh mesh = AtlasMeshFromNetCDFMinimal(argv[1]).value();
  atlas::mesh::actions::build_edges(mesh, atlas::util::Config("pole_edges", false));
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
  atlas::mesh::actions::build_element_to_edge_connectivity(mesh);

  std::vector<double> nxNetcdf = LoadFieldNetCDF(argv[1], "zonal_normal_primal_edge");
  std::vector<double> nyNetcdf = LoadFieldNetCDF(argv[1], "meridional_normal_primal_edge");
  auto [edge_orientation_cellNetcdf, vertexPerCell, ncells] =
      Load2DFieldNetCDF(argv[1], "orientation_of_normal");

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

  dumpMesh4Triplot(mesh, "shallow", xyz);

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
  auto [nz_F, nz] = MakeAtlasField("nz", mesh.edges().size());
  auto [nxLC_F, nxLC] = MakeAtlasField("nx", mesh.edges().size()); // normals in local coordinates
  auto [nyLC_F, nyLC] = MakeAtlasField("ny", mesh.edges().size()); // (centered at edge midpoint)
  auto [alpha_F, alpha] = MakeAtlasField("alpha", mesh.edges().size());

  // Geometrical factors on cells
  auto [A_F, A] = MakeAtlasField("A", mesh.cells().size());
  auto [edge_orientation_cell_F, edge_orientation_cell] =
      MakeAtlasSparseField("edge_orientation_cell", mesh.cells().size(), edgesPerCell);

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  // {
  //   FILE* fp = fopen("dbgNorm.txt", "w+");
  //   for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
  //     // if(edgeIdx == 23120) {
  //     //   printf("!\n");
  //     // }

  //     L(edgeIdx, level) = edgeLength(mesh, xyz, edgeIdx);
  //     glm::dvec3 n = primalNormal(mesh, xyz, edgeIdx);
  //     nx(edgeIdx, level) = n.x;
  //     ny(edgeIdx, level) = n.y;
  //     nz(edgeIdx, level) = n.z;

  //     const auto& conn = mesh.edges().cell_connectivity();
  //     glm::dvec3 c0 = cellCircumcenter(mesh, xyz, conn(edgeIdx, 0));

  //     glm::dvec3 p = edgeMidpoint(mesh, xyz, edgeIdx);
  //     auto toLonLat = [](glm::dvec3 cart) -> glm::dvec2 {
  //       return {atan2(cart.y, cart.x), asin(cart.z / RADIUS_SPHERE)};
  //     };
  //     glm::dvec2 pLL = toLonLat(p); // <- ICON version
  //     // glm::dvec2 pLL = toLonLat(c0);

  //     glm::dvec2 nrmLC = glm::normalize(cvec2gvec(pLL, n));
  //     auto [n1, n2] = cartEdge(mesh, xyz, edgeIdx);
  //     glm::dvec2 n1LC = cvec2gvec(pLL, n1);
  //     glm::dvec2 n2LC = cvec2gvec(pLL, n2);
  //     glm::dvec2 eLC = glm::normalize(n1LC - n2LC);
  //     double check = glm::dot(eLC, nrmLC);
  //     if(check > 1e-6) {
  //       printf("%f %f %f: %e\n", p.x, p.y, p.z, check);
  //     }
  //     nxLC(edgeIdx, level) = nrmLC.x;
  //     nyLC(edgeIdx, level) = nrmLC.y;
  //     fprintf(fp, "%f %f %f %e\n", n.x, n.y, n.z, check);
  //   }
  //   fclose(fp);
  // }

  {
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      L(edgeIdx, level) = edgeLength(mesh, xyz, edgeIdx);
      glm::dvec3 n = primalNormal(mesh, xyz, edgeIdx);
      nx(edgeIdx, level) = n.x;
      ny(edgeIdx, level) = n.y;
      nz(edgeIdx, level) = n.z;

      glm::dvec3 p = edgeMidpoint(mesh, xyz, edgeIdx);

      glm::dvec2 nrmLC = myprojLTPC_ENU(p + n, p);
      // nxLC(edgeIdx, level) = nrmLC.x;
      // nyLC(edgeIdx, level) = nrmLC.y;

      nxLC(edgeIdx, level) = nxNetcdf[edgeIdx];
      nyLC(edgeIdx, level) = nyNetcdf[edgeIdx];
    }
  }

  {
    FILE* fp = fopen("cc.txt", "w+");
    int printIter = 1;
    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      glm::dvec3 pos = cellCircumcenter(mesh, xyz, cellIdx);
      if(pos.z < -0.8) {
        fprintf(fp, "%e %e %e\n", pos.x, pos.y, pos.z);
        printIter++;
      }
      if(cellIdx == 15360) {
        printf("here!\n");
      }
    }
    fclose(fp);
  }

  {
    FILE* fp = fopen("nrm.txt", "w+");
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      glm::dvec3 pos = edgeMidpoint(mesh, xyz, edgeIdx);
      fprintf(fp, "%f %f %f %f %f %f\n", pos.x, pos.y, pos.z, nx(edgeIdx, level),
              ny(edgeIdx, level), nz(edgeIdx, level));
    }
    fclose(fp);
  }

  {
    const auto& conn = mesh.edges().cell_connectivity();
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      int cellIdx1 = conn(edgeIdx, 0);
      int cellIdx2 = conn(edgeIdx, 1);
      assert(cellIdx1 != conn.missing_value());
      assert(cellIdx2 != conn.missing_value());
      double d1 = distanceToCircumcenter(mesh, xyz, cellIdx1, edgeIdx);
      double d2 = distanceToCircumcenter(mesh, xyz, cellIdx2, edgeIdx);
      alpha(edgeIdx, level) = d2 / (d1 + d2);
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on cells
  //===------------------------------------------------------------------------------------------===//
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    A(cellIdx, level) = cellArea(mesh, xyz, cellIdx);
  }

  dumpCellField("A", mesh, A, level);

  // for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
  //   const atlas::mesh::HybridElements::Connectivity& cellEdgeConnectivity =
  //       mesh.cells().edge_connectivity();

  //   glm::dvec3 c = cellCircumcenter(mesh, xyz, cellIdx);
  //   auto toLonLat = [](glm::dvec3 cart) -> glm::dvec2 {
  //     return {atan2(cart.y, cart.x), asin(cart.z / RADIUS_SPHERE)};
  //   };
  //   glm::dvec2 cLL = toLonLat(c);
  //   glm::dvec2 cLC = cvec2gvec(cLL, c);

  //   const int missingVal = cellEdgeConnectivity.missing_value();
  //   int numNbh = cellEdgeConnectivity.cols(cellIdx);
  //   assert(numNbh == edgesPerCell);

  //   if(cellIdx == 300) {
  //     const auto& conn = mesh.cells().node_connectivity();
  //     glm::dvec3 p1 = xyz[conn(cellIdx, 0)];
  //     glm::dvec3 p2 = xyz[conn(cellIdx, 1)];
  //     glm::dvec3 p3 = xyz[conn(cellIdx, 2)];

  //     glm::dvec2 p1LC = cvec2gvec(cLL, p1);
  //     glm::dvec2 p2LC = cvec2gvec(cLL, p2);
  //     glm::dvec2 p3LC = cvec2gvec(cLL, p3);

  //     printf("%f %f\n", cLC.x, cLC.y);
  //     printf("%f %f %f %f\n", p1LC.x, p2LC.x, p3LC.x, p1LC.x);
  //     printf("%f %f %f %f\n", p1LC.y, p2LC.y, p3LC.y, p1LC.y);

  //     // printf("%f %f %f %f\n", p1.x, p2.x, p3.x, p1.x);
  //     // printf("%f %f %f %f\n", p1.y, p2.y, p3.y, p1.y);
  //     // printf("%f %f %f %f\n", p1.z, p2.z, p3.z, p1.z);

  //     // cvec2gvec(em, em - c)
  //   }

  //   for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
  //     int edgeIdx = cellEdgeConnectivity(cellIdx, nbhIdx);

  //     glm::dvec3 em = edgeMidpoint(mesh, xyz, edgeIdx);
  //     glm::dvec2 emLC = cvec2gvec(cLL, em);
  //     // glm::dvec2 emLC =
  //     //     cc2tv(cLL, em); // guaranteed to point to outside since circumcenter is origin
  //     // glm::dvec2 nrmLC = cc2tv(cLL, {nx(edgeIdx, level), ny(edgeIdx, level), nz(edgeIdx,
  //     // level)});

  //     // glm::dvec2 toOutsideLC =
  //     //     cvec2gvec(em, em - c); // guaranteed to point to outside since circumcenter is
  //     //     origin
  //     // glm::dvec2 nrmLC =
  //     //     cvec2gvec(em, {nx(edgeIdx, level), ny(edgeIdx, level), nz(edgeIdx, level)});

  //     // edge_orientation_cell(cellIdx, nbhIdx, level) =
  //     //     sgn(glm::dot(glm::normalize(toOutsideLC), glm::normalize(nrmLC)));

  //     edge_orientation_cell(cellIdx, nbhIdx, level) = sgn(glm::dot(
  //         glm::normalize(em - c),
  //         glm::normalize(glm::dvec3{nx(edgeIdx, level), ny(edgeIdx, level), nz(edgeIdx,
  //         level)})));

  //     if(cellIdx == 300) {
  //       // printf("%f %f %f %f %f %f\n", em.x, em.y, em.z,
  //       //        edge_orientation_cell(cellIdx, nbhIdx, level) * nx(edgeIdx, level),
  //       //        edge_orientation_cell(cellIdx, nbhIdx, level) * ny(edgeIdx, level),
  //       //        edge_orientation_cell(cellIdx, nbhIdx, level) * nz(edgeIdx, level));

  //       printf("%f %f %f %f\n", emLC.x, emLC.y, nxLC(edgeIdx, level), nyLC(edgeIdx, level));

  //       // printf("%f %f %f %f %f %f\n", em.x, em.y, em.z, nx(edgeIdx, level), ny(edgeIdx,
  //       level),
  //       //        nz(edgeIdx, level));
  //     }
  //   }
  //   // explanation: the vector cellMidpoint -> edgeMidpoint is guaranteed to point outside. The
  //   // dot product checks if the edge normal has the same orientation. edgeMidpoint is arbitrary,
  //   // any point on e would work just as well
  // }

  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    const atlas::mesh::HybridElements::Connectivity& cellEdgeConnectivity =
        mesh.cells().edge_connectivity();

    glm::dvec3 c = cellCircumcenter(mesh, xyz, cellIdx);

    const int missingVal = cellEdgeConnectivity.missing_value();
    int numNbh = cellEdgeConnectivity.cols(cellIdx);
    assert(numNbh == edgesPerCell);

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = cellEdgeConnectivity(cellIdx, nbhIdx);

      glm::dvec3 em = edgeMidpoint(mesh, xyz, edgeIdx);
      if(cellIdx == 15360) {
        printf("%f %f %f\n", em.x, em.y, em.z);
      }

      // taking em as reference => normal and toInside have same reference / are in the same coord
      // system
      glm::dvec2 toInsideLC = myprojLTPC_ENU(c, em);
      glm::dvec2 n{nxLC(edgeIdx, level), nyLC(edgeIdx, level)};

      edge_orientation_cell(cellIdx, nbhIdx, level) =
          edge_orientation_cellNetcdf[nbhIdx * ncells + cellIdx];
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize height and other fields
  //===------------------------------------------------------------------------------------------===//
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    // auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
    // xm -= 1;
    // ym -= 1;
    // double v = sqrt(xm * xm + ym * ym);
    // h(cellIdx, level) = exp(-5 * v * v) + refHeight;
    h(cellIdx, level) = refHeight;
    // h(cellIdx, level) = sin(xm) * sin(ym) + refHeight;
  }

  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    qx(cellIdx, level) = 0.;
    qy(cellIdx, level) = 0.;
  }

  // {
  //   FILE* fp = fopen("ccExp.txt", "w+");
  //   for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
  //     glm::dvec2 c = cellCircumcenter2(mesh, cellIdx);
  //     double lon = c.x;
  //     double lat = c.y;
  //     double x = RADIUS_SPHERE * cos(lat) * cos(lon);
  //     double y = RADIUS_SPHERE * cos(lat) * sin(lon);
  //     double z = RADIUS_SPHERE * sin(lat);
  //     fprintf(fp, "%f %f %f\n", x, y, z);
  //   }
  //   fclose(fp);
  // }

  // {
  //   FILE* fp = fopen("cvec2gvec.txt", "w+");
  //   glm::dvec2 pole{0., 0.5 * M_PI};
  //   for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {

  //     glm::dvec3 pCart = edgeMidpoint(mesh, xyz, edgeIdx);
  //     glm::dvec3 nrmCart{nx(edgeIdx, level), ny(edgeIdx, level), nz(edgeIdx, level)};

  //     glm::dvec2 pos = cvec2gvec({0., 0.}, pCart);
  //     glm::dvec2 nrm = cvec2gvec({0., 0.}, nrmCart);

  //     fprintf(fp, "%f %f %f %f\n", pos.x, pos.y, nrm.x, nrm.y);
  //   }
  //   fclose(fp);
  // }

  // {
  //   FILE* fp = fopen("expP.txt", "w+");
  //   auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
  //   for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {

  //     glm::dvec2 pos = cvec2gvec({0., 0.}, xyz[nodeIdx]);

  //     fprintf(fp, "%f %f\n", pos.x, pos.y);
  //   }
  //   fclose(fp);
  // }

  double t = 0.;
  double dt = 0;
  double t_final = 16.;
  int step = 0;

  // writing this intentionally close to generated code
  while(t < t_final) {
    // if(step > 0 && step % 1000 == 0) {
    //   for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    //     auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
    //     xm -= 1;
    //     ym -= 1;
    //     double v = sqrt(xm * xm + ym * ym);
    //     h(cellIdx, level) += exp(-5 * v * v);
    //   }
    // }

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
    // exit(0);

    // normal edge velocity
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      lambda(edgeIdx, level) =
          nxLC(edgeIdx, level) * Ux(edgeIdx, level) + nyLC(edgeIdx, level) * Uy(edgeIdx, level);
    }

    // dumpEdgeField("lambda", mesh, xyz, lambda, level);
    // dumpEdgeField("hs", mesh, xyz, hs, level);

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
        assert(innerCell);
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
        dhdt(cellIdx, level) = lhs / A(cellIdx, level);
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
          lhs -= hs(edgeIdx, level) * nxLC(edgeIdx, level) *
                 edge_orientation_cell(cellIdx, nbhIdx, level) * L(cellIdx, level);
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
          // if(cellIdx == 15360) {
          //   double a = hs(edgeIdx, level);
          //   double b = nyLC(edgeIdx, level);
          //   double c = edge_orientation_cell(cellIdx, nbhIdx, level);
          //   printf("!\n");
          // }
          lhs -= hs(edgeIdx, level) * nyLC(edgeIdx, level) *
                 edge_orientation_cell(cellIdx, nbhIdx, level) * L(cellIdx, level);
        }
        Sy(cellIdx, level) = lhs / A(cellIdx, level);
      }
    }
    // dumpEdgeField("hs", mesh, wrapper, hs, level);
    // dumpCellField("Sx", mesh, wrapper, Sx, level);
    // dumpCellField("Sy", mesh, wrapper, Sy, level);
    // exit(0);

    // dumpCellField("dhdt", mesh, dhdt, level);
    // dumpCellField("dqxdt", mesh, dqxdt, level);
    // dumpCellField("dqydt", mesh, dqydt, level);
    // dumpCellField("Sx", mesh, Sx, level);
    // dumpCellField("Sy", mesh, Sy, level);
    // exit(0);

    for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      dhdt(cellIdx, level) = dhdt(cellIdx, level) * dt;
      dqxdt(cellIdx, level) =
          (dqxdt(cellIdx, level) - Grav * (h(cellIdx, level)) * Sx(cellIdx, level)) * dt;
      dqydt(cellIdx, level) =
          (dqydt(cellIdx, level) - Grav * (h(cellIdx, level)) * Sy(cellIdx, level)) * dt;
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

    if(step % 50 == 0) {
      char buf[256];
      // sprintf(buf, "out/step_%04d.txt", step);
      sprintf(buf, "out/stepH_%04d.txt", step);
      dumpCellField(buf, mesh, h, level);
    }
    std::cout << "time " << t << " timestep " << step++ << " dt " << dt << "\n";
  }
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
      // if(xyz[nodeIdx0 - 1].z < -0.8 && xyz[nodeIdx1 - 1].z < -0.8 && xyz[nodeIdx2 - 1].z < -0.8)
      // {
      //   fprintf(fp, "%d %d %d\n", nodeIdx0, nodeIdx1, nodeIdx2);
      // }
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
