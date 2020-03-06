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
// mo_math_divrot.f90) in ICON, using the Atlas, dawn, and the atlas interface therein. Some notes:
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

#include <cmath>
#include <cstdio>
#include <fenv.h>
#include <optional>
#include <vector>

// atlas functions
#include <atlas/array.h>
#include <atlas/grid.h>
#include <atlas/mesh.h>
#include <atlas/mesh/actions/BuildEdges.h>
// #include <atlas/mesh/actions/BuildEdges.h>
#include <atlas/util/CoordinateEnums.h>

// atlas interface for dawn generated code
#include "atlas_interface.hpp"

// icon stencil
#include "generated_iconLaplace.hpp"

// atlas utilities
#include "AtlasCartesianWrapper.h"
#include "AtlasFromNetcdf.h"
#include "GenerateRectAtlasMesh.h"

namespace {
//===------------------------------------------------------------------------------------------===//
// output (debugging)
//===------------------------------------------------------------------------------------------===//
template <typename T>
static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

void dumpMesh(const atlas::Mesh& m, AtlasToCartesian& wrapper, const std::string& fname);
void dumpDualMesh(const atlas::Mesh& m, AtlasToCartesian& wrapper, const std::string& fname);
void dumpMesh4Triplot(const atlas::Mesh& mesh, const std::string prefix,
                      std::optional<AtlasToCartesian> wrapper = std::nullopt);
void dumpNodeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level);
void dumpCellField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level,
                   std::optional<Orientation> color = std::nullopt);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level, std::vector<int> edgeList,
                   std::optional<Orientation> color = std::nullopt);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field_x, atlasInterface::Field<double>& field_y,
                   int level, std::optional<Orientation> color = std::nullopt);
//===------------------------------------------------------------------------------------------===//
// error reporting
//===------------------------------------------------------------------------------------------===//
std::tuple<double, double, double> MeasureErrors(std::vector<int> indices,
                                                 const atlasInterface::Field<double>& ref,
                                                 const atlasInterface::Field<double>& sol,
                                                 int level) {
  double Linf = 0.;
  double L1 = 0.;
  double L2 = 0.;
  for(int idx : indices) {
    double dif = ref(idx, level) - sol(idx, level);
    Linf = fmax(fabs(dif), Linf);
    L1 += fabs(dif);
    L2 += dif * dif;
  }
  L1 /= indices.size();
  L2 = sqrt(L2) / sqrt(indices.size());
  return {Linf, L1, L2};
}

} // namespace

int main(int argc, char const* argv[]) {
  // enable floating point exception
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  if(argc != 2) {
    std::cout << "intended use is\n" << argv[0] << " ny" << std::endl;
    return -1;
  }
  int w = atoi(argv[1]);
  int k_size = 1;
  const int level = 0;
  double lDomain = M_PI;

  // dump a whole bunch of debug output (meant to be visualized using Octave, but gnuplot and the
  // like will certainly work too)
  const bool dbg_out = false;
  const bool readMeshFromDisk = false;

  atlas::Mesh mesh;
  if(!readMeshFromDisk) {
    mesh = AtlasMeshRect(w);
    atlas::mesh::actions::build_edges(mesh, atlas::util::Config("pole_edges", false));
    atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
    atlas::mesh::actions::build_element_to_edge_connectivity(mesh);
  } else {
    mesh = AtlasMeshFromNetCDFComplete("testCaseMesh.nc").value();
    {
      auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
      auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
      for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
        xy(nodeIdx, atlas::LON) = lonlat(nodeIdx, atlas::LON);
        xy(nodeIdx, atlas::LAT) = lonlat(nodeIdx, atlas::LAT);
      }
    }
  }

  // wrapper with various atlas helper functions
  AtlasToCartesian wrapper(mesh);

  if(dbg_out) {
    dumpMesh4Triplot(mesh, "laplICONatlas_Mesh", wrapper);
  }

  const int edgesPerVertex = 6;
  const int edgesPerCell = 3;

  // current atlas mesh is not compatible with parallel computing
  // atlas::functionspace::CellColumns fs_cells(mesh, atlas::option::levels(k_size));
  // atlas::functionspace::NodeColumns fs_nodes(mesh, atlas::option::levels(k_size));
  // atlas::functionspace::EdgeColumns fs_edges(mesh, atlas::option::levels(k_size));

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

  //===------------------------------------------------------------------------------------------===//
  // input field (field we want to take the laplacian of)
  //===------------------------------------------------------------------------------------------===//
  auto [vec_F, vec] = MakeAtlasField("vec", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // control field holding the analytical solution for the divergence
  //===------------------------------------------------------------------------------------------===//
  auto [divVecSol_F, divVecSol] = MakeAtlasField("divVecSol", mesh.cells().size());

  //===------------------------------------------------------------------------------------------===//
  // control field holding the analytical solution for the curl
  //===------------------------------------------------------------------------------------------===//
  auto [rotVecSol_F, rotVecSol] = MakeAtlasField("rotVecSol", mesh.nodes().size());

  //===------------------------------------------------------------------------------------------===//
  // control field holding the analytical solution for Laplacian
  //===------------------------------------------------------------------------------------------===//
  auto [lapVecSol_F, lapVecSol] = MakeAtlasField("lapVecSol", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // output field (field containing the computed laplacian)
  //===------------------------------------------------------------------------------------------===//
  auto [nabla2_vec_F, nabla2_vec] = MakeAtlasField("nabla2_vec", mesh.edges().size());
  // term 1 and term 2 of nabla for debugging
  auto [nabla2t1_vec_F, nabla2t1_vec] = MakeAtlasField("nabla2t1_vec", mesh.edges().size());
  auto [nabla2t2_vec_F, nabla2t2_vec] = MakeAtlasField("nabla2t2_vec", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // intermediary fields (curl/rot and div of vec_e)
  //===------------------------------------------------------------------------------------------===//

  // rotation (more commonly curl) of vec_e on vertices
  auto [rot_vec_F, rot_vec] = MakeAtlasField("nabla2t2_vec", mesh.nodes().size());

  // divergence of vec_e on cells
  auto [div_vec_F, div_vec] = MakeAtlasField("nabla2t2_vec", mesh.cells().size());

  //===------------------------------------------------------------------------------------------===//
  // sparse dimensions for computing intermediary fields
  //===------------------------------------------------------------------------------------------===//

  // needed for the computation of the curl/rotation. according to documentation this needs to be:
  //
  // ! the appropriate dual cell based verts%edge_orientation
  // ! is required to obtain the correct value for the
  // ! application of Stokes theorem (which requires the scalar
  // ! product of the vector field with the tangent unit vectors
  // ! going around dual cell jv COUNTERCLOKWISE;
  // ! since the positive direction for the vec_e components is
  // ! not necessarily the one yelding counterclockwise rotation
  // ! around dual cell jv, a correction coefficient (equal to +-1)
  // ! is necessary, given by g%verts%edge_orientation
  auto [geofac_rot_F, geofac_rot] =
      MakeAtlasSparseField("geofac_rot", mesh.nodes().size(), edgesPerVertex);

  auto [edge_orientation_vertex_F, edge_orientation_vertex] =
      MakeAtlasSparseField("edge_orientation_vertex", mesh.nodes().size(), edgesPerVertex);

  // needed for the computation of the curl/rotation. according to documentation this needs to be:
  //
  //   ! ...the appropriate cell based edge_orientation is required to
  //   ! obtain the correct value for the application of Gauss theorem
  //   ! (which requires the scalar product of the vector field with the
  //   ! OUTWARD pointing unit vector with respect to cell jc; since the
  //   ! positive direction for the vector components is not necessarily
  //   ! the outward pointing one with respect to cell jc, a correction
  //   ! coefficient (equal to +-1) is necessary, given by
  //   ! ptr_patch%grid%cells%edge_orientation)
  auto [geofac_div_F, geofac_div] =
      MakeAtlasSparseField("geofac_div", mesh.cells().size(), edgesPerVertex);

  auto [edge_orientation_cell_F, edge_orientation_cell] =
      MakeAtlasSparseField("edge_orientation_cell", mesh.cells().size(), edgesPerCell);

  //===------------------------------------------------------------------------------------------===//
  // fields containing geometric information
  //===------------------------------------------------------------------------------------------===//
  auto [tangent_orientation_F, tangent_orientation] =
      MakeAtlasField("tangent_orientation", mesh.edges().size());
  auto [primal_edge_length_F, primal_edge_length] =
      MakeAtlasField("primal_edge_length", mesh.edges().size());
  auto [dual_edge_length_F, dual_edge_length] =
      MakeAtlasField("dual_edge_length", mesh.edges().size());
  auto [primal_normal_x_F, primal_normal_x] =
      MakeAtlasField("primal_normal_x", mesh.edges().size());
  auto [primal_normal_y_F, primal_normal_y] =
      MakeAtlasField("primal_normal_y", mesh.edges().size());
  auto [dual_normal_x_F, dual_normal_x] = MakeAtlasField("dual_normal_x", mesh.edges().size());
  auto [dual_normal_y_F, dual_normal_y] = MakeAtlasField("dual_normal_y", mesh.edges().size());
  auto [cell_area_F, cell_area] = MakeAtlasField("cell_area", mesh.cells().size());
  auto [dual_cell_area_F, dual_cell_area] = MakeAtlasField("dual_cell_area", mesh.nodes().size());

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    primal_edge_length(edgeIdx, level) = wrapper.edgeLength(mesh, edgeIdx);
    dual_edge_length(edgeIdx, level) = wrapper.dualEdgeLength(mesh, edgeIdx);
    tangent_orientation(edgeIdx, level) = wrapper.tangentOrientation(mesh, edgeIdx);
    auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);
    primal_normal_x(edgeIdx, level) = nx;
    primal_normal_y(edgeIdx, level) = ny;
    // The primal normal, dual normal
    // forms a left-handed coordinate system
    dual_normal_x(edgeIdx, level) = ny;
    dual_normal_y(edgeIdx, level) = -nx;
  }

  if(dbg_out) {
    dumpEdgeField("laplICONatlas_tangentOrientation.txt", mesh, wrapper, tangent_orientation,
                  level);
    dumpEdgeField("laplICONatlas_EdgeLength.txt", mesh, wrapper, primal_edge_length, level);
    dumpEdgeField("laplICONatlas_dualEdgeLength.txt", mesh, wrapper, dual_edge_length, level);
    dumpEdgeField("laplICONatlas_nrm.txt", mesh, wrapper, primal_normal_x, primal_normal_y, level);
    dumpEdgeField("laplICONatlas_dnrm.txt", mesh, wrapper, dual_normal_x, dual_normal_y, level);
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on cells
  //===------------------------------------------------------------------------------------------===//
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    cell_area(cellIdx, level) = wrapper.cellArea(mesh, cellIdx);
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on vertices
  //===------------------------------------------------------------------------------------------===//
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    dual_cell_area(nodeIdx, level) = wrapper.dualCellArea(mesh, nodeIdx);
  }

  if(dbg_out) {
    dumpCellField("laplICONatlas_areaCell.txt", mesh, wrapper, cell_area, level);
    dumpNodeField("laplICONatlas_areaCellDual.txt", mesh, wrapper, dual_cell_area, level);
  }

  //===------------------------------------------------------------------------------------------===//
  // input (spherical harmonics) and analytical solutions for div, curl and Laplacian
  //===------------------------------------------------------------------------------------------===//

  auto sphericalHarmonic = [](double x, double y) -> std::tuple<double, double> {
    return {0.25 * sqrt(105. / (2 * M_PI)) * cos(2 * x) * cos(y) * cos(y) * sin(y),
            0.5 * sqrt(15. / (2 * M_PI)) * cos(x) * cos(y) * sin(y)};
  };
  auto analyticalDivergence = [](double x, double y) {
    return -0.5 * (sqrt(105. / (2 * M_PI))) * sin(2 * x) * cos(y) * cos(y) * sin(y) +
           0.5 * sqrt(15. / (2 * M_PI)) * cos(x) * (cos(y) * cos(y) - sin(y) * sin(y));
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

  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    auto [xm, ym] = wrapper.edgeMidpoint(mesh, edgeIdx);
    auto [u, v] = sphericalHarmonic(xm, ym);
    auto [lu, lv] = analyticalLaplacian(xm, ym);
    vec(edgeIdx, level) = primal_normal_x(edgeIdx, level) * u + primal_normal_y(edgeIdx, level) * v;
    lapVecSol(edgeIdx, level) =
        primal_normal_x(edgeIdx, level) * lu + primal_normal_y(edgeIdx, level) * lv;
  }
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    auto [xm, ym] = wrapper.cellMidpoint(mesh, cellIdx);
    divVecSol(cellIdx, level) = analyticalDivergence(xm, ym);
  }
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    auto [xm, ym] = wrapper.nodeLocation(nodeIdx);
    rotVecSol(nodeIdx, level) = analyticalCurl(xm, ym);
  }

  //===------------------------------------------------------------------------------------------===//
  // Init geometrical factors (sparse fields)
  //===------------------------------------------------------------------------------------------===//

  // init edge orientations for vertices and cells
  auto dot = [](const Vector& v1, const Vector& v2) {
    return std::get<0>(v1) * std::get<0>(v2) + std::get<1>(v1) * std::get<1>(v2);
  };

  // Here, the ICON documentation states confusingly enough:
  //
  // +1 when the vector from this to the neigh-
  // bor vertex has the same orientation as the
  // tangent unit vector of the connecting edge.
  // -1 otherwise
  //
  // what this is really supposed to achieve is to compute a sparse dimension which ensures
  // that the normal and the orientation of the edge form a left handed coordinate system
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    const auto& nodeEdgeConnectivity = mesh.nodes().edge_connectivity();
    const auto& edgeNodeConnectivity = mesh.edges().node_connectivity();

    const int missingVal = nodeEdgeConnectivity.missing_value();
    int numNbh = nodeEdgeConnectivity.cols(nodeIdx);

    // arbitrary val at boundary
    bool anyMissing = false;
    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      anyMissing |= nodeEdgeConnectivity(nodeIdx, nbhIdx) == missingVal;
    }
    if(numNbh != 6 || anyMissing) {
      for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
        edge_orientation_vertex(nodeIdx, nbhIdx, level) = -1;
      }
      continue;
    }

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = nodeEdgeConnectivity(nodeIdx, nbhIdx);

      int n0 = edgeNodeConnectivity(edgeIdx, 0);
      int n1 = edgeNodeConnectivity(edgeIdx, 1);

      int centerIdx = (n0 == nodeIdx) ? n0 : n1;
      int farIdx = (n0 == nodeIdx) ? n1 : n0;

      auto [xLo, yLo] = wrapper.nodeLocation(centerIdx);
      auto [xHi, yHi] = wrapper.nodeLocation(farIdx);

      Vector edge = {xHi - xLo, yHi - yLo};
      Vector dualNormal = {dual_normal_x(edgeIdx, level), dual_normal_y(edgeIdx, level)};

      double dbg = dot(edge, dualNormal);
      int systemSign = sgn(dot(edge, dualNormal)); // geometrical factor "corrects" normal such that
                                                   // the resulting system is left handed
      edge_orientation_vertex(nodeIdx, nbhIdx, level) = systemSign;
    }
  }

  // ICON documentation states
  //
  // The orientation of the edge normal vector
  // (the variable primal normal in the edges ta-
  // ble) for the cell according to Gauss formula.
  // It is equal to +1 if the normal to the edge
  // is outwards from the cell, otherwise is -1.
  //
  // which is quite clear, here goes:
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
      Vector primal = {primal_normal_x(edgeIdx, level), primal_normal_y(edgeIdx, level)};
      edge_orientation_cell(cellIdx, nbhIdx, level) = sgn(dot(toOutsdie, primal));
    }
    // explanation: the vector cellMidpoint -> edgeMidpoint is guaranteed to point outside. The
    // dot product checks if the edge normal has the same orientation. edgeMidpoint is arbitrary,
    // any point on e would work just as well
  }

  // now, consume these two "orientation" fields to form the actual geometrical factors (which
  // include information about the meshes edge lengths and cell areas)
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    const atlas::mesh::Nodes::Connectivity& nodeEdgeConnectivity = mesh.nodes().edge_connectivity();

    int numNbh = nodeEdgeConnectivity.cols(nodeIdx);

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = nodeEdgeConnectivity(nodeIdx, nbhIdx);
      geofac_rot(nodeIdx, nbhIdx, level) =
          (dual_cell_area(nodeIdx, level) == 0.)
              ? 0
              : dual_edge_length(edgeIdx, level) * edge_orientation_vertex(nodeIdx, nbhIdx, level) /
                    dual_cell_area(nodeIdx, level);
    }
    // Original ICON code
    //
    // ptr_int%geofac_rot(jv,je,jb) =                &
    //    & ptr_patch%edges%dual_edge_length(ile,ibe) * &
    //    & ptr_patch%verts%edge_orientation(jv,jb,je)/ &
    //    & ptr_patch%verts%dual_area(jv,jb) * REAL(ifac,wp)
  }

  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    const atlas::mesh::HybridElements::Connectivity& cellEdgeConnectivity =
        mesh.cells().edge_connectivity();

    int numNbh = cellEdgeConnectivity.cols(cellIdx);
    assert(numNbh == edgesPerCell);

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = cellEdgeConnectivity(cellIdx, nbhIdx);
      geofac_div(cellIdx, nbhIdx, level) = primal_edge_length(edgeIdx, level) *
                                           edge_orientation_cell(cellIdx, nbhIdx, level) /
                                           cell_area(cellIdx, level);
    }
    // Original ICON code
    //
    //  ptr_int%geofac_div(jc,je,jb) = &
    //    & ptr_patch%edges%primal_edge_length(ile,ibe) * &
    //    & ptr_patch%cells%edge_orientation(jc,jb,je)  / &
    //    & ptr_patch%cells%area(jc,jb)
  }

  //===------------------------------------------------------------------------------------------===//
  // stencil call
  //===------------------------------------------------------------------------------------------===/
  dawn_generated::cxxnaiveico::icon<atlasInterface::atlasTag>(
      mesh, k_size, vec, div_vec, rot_vec, nabla2t1_vec, nabla2t2_vec, nabla2_vec,
      primal_edge_length, dual_edge_length, tangent_orientation, geofac_rot, geofac_div)
      .run();

  if(dbg_out) {
    dumpEdgeField("laplICONatlas_nabla2t1.txt", mesh, wrapper, nabla2t1_vec, level,
                  wrapper.innerEdges(mesh));
    dumpEdgeField("laplICONatlas_nabla2t2.txt", mesh, wrapper, nabla2t1_vec, level,
                  wrapper.innerEdges(mesh));
  }

  //===------------------------------------------------------------------------------------------===//
  // dumping a hopefully nice colorful divergence, curl & laplacian
  //===------------------------------------------------------------------------------------------===//
  dumpCellField("laplICONatlas_div.txt", mesh, wrapper, div_vec, level);
  dumpNodeField("laplICONatlas_rot.txt", mesh, wrapper, rot_vec, level);
  dumpEdgeField("laplICONatlas_out.txt", mesh, wrapper, nabla2_vec, level,
                wrapper.innerEdges(mesh));

  //===------------------------------------------------------------------------------------------===//
  // measuring errors
  //===------------------------------------------------------------------------------------------===//
  {
    auto [Linf, L1, L2] = MeasureErrors(wrapper.innerCells(mesh), divVecSol, div_vec, level);
    printf("[div] dx: %e L_inf: %e L_1: %e L_2: %e\n", 180. / w, Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureErrors(wrapper.innerNodes(mesh), rotVecSol, rot_vec, level);
    printf("[rot] dx: %e L_inf: %e L_1: %e L_2: %e\n", 180. / w, Linf, L1, L2);
  }
  {
    auto [Linf, L1, L2] = MeasureErrors(wrapper.innerEdges(mesh), lapVecSol, nabla2_vec, level);
    printf("[lap] dx: %e L_inf: %e L_1: %e L_2: %e\n", 180. / w, Linf, L1, L2);
  }

  printf("----\n");

  return 0;
}

namespace {

void dumpMesh4Triplot(const atlas::Mesh& mesh, const std::string prefix,
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
}

void dumpMesh(const atlas::Mesh& mesh, AtlasToCartesian& wrapper, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  const atlas::mesh::HybridElements::Connectivity& edgeNodeConnectivity =
      mesh.edges().node_connectivity();
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    int numNbh = edgeNodeConnectivity.cols(edgeIdx);
    assert(numNbh == 2);

    int nbhLo = edgeNodeConnectivity(edgeIdx, 0);
    int nbhHi = edgeNodeConnectivity(edgeIdx, 1);

    auto [xLo, yLo] = wrapper.nodeLocation(nbhLo);
    auto [xHi, yHi] = wrapper.nodeLocation(nbhHi);

    fprintf(fp, "%f %f %f %f\n", xLo, yLo, xHi, yHi);
  }
  fclose(fp);
}

void dumpDualMesh(const atlas::Mesh& mesh, AtlasToCartesian& wrapper, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  const atlas::mesh::HybridElements::Connectivity& edgeCellConnectivity =
      mesh.edges().cell_connectivity();
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {

    int nbhLo = edgeCellConnectivity(edgeIdx, 0);
    int nbhHi = edgeCellConnectivity(edgeIdx, 1);

    if(nbhLo == edgeCellConnectivity.missing_value() ||
       nbhHi == edgeCellConnectivity.missing_value()) {
      continue;
    }

    auto [xm1, ym1] = wrapper.cellCircumcenter(mesh, nbhLo);
    auto [xm2, ym2] = wrapper.cellCircumcenter(mesh, nbhHi);

    fprintf(fp, "%f %f %f %f\n", xm1, ym1, xm2, ym2);
  }
  fclose(fp);
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

void dumpCellField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
    fprintf(fp, "%f %f %f\n", xm, ym, field(cellIdx, level));
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
} // namespace