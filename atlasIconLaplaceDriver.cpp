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
//
//    WARNING! THIS IS A PROTOTYPE! DOES NOT YET PRODUCE CORRECT RESULTS! WARNING!
//
//
//===------------------------------------------------------------------------------------------===//

#include <cmath>
#include <vector>

#include "AtlasCartesianWrapper.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/library/Library.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/meshgenerator.h"

#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"

#include "atlas_interface.hpp"

#include "generated_iconLaplace.hpp"

// remove later
#include "atlas/output/Gmsh.h"

#include "AtlasFromNetcdf.h"
#include "GenerateRectAtlasMesh.h"

template <typename T>
static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

//===------------------------------------------------------------------------------------------===//
// output (debugging)
//===------------------------------------------------------------------------------------------===//
void dumpMesh(const atlas::Mesh& m, AtlasToCartesian& wrapper, const std::string& fname);
void dumpDualMesh(const atlas::Mesh& m, AtlasToCartesian& wrapper, const std::string& fname);

void dumpNodeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level);
void dumpCellField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field, int level,
                   std::optional<Orientation> color = std::nullopt);
void dumpEdgeField(const std::string& fname, const atlas::Mesh& mesh, AtlasToCartesian wrapper,
                   atlasInterface::Field<double>& field_x, atlasInterface::Field<double>& field_y,
                   int level, std::optional<Orientation> color = std::nullopt);

void debugDump(const atlas::Mesh& mesh, const std::string prefix) {
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

int main() {
  int w = 64;
  int k_size = 1;
  const int level = 0;
  double lDomain = M_PI;

  const bool dbg_out = false;
  const bool readMeshFromDisk = false;

  atlas::Mesh mesh;
  if(!readMeshFromDisk) {
    mesh = AtlasMeshRect(w);
    atlas::mesh::actions::build_edges(mesh, atlas::util::Config("pole_edges", false));
    atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
    atlas::mesh::actions::build_element_to_edge_connectivity(mesh);
    debugDump(mesh, "atlasMesh");
  } else {
    mesh = AtlasMeshFromNetCDFComplete("testCaseMesh.nc").value();
    {
      // auto radToLat = [](double rad) { return rad / (0.5 * M_PI) * 90; };
      // auto radToLon = [](double rad) { return rad / (M_PI)*180; };
      auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
      auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
      for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
        xy(nodeIdx, atlas::LON) = lonlat(nodeIdx, atlas::LON);
        xy(nodeIdx, atlas::LAT) = lonlat(nodeIdx, atlas::LAT);
      }
    }
    debugDump(mesh, "testCaseMesh");
  }

  const bool skewToEquilateral = true;
  AtlasToCartesian wrapper(mesh);

  if(dbg_out) {
    dumpMesh(mesh, wrapper, "laplICONatlas_mesh.txt");
    dumpDualMesh(mesh, wrapper, "laplICONatlas_dualMesh.txt");
  }

  const int edgesPerVertex = 6;
  const int edgesPerCell = 3;

  // current atlas mesh is not compatible with parallel computing
  // atlas::functionspace::CellColumns fs_cells(mesh, atlas::option::levels(k_size));
  // atlas::functionspace::NodeColumns fs_nodes(mesh, atlas::option::levels(k_size));
  // atlas::functionspace::EdgeColumns fs_edges(mesh, atlas::option::levels(k_size));

  //===------------------------------------------------------------------------------------------===//
  // input field (field we want to take the laplacian of)
  //===------------------------------------------------------------------------------------------===//
  atlas::Field vec_F{"vec", atlas::array::DataType::real64(),
                     atlas::array::make_shape(mesh.edges().size(), 1)};
  atlasInterface::Field<double> vec = atlas::array::make_view<double, 2>(vec_F);
  //    this is, confusingly, called vec_e, even though it is a scalar
  //    _conceptually_, this can be regarded as a vector with implicit direction (presumably
  //    normal to edge direction)

  //===------------------------------------------------------------------------------------------===//
  // control field holding the analytical solution for the divergence
  //===------------------------------------------------------------------------------------------===//
  atlas::Field divVecSol_F{"divVecSol", atlas::array::DataType::real64(),
                           atlas::array::make_shape(mesh.edges().size(), 1)};
  atlasInterface::Field<double> divVecSol = atlas::array::make_view<double, 2>(divVecSol_F);

  //===------------------------------------------------------------------------------------------===//
  // output field (field containing the computed laplacian)
  //===------------------------------------------------------------------------------------------===//
  atlas::Field nabla2_vec_F{"nabla2_vec", atlas::array::DataType::real64(),
                            atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> nabla2_vec = atlas::array::make_view<double, 2>(nabla2_vec_F);
  // term 1 and term 2 of nabla for debugging
  atlas::Field nabla2t1_vec_F{"nabla2t1_vec", atlas::array::DataType::real64(),
                              atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> nabla2t1_vec = atlas::array::make_view<double, 2>(nabla2t1_vec_F);
  atlas::Field nabla2t2_vec_F{"nabla2t2_vec", atlas::array::DataType::real64(),
                              atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> nabla2t2_vec = atlas::array::make_view<double, 2>(nabla2t2_vec_F);
  //    again, surprisingly enough, this is a scalar quantity even though the vector laplacian is
  //    a laplacian.

  //===------------------------------------------------------------------------------------------===//
  // intermediary fields (curl/rot and div of vec_e)
  //===------------------------------------------------------------------------------------------===//

  // rotation (more commonly curl) of vec_e on vertices
  //    I'm not entirely positive how one can take the curl of a scalar field (commonly a
  //    undefined operation), however, since vec_e is _conceptually_ a vector this works out.
  //    somehow.
  atlas::Field rot_vec_F{"rot_vec", atlas::array::DataType::real64(),
                         atlas::array::make_shape(mesh.nodes().size(), k_size)};
  atlasInterface::Field<double> rot_vec = atlas::array::make_view<double, 2>(rot_vec_F);

  // divergence of vec_e on cells
  //    Again, not entirely sure how one can measure the divergence of scalars, but again, vec_e
  //    is _conceptually_ a vector, so...
  atlas::Field div_vec_F{"div_vec", atlas::array::DataType::real64(),
                         atlas::array::make_shape(mesh.cells().size(), k_size)};
  atlasInterface::Field<double> div_vec = atlas::array::make_view<double, 2>(div_vec_F);

  //===------------------------------------------------------------------------------------------===//
  // sparse dimensions for computing intermediary fields
  //===------------------------------------------------------------------------------------------===//

  // needed for the computation of the curl/rotation. according to documentation this needs to be:
  // ! the appropriate dual cell based verts%edge_orientation
  // ! is required to obtain the correct value for the
  // ! application of Stokes theorem (which requires the scalar
  // ! product of the vector field with the tangent unit vectors
  // ! going around dual cell jv COUNTERCLOKWISE;
  // ! since the positive direction for the vec_e components is
  // ! not necessarily the one yelding counterclockwise rotation
  // ! around dual cell jv, a correction coefficient (equal to +-1)
  // ! is necessary, given by g%verts%edge_orientation

  atlas::Field geofac_rot_F{"geofac_rot", atlas::array::DataType::real64(),
                            atlas::array::make_shape(mesh.nodes().size(), k_size, edgesPerVertex)};
  atlasInterface::SparseDimension<double> geofac_rot =
      atlas::array::make_view<double, 3>(geofac_rot_F);
  atlas::Field edge_orientation_vertex_F{
      "edge_orientation_vertex", atlas::array::DataType::real64(),
      atlas::array::make_shape(mesh.nodes().size(), k_size, edgesPerVertex)};
  atlasInterface::SparseDimension<double> edge_orientation_vertex =
      atlas::array::make_view<double, 3>(edge_orientation_vertex_F);

  // needed for the computation of the curl/rotation. according to documentation this needs to be:
  //   ! ...the appropriate cell based edge_orientation is required to
  //   ! obtain the correct value for the application of Gauss theorem
  //   ! (which requires the scalar product of the vector field with the
  //   ! OUTWARD pointing unit vector with respect to cell jc; since the
  //   ! positive direction for the vector components is not necessarily
  //   ! the outward pointing one with respect to cell jc, a correction
  //   ! coefficient (equal to +-1) is necessary, given by
  //   ! ptr_patch%grid%cells%edge_orientation)

  atlas::Field geofac_div_F{"geofac_div", atlas::array::DataType::real64(),
                            atlas::array::make_shape(mesh.cells().size(), k_size, edgesPerCell)};
  atlasInterface::SparseDimension<double> geofac_div =
      atlas::array::make_view<double, 3>(geofac_div_F);
  atlas::Field edge_orientation_cell_F{
      "edge_orientation_cell", atlas::array::DataType::real64(),
      atlas::array::make_shape(mesh.cells().size(), k_size, edgesPerCell)};
  atlasInterface::SparseDimension<double> edge_orientation_cell =
      atlas::array::make_view<double, 3>(edge_orientation_cell_F);

  // //===------------------------------------------------------------------------------------------===//
  // // fields containing geometric information
  // //===------------------------------------------------------------------------------------------===//
  atlas::Field tangent_orientation_F{"tangent_orientation", atlas::array::DataType::real64(),
                                     atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> tangent_orientation =
      atlas::array::make_view<double, 2>(tangent_orientation_F);
  atlas::Field primal_edge_length_F{"primal_edge_length", atlas::array::DataType::real64(),
                                    atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> primal_edge_length =
      atlas::array::make_view<double, 2>(primal_edge_length_F);
  atlas::Field dual_edge_length_F{"dual_edge_length", atlas::array::DataType::real64(),
                                  atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> dual_edge_length =
      atlas::array::make_view<double, 2>(dual_edge_length_F);
  atlas::Field dual_normal_x_F{"dual_normal_x", atlas::array::DataType::real64(),
                               atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> dual_normal_x = atlas::array::make_view<double, 2>(dual_normal_x_F);
  atlas::Field dual_normal_y_F{"dual_normal_y", atlas::array::DataType::real64(),
                               atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> dual_normal_y = atlas::array::make_view<double, 2>(dual_normal_y_F);
  atlas::Field primal_normal_x_F{"primal_normal_x", atlas::array::DataType::real64(),
                                 atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> primal_normal_x =
      atlas::array::make_view<double, 2>(primal_normal_x_F);
  atlas::Field primal_normal_y_F{"primal_normal_y", atlas::array::DataType::real64(),
                                 atlas::array::make_shape(mesh.edges().size(), k_size)};
  atlasInterface::Field<double> primal_normal_y =
      atlas::array::make_view<double, 2>(primal_normal_y_F);

  atlas::Field cell_area_F{"cell_area", atlas::array::DataType::real64(),
                           atlas::array::make_shape(mesh.cells().size(), k_size)};
  atlasInterface::Field<double> cell_area = atlas::array::make_view<double, 2>(cell_area_F);

  atlas::Field dual_cell_area_F{"dual_cell_area", atlas::array::DataType::real64(),
                                atlas::array::make_shape(mesh.nodes().size(), k_size)};
  atlasInterface::Field<double> dual_cell_area =
      atlas::array::make_view<double, 2>(dual_cell_area_F);

  //===------------------------------------------------------------------------------------------===//
  // initialize fields
  //===------------------------------------------------------------------------------------------===//

  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    rot_vec(nodeIdx, level) = 0;
  }

  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    div_vec(cellIdx, level) = 0;
  }

  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    primal_edge_length(edgeIdx, level) = wrapper.edgeLength(mesh, edgeIdx);
    dual_edge_length(edgeIdx, level) = wrapper.dualEdgeLength(mesh, edgeIdx);
    tangent_orientation(edgeIdx, level) = wrapper.tangentOrientation(mesh, edgeIdx);
    auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);
    primal_normal_x(edgeIdx, level) = nx * tangent_orientation(edgeIdx, level);
    primal_normal_y(edgeIdx, level) = ny * tangent_orientation(edgeIdx, level);
    // The primal normal, dual normal
    // forms a left-handed coordinate system
    dual_normal_x(edgeIdx, level) = ny;
    dual_normal_y(edgeIdx, level) = -nx;
  }

  if(dbg_out) {
    dumpEdgeField("laplICONatlas_EdgeLength.txt", mesh, wrapper, primal_edge_length, level);
    dumpEdgeField("laplICONatlas_dualEdgeLength.txt", mesh, wrapper, dual_edge_length, level);
    dumpEdgeField("laplICONatlas_nrm.txt", mesh, wrapper, primal_normal_x, primal_normal_y, level);
    dumpEdgeField("laplICONatlas_dnrm.txt", mesh, wrapper, dual_normal_x, dual_normal_y, level);
  }

  auto wave = [](double px, double py) { return sin(px) * sin(py); };
  auto constant = [](double px, double py) { return 1.; };
  auto lin = [](double px, double py) { return px; };

  auto sphericalHarmonic = [](double x, double y) -> std::tuple<double, double> {
    return {0.25 * sqrt(105. / (2 * M_PI)) * cos(2 * x) * cos(y) * cos(y) * sin(y),
            0.5 * sqrt(15. / (2 * M_PI)) * cos(x) * cos(y) * sin(y)};
  };

  auto analyticalDivergence = [](double x, double y) {
    return 1. / (2 * sqrt(2 * M_PI)) *
           (sqrt(105) * sin(2 * x) * cos(y) * cos(y) * sin(y) + sqrt(15.) * cos(x) * cos(2 * y));
  };

  // init zero and test function
  FILE* fp = fopen("sphericalHarmonics.txt", "w+");
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    auto [xm, ym] = wrapper.edgeMidpoint(mesh, edgeIdx);
    auto [u, v] = sphericalHarmonic(xm, ym);
    double normalWind = primal_normal_x(edgeIdx, level) * u + primal_normal_y(edgeIdx, level) * v;
    vec(edgeIdx, level) = normalWind;
    divVecSol(edgeIdx, level) = analyticalDivergence(xm, ym);
    nabla2_vec(edgeIdx, level) = 0;
    fprintf(fp, "%f %f %f %f\n", xm, ym, u, v);
  }
  fclose(fp);

  // init geometric info for cells
  for(int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    cell_area(cellIdx, level) = wrapper.cellArea(mesh, cellIdx);
  }
  // init geometric info for vertices
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    dual_cell_area(nodeIdx, level) = wrapper.dualCellArea(mesh, nodeIdx);
  }

  if(dbg_out) {
    dumpEdgeField("laplICONatlas_divAnalytical.txt", mesh, wrapper, divVecSol, level);
    dumpCellField("laplICONatlas_areaCell.txt", mesh, wrapper, cell_area, level);
    dumpNodeField("laplICONatlas_areaCellDual.txt", mesh, wrapper, dual_cell_area, level);
  }

  // init edge orientations for vertices and cells
  auto dot = [](const Vector& v1, const Vector& v2) {
    return std::get<0>(v1) * std::get<0>(v2) + std::get<1>(v1) * std::get<1>(v2);
  };

  // +1 when the vector from this to the neigh-
  // bor vertex has the same orientation as the
  // tangent unit vector of the connecting edge.
  // -1 otherwise

  auto nodeNeighboursOfNode = [](atlas::Mesh const& m, int const& idx) {
    const auto& conn_nodes_to_edge = m.nodes().edge_connectivity();
    auto neighs = std::vector<std::tuple<int, int>>{};
    for(int ne = 0; ne < conn_nodes_to_edge.cols(idx); ++ne) {
      int nbh_edge_idx = conn_nodes_to_edge(idx, ne);
      const auto& conn_edge_to_nodes = m.edges().node_connectivity();
      for(int nn = 0; nn < conn_edge_to_nodes.cols(nbh_edge_idx); ++nn) {
        int nbhNode = conn_edge_to_nodes(idx, nn);
        if(nbhNode != idx) {
          neighs.emplace_back(std::tuple<int, int>(nbh_edge_idx, nbhNode));
        }
      }
    }
    return neighs;
  };

  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    const atlas::mesh::Nodes::Connectivity& nodeEdgeConnectivity = mesh.nodes().edge_connectivity();
    const atlas::mesh::HybridElements::Connectivity& edgeNodeConnectivity =
        mesh.edges().node_connectivity();

    const int missingVal = nodeEdgeConnectivity.missing_value();
    int numNbh = nodeEdgeConnectivity.cols(nodeIdx);

    if(numNbh != 6) {
      continue;
    }

    auto nbh = nodeNeighboursOfNode(mesh, nodeIdx);
    auto [cx, cy] = wrapper.nodeLocation(nodeIdx);

    if(nbh.size() != 6) {
      continue;
    }

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = std::get<0>(nbh[nbhIdx]);
      int farNodeIdx = std::get<1>(nbh[nbhIdx]);

      int nodeIdxLo = edgeNodeConnectivity(edgeIdx, 0);
      int nodeIdxHi = edgeNodeConnectivity(edgeIdx, 1);

      auto [xLo, yLo] = wrapper.nodeLocation(nodeIdxLo);
      auto [xHi, yHi] = wrapper.nodeLocation(nodeIdxHi);

      auto [farX, farY] = wrapper.nodeLocation(farNodeIdx);

      Vector edgeVec{xHi - xLo, yHi - yLo};
      // its not quite clear how to implement this properly in Atlas
      //        nodes have no node neighbors on an atlas grid
      //        Vector dualNrm{cx - farX, cy - farY}; <- leads to oscillations in rot field
      Vector dualNrm{dual_normal_x(edgeIdx, level), dual_normal_y(edgeIdx, level)};
      edge_orientation_vertex(nodeIdx, nbhIdx, level) = sgn(dot(edgeVec, dualNrm));
    }
  }

  // The orientation of the edge normal vector
  // (the variable primal normal in the edges ta-
  // ble) for the cell according to Gauss formula.
  // It is equal to +1 if the normal to the edge
  // is outwards from the cell, otherwise is -1.
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

  // init sparse quantities for div and rot
  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    const atlas::mesh::Nodes::Connectivity& nodeEdgeConnectivity = mesh.nodes().edge_connectivity();

    int numNbh = nodeEdgeConnectivity.cols(nodeIdx);
    // assert(numNbh == edgesPerVertex);

    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = nodeEdgeConnectivity(nodeIdx, nbhIdx);
      geofac_rot(nodeIdx, nbhIdx, level) = dual_edge_length(edgeIdx, level) *
                                           edge_orientation_vertex(nodeIdx, nbhIdx, level) /
                                           dual_cell_area(nodeIdx, level);
    }
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
    //  ptr_int%geofac_div(jc,je,jb) = &
    //    & ptr_patch%edges%primal_edge_length(ile,ibe) * &
    //    & ptr_patch%cells%edge_orientation(jc,jb,je)  / &
    //    & ptr_patch%cells%area(jc,jb)
  }

  //===------------------------------------------------------------------------------------------===//
  // stencil call
  //===------------------------------------------------------------------------------------------===//

  dawn_generated::cxxnaiveico::icon<atlasInterface::atlasTag>(
      mesh, k_size, vec, div_vec, rot_vec, nabla2t1_vec, nabla2t2_vec, nabla2_vec,
      primal_edge_length, dual_edge_length, tangent_orientation, geofac_rot, geofac_div)
      .run();

  if(dbg_out) {
    dumpCellField("laplICONatlas_div.txt", mesh, wrapper, div_vec, level);
    dumpNodeField("laplICONatlas_rot.txt", mesh, wrapper, rot_vec, level);

    dumpEdgeField("laplICONatlas_rotH.txt", mesh, wrapper, nabla2t1_vec, level,
                  Orientation::Horizontal);
    dumpEdgeField("laplICONatlas_rotV.txt", mesh, wrapper, nabla2t1_vec, level,
                  Orientation::Vertical);
    dumpEdgeField("laplICONatlas_rotD.txt", mesh, wrapper, nabla2t1_vec, level,
                  Orientation::Diagonal);

    dumpEdgeField("laplICONatlas_divH.txt", mesh, wrapper, nabla2t2_vec, level,
                  Orientation::Horizontal);
    dumpEdgeField("laplICONatlas_divV.txt", mesh, wrapper, nabla2t2_vec, level,
                  Orientation::Vertical);
    dumpEdgeField("laplICONatlas_divD.txt", mesh, wrapper, nabla2t2_vec, level,
                  Orientation::Diagonal);
  }

  //===------------------------------------------------------------------------------------------===//
  // dumping a hopefully nice colorful laplacian
  //===------------------------------------------------------------------------------------------===//
  dumpEdgeField("laplICONatlas_out.txt", mesh, wrapper, nabla2_vec, level);

  //===------------------------------------------------------------------------------------------===//
  // measuring errors
  //===------------------------------------------------------------------------------------------===//
  {
    double Linf = 0;
    double L1 = 0;
    double L2 = 0;
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      double diff = div_vec(edgeIdx, level) - divVecSol(edgeIdx, level);
      if(diff > 10.) { // check outliers
        continue;
      }
      Linf = fmax(Linf, diff);
      L1 += fabs(diff);
      L2 += diff * diff;
      // printf("%f\n", diff);
    }
    L1 /= mesh.edges().size();
    L2 /= mesh.edges().size();
    printf("errors divergence: Linf: %e, L1 %e, L2 %e\n", Linf, L1, L2);
  }

  return 0;
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
    fprintf(fp, "%f %f %f\n", xm, ym, field(edgeIdx, level));
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

//   auto nodeNeighboursOfNode = [](atlas::Mesh const& m, int const& idx) {
//     const auto& conn_nodes_to_edge = m.nodes().edge_connectivity();
//     auto neighs = std::vector<std::tuple<int, int>>{};
//     for(int ne = 0; ne < conn_nodes_to_edge.cols(idx); ++ne) {
//       int nbh_edge_idx = conn_nodes_to_edge(idx, ne);
//       const auto& conn_edge_to_nodes = m.edges().node_connectivity();
//       for(int nn = 0; nn < conn_edge_to_nodes.cols(nbh_edge_idx); ++nn) {
//         int nbhNode = conn_edge_to_nodes(idx, nn);
//         if(nbhNode != idx) {
//           neighs.emplace_back(std::tuple<int, int>(nbh_edge_idx, nbhNode));
//         }
//       }
//     }
//     return neighs;
//   };

//   for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
//     const atlas::mesh::Nodes::Connectivity& nodeEdgeConnectivity =
//     mesh.nodes().edge_connectivity(); const atlas::mesh::HybridElements::Connectivity&
//     edgeNodeConnectivity =
//         mesh.edges().node_connectivity();

//     const int missingVal = nodeEdgeConnectivity.missing_value();
//     int numNbh = nodeEdgeConnectivity.cols(nodeIdx);

//     if(numNbh != 6) {
//       continue;
//     }

//     auto [vx, vy] = wrapper.nodeLocation(nodeIdx);
//     auto nodeNbh = nodeNeighboursOfNode(mesh, nodeIdx);

//     int nbh_idx = 0;
//     for(const auto& neighbor : nodeNbh) {
//       int edgeIdx = std::get<0>(neighbor);
//       int nodeIdx = std::get<1>(neighbor);

//       auto [farVx, farVy] = wrapper.nodeLocation(nodeIdx);

//       Vector testVector{farVx - vx, farVy - vy};

//       Vector dual = {dual_normal_x(edgeIdx, level), dual_normal_y(edgeIdx, level)};
//       edge_orientation_vertex(nodeIdx, nbh_idx, level) = sgn(dot(testVector, dual));
//       nbh_idx++;
//     }
//   }