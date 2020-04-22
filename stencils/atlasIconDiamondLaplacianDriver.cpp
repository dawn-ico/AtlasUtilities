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
#include "interfaces/atlas_interface.hpp"

// icon stencil
#include "generated_iconDiamondLaplace.hpp"

// atlas utilities
#include "../utils/AtlasCartesianWrapper.h"
#include "../utils/AtlasFromNetcdf.h"
#include "../utils/GenerateRectAtlasMesh.h"

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
  mesh = AtlasMeshRect(w);
  atlas::mesh::actions::build_edges(mesh, atlas::util::Config("pole_edges", false));
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
  atlas::mesh::actions::build_element_to_edge_connectivity(mesh);

  // wrapper with various atlas helper functions
  AtlasToCartesian wrapper(mesh, true);

  const int edgesPerVertex = 6;
  const int edgesPerCell = 3;
  const int verticesInDiamond = 4;

  //===------------------------------------------------------------------------------------------===//
  // helper lambdas to readily construct atlas fields and views on one line
  //===------------------------------------------------------------------------------------------===//
  auto MakeAtlasField = [&](const std::string& name,
                            int size) -> std::tuple<atlas::Field, atlasInterface::Field<double>> {
    atlas::Field field_F{name, atlas::array::DataType::real64(),
                         atlas::array::make_shape(size, k_size)};
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
  //  in the ICON stencil the velocity is reconstructed at vertices (from edges)
  //  for this test, we simply assign an analytical function
  //===------------------------------------------------------------------------------------------===//
  auto [u_F, u] = MakeAtlasField("u", mesh.nodes().size());
  auto [v_F, v] = MakeAtlasField("v", mesh.nodes().size());

  //===------------------------------------------------------------------------------------------===//
  // input field (field we want to take the laplacian of)
  //  normal velocity on edges
  //===------------------------------------------------------------------------------------------===//
  auto [vn_F, vn] = MakeAtlasField("vn", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // output fields (kh_smag(1|2) are "helper" fields to store intermediary results)
  //===------------------------------------------------------------------------------------------===//
  auto [nabla2_F, nabla2] = MakeAtlasField("nabla2", mesh.edges().size());
  auto [kh_smag_1_F, kh_smag_1] = MakeAtlasField("kh_smag_1", mesh.edges().size());
  auto [kh_smag_2_F, kh_smag_2] = MakeAtlasField("kh_smag_2", mesh.edges().size());
  auto [kh_smag_F, kh_smag] = MakeAtlasField("kh_smag", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // control field
  //===------------------------------------------------------------------------------------------===//
  auto [nabla2_sol_F, nabla2_sol] = MakeAtlasField("nabla2", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // geometrical quantities on edges (vert_vert_lenght is distance between far vertices of diamond)
  //===------------------------------------------------------------------------------------------===//
  auto [inv_primal_edge_length_F, inv_primal_edge_length] =
      MakeAtlasField("inv_primal_edge_length", mesh.edges().size());
  auto [inv_vert_vert_length_F, inv_vert_vert_length] =
      MakeAtlasField("inv_vert_vert_length", mesh.edges().size());
  auto [tangent_orientation_F, tangent_orientation] =
      MakeAtlasField("tangent_orientation", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // smagorinsky coefficient stored on edges (=1 for us, simply there to force the same number of
  // reads in both ICON and our version)
  //===------------------------------------------------------------------------------------------===//
  auto [diff_multfac_smag_F, diff_multfac_smag] =
      MakeAtlasField("diff_multfac_smag", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // tangential and normal components for smagorinsky diffusion
  //===------------------------------------------------------------------------------------------===//
  auto [dvt_norm_F, dvt_norm] = MakeAtlasField("dvt_norm", mesh.edges().size());
  auto [dvt_tang_F, dvt_tang] = MakeAtlasField("dvt_tang", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // primal and dual normals at vertices (!)
  //  supposedly simply a copy of the edge normal in planar geometry (to be checked)
  //===------------------------------------------------------------------------------------------===//
  auto [primal_normal_x_F, primal_normal_x] =
      MakeAtlasSparseField("primal_normal_x", mesh.nodes().size(), verticesInDiamond);
  auto [primal_normal_y_F, primal_normal_y] =
      MakeAtlasSparseField("primal_normal_y", mesh.nodes().size(), verticesInDiamond);
  auto [dual_normal_x_F, dual_normal_x] =
      MakeAtlasSparseField("dual_normal_x", mesh.nodes().size(), verticesInDiamond);
  auto [dual_normal_y_F, dual_normal_y] =
      MakeAtlasSparseField("dual_normal_y", mesh.nodes().size(), verticesInDiamond);

  //===------------------------------------------------------------------------------------------===//
  // sparse dimension intermediary field for diamond
  //===------------------------------------------------------------------------------------------===//
  auto [vn_vert_F, vn_vert] =
      MakeAtlasSparseField("vn_vert", mesh.nodes().size(), verticesInDiamond);

  //===------------------------------------------------------------------------------------------===//
  // input (spherical harmonics) and analytical solutions for div, curl and Laplacian
  //===------------------------------------------------------------------------------------------===//

  auto sphericalHarmonic = [](double x, double y) -> std::tuple<double, double> {
    return {0.25 * sqrt(105. / (2 * M_PI)) * cos(2 * x) * cos(y) * cos(y) * sin(y),
            0.5 * sqrt(15. / (2 * M_PI)) * cos(x) * cos(y) * sin(y)};
  };
  auto analyticalLaplacian = [](double x, double y) -> std::tuple<double, double> {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return {-4 * c1 * cos(2 * x) * cos(y) * cos(y) * sin(y), -4 * c2 * cos(x) * sin(y) * cos(y)};
  };

  for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    auto [x, y] = wrapper.nodeLocation(nodeIdx);
    auto [ui, vi] = sphericalHarmonic(x, y);
    u(nodeIdx, level) = ui;
    v(nodeIdx, level) = vi;
  }
  for(int edgeIdx = 0; edgeIdx < mesh.nodes().size(); edgeIdx++) {
    auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);
    auto [x, y] = wrapper.edgeMidpoint(mesh, edgeIdx);
    auto [ui, vi] = sphericalHarmonic(x, y);
    vn(edgeIdx, level) = ui * nx + vi * ny;
  }

  for(int edgeIdx = 0; edgeIdx < mesh.nodes().size(); edgeIdx++) {
    auto [x, y] = wrapper.edgeMidpoint(mesh, edgeIdx);
    auto [ui, vi] = analyticalLaplacian(x, y);
    auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);
    nabla2_sol(edgeIdx, level) = ui * nx + vi * ny;
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//

  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    inv_primal_edge_length(edgeIdx, level) = 1. / wrapper.edgeLength(mesh, edgeIdx);
    double vert_vert_length = wrapper.vertVertLength(mesh, edgeIdx);
    inv_vert_vert_length(edgeIdx, level) = (vert_vert_length == 0) ? 0 : 1. / vert_vert_length;
    tangent_orientation(edgeIdx, level) = wrapper.tangentOrientation(mesh, edgeIdx);
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize sparse geometrical info
  //===------------------------------------------------------------------------------------------===//
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);

    for(int nbhIdx = 0; nbhIdx < verticesInDiamond; nbhIdx++) {
      primal_normal_x(edgeIdx, nbhIdx, level) = nx;
      primal_normal_y(edgeIdx, nbhIdx, level) = ny;
      dual_normal_x(edgeIdx, nbhIdx, level) = ny;
      dual_normal_y(edgeIdx, nbhIdx, level) = -nx;
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // smagorinsky fac (dummy)
  //===------------------------------------------------------------------------------------------===//
  for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    diff_multfac_smag(edgeIdx, level) = 1.;
  }
  dawn_generated::cxxnaiveico::ICON_laplacian_diamond_stencil<atlasInterface::atlasTag>(
      mesh, k_size, diff_multfac_smag, tangent_orientation, inv_primal_edge_length,
      inv_vert_vert_length, u, v, primal_normal_x, primal_normal_y, dual_normal_x, dual_normal_y,
      vn_vert, vn, dvt_tang, dvt_norm, kh_smag_1, kh_smag_2, kh_smag, nabla2)
      .run();
}