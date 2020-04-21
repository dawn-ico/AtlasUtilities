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
#include "generated_iconLaplace.hpp"

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
  //===------------------------------------------------------------------------------------------===//
  auto [u_F, u] = MakeAtlasField("u", mesh.nodes().size());
  auto [v_F, v] = MakeAtlasField("v", mesh.nodes().size());

  //===------------------------------------------------------------------------------------------===//
  // output fields
  //===------------------------------------------------------------------------------------------===//
  auto [nabla2_F, nabla2] = MakeAtlasField("nabla2", mesh.edges().size());
  auto [kh_smag_F, kh_smag] = MakeAtlasField("kh_smag", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // control field
  //===------------------------------------------------------------------------------------------===//
  auto [nabla2_sol_F, nabla2_sol] = MakeAtlasField("nabla2", mesh.edges().size());

  //===------------------------------------------------------------------------------------------===//
  // primal and dual normals at vertices (!)
  //  it is not yet entirely clear how they are defined
  //===------------------------------------------------------------------------------------------===//
  auto [primal_normal_x_F, primal_normal_x] =
      MakeAtlasField("primal_normal_x", mesh.nodes().size());
  auto [primal_normal_y_F, primal_normal_y] =
      MakeAtlasField("primal_normal_y", mesh.nodes().size());
  auto [dual_normal_x_F, dual_normal_x] = MakeAtlasField("dual_normal_x", mesh.nodes().size());
  auto [dual_normal_y_F, dual_normal_y] = MakeAtlasField("dual_normal_y", mesh.nodes().size());
}