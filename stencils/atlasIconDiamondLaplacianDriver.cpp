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
#include "interfaces/unstructured_interface.hpp"

// io
#include "io/atlasIO.h"

std::tuple<double, double, double> MeasureErrors(std::vector<int> indices,
                                                 const atlasInterface::Field<double>& ref,
                                                 const atlasInterface::Field<double>& sol,
                                                 int level);

std::tuple<double, double, double> MeasureErrors(std::vector<int> indices,
                                                 const std::vector<double>& ref,
                                                 const atlasInterface::Field<double>& sol,
                                                 int level);

int main(int argc, char const* argv[]) {
  // enable floating point exception
  // feenableexcept(FE_INVALID | FE_OVERFLOW);

  if(argc != 2) {
    std::cout << "intended use is\n" << argv[0] << " ny" << std::endl;
    return -1;
  }
  int w = atoi(argv[1]);
  int k_size = 1;
  double lDomain = M_PI;

  // dump a whole bunch of debug output (meant to be visualized using Octave, but gnuplot and the
  // like will certainly work too)
  const bool dbg_out = true;
  const bool readMeshFromDisk = true;

  atlas::Mesh mesh;
  if(!readMeshFromDisk) {
    mesh = AtlasMeshSquare(w);
    atlas::mesh::actions::build_edges(mesh, atlas::util::Config("pole_edges", false));
    atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
    atlas::mesh::actions::build_element_to_edge_connectivity(mesh);
  } else {
    mesh = *AtlasMeshFromNetCDFComplete("grid.nc");
    auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
    auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());

    auto lonToRad = [](double rad) { return rad * (0.5 * M_PI) / 90; };
    auto latToRad = [](double rad) { return rad * (M_PI) / 180; };

    for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      xy(nodeIdx, atlas::LON) = lonToRad(lonlat(nodeIdx, atlas::LON));
      xy(nodeIdx, atlas::LAT) = latToRad(lonlat(nodeIdx, atlas::LAT));
    }
  }

  // wrapper with various atlas helper functions
  AtlasToCartesian wrapper(mesh, !readMeshFromDisk);

  const int edgesPerVertex = 6;
  const int edgesPerCell = 3;
  const int verticesInDiamond = 4;

  if(dbg_out) {
    dumpMesh4Triplot(mesh, "atlas", wrapper);
  }

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
                         atlas::array::make_shape(size, k_size, sparseSize)};
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
      MakeAtlasSparseField("primal_normal_x", mesh.edges().size(), verticesInDiamond);
  auto [primal_normal_y_F, primal_normal_y] =
      MakeAtlasSparseField("primal_normal_y", mesh.edges().size(), verticesInDiamond);
  auto [dual_normal_x_F, dual_normal_x] =
      MakeAtlasSparseField("dual_normal_x", mesh.edges().size(), verticesInDiamond);
  auto [dual_normal_y_F, dual_normal_y] =
      MakeAtlasSparseField("dual_normal_y", mesh.edges().size(), verticesInDiamond);

  //===------------------------------------------------------------------------------------------===//
  // sparse dimension intermediary field for diamond
  //===------------------------------------------------------------------------------------------===//
  auto [vn_vert_F, vn_vert] =
      MakeAtlasSparseField("vn_vert", mesh.edges().size(), verticesInDiamond);

  //===------------------------------------------------------------------------------------------===//
  // input (spherical harmonics) and analytical solutions for div, curl and Laplacian
  //===------------------------------------------------------------------------------------------===//

  auto sphericalHarmonic = [](double x, double y) -> std::tuple<double, double> {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return {c1 * cos(2 * x) * cos(y) * cos(y) * sin(y), c2 * cos(x) * cos(y) * sin(y)};
  };
  auto analyticalLaplacian = [](double x, double y) -> std::tuple<double, double> {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return {-4 * c1 * cos(2 * x) * cos(y) * cos(y) * sin(y), -4 * c2 * cos(x) * sin(y) * cos(y)};
  };
  auto analyticalScalarLaplacian = [](double x, double y) -> double {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return sin(y) * (-1. / 2. * c1 * cos(2 * x) * (13 * cos(2 * y) + 9) - 5 * c2 * cos(x) * cos(y));
  };

  auto wave = [](double x, double y) -> double { return sin(x) * sin(y); };
  auto analyticalLaplacianWave = [](double x, double y) -> double { return -2 * sin(x) * sin(y); };

  for(int level = 0; level < k_size; level++) {
    for(int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      auto [x, y] = wrapper.nodeLocation(nodeIdx);
      auto [ui, vi] = sphericalHarmonic(x, y);
      u(nodeIdx, level) = ui;
      v(nodeIdx, level) = vi;
    }
  }
  for(int level = 0; level < k_size; level++) {
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);
      auto [x, y] = wrapper.edgeMidpoint(mesh, edgeIdx);
      auto [ui, vi] = sphericalHarmonic(x, y);
      // vn(edgeIdx, level) = ui * nx + vi * ny;
      vn(edgeIdx, level) = ui * 1. + vi * 1.;
    }
  }

  if(dbg_out) {
    dumpEdgeField("diamondLaplICONatlas_in.txt", mesh, wrapper, vn, 0, wrapper.innerEdges(mesh));
  }

  for(int level = 0; level < k_size; level++) {
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      auto [x, y] = wrapper.edgeMidpoint(mesh, edgeIdx);
      auto [ui, vi] = analyticalLaplacian(x, y);
      auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);
      nabla2_sol(edgeIdx, level) = analyticalScalarLaplacian(x, y);
      // nabla2_sol(edgeIdx, level) = ui * nx + vi * ny;
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  for(int level = 0; level < k_size; level++) {
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      inv_primal_edge_length(edgeIdx, level) = 1. / wrapper.edgeLength(mesh, edgeIdx);
      double vert_vert_length = wrapper.vertVertLength(mesh, edgeIdx);
      inv_vert_vert_length(edgeIdx, level) = (vert_vert_length == 0) ? 0 : 1. / vert_vert_length;
      tangent_orientation(edgeIdx, level) = wrapper.tangentOrientation(mesh, edgeIdx);
    }
  }

  if(dbg_out) {
    dumpEdgeField("diamondLaplICONatlas_edgeLength.txt", mesh, wrapper, inv_primal_edge_length, 0,
                  wrapper.innerEdges(mesh));
    dumpEdgeField("diamondLaplICONatlas_vertVertLength.txt", mesh, wrapper, inv_vert_vert_length, 0,
                  wrapper.innerEdges(mesh));
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize sparse geometrical info
  //===------------------------------------------------------------------------------------------===//
  for(int level = 0; level < k_size; level++) {
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);

      for(int nbhIdx = 0; nbhIdx < verticesInDiamond; nbhIdx++) {
        primal_normal_x(edgeIdx, nbhIdx, level) = 1.;
        primal_normal_y(edgeIdx, nbhIdx, level) = 1.;
        dual_normal_x(edgeIdx, nbhIdx, level) = 1.;
        dual_normal_y(edgeIdx, nbhIdx, level) = 1.;

        // primal_normal_x(edgeIdx, nbhIdx, level) = nx;
        // primal_normal_y(edgeIdx, nbhIdx, level) = ny;
        // dual_normal_x(edgeIdx, nbhIdx, level) = ny;
        // dual_normal_y(edgeIdx, nbhIdx, level) = -nx;
      }
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // smagorinsky fac (dummy)
  //===------------------------------------------------------------------------------------------===//
  for(int level = 0; level < k_size; level++) {
    for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      diff_multfac_smag(edgeIdx, level) = 1.;
    }
  }

  dawn_generated::cxxnaiveico::ICON_laplacian_diamond_stencil<atlasInterface::atlasTag>(
      mesh, k_size, diff_multfac_smag, tangent_orientation, inv_primal_edge_length,
      inv_vert_vert_length, u, v, primal_normal_x, primal_normal_y, dual_normal_x, dual_normal_y,
      vn_vert, vn, dvt_tang, dvt_norm, kh_smag_1, kh_smag_2, kh_smag, nabla2)
      .run();

  //===------------------------------------------------------------------------------------------===//
  // dumping a hopefully nice colorful laplacian
  //===------------------------------------------------------------------------------------------===//
  dumpEdgeField("diamondLaplICONatlas_out.txt", mesh, wrapper, nabla2, 0, wrapper.innerEdges(mesh));
  dumpEdgeField("diamondLaplICONatlas_sol.txt", mesh, wrapper, nabla2_sol, 0,
                wrapper.innerEdges(mesh));

  {
    FILE* fp = fopen("kh_smag_ref.txt", "w+");
    for(int level = 0; level < k_size; level++) {
      for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
        fprintf(fp, "%f ", kh_smag(edgeIdx, level));
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  //===------------------------------------------------------------------------------------------===//
  // measuring errors
  //===------------------------------------------------------------------------------------------===//

  // the current test case, or rather, the analytical solution (to the current test case) is only
  // valid if the equilat triangles are aligned with the x/y axis. this is the case for the atlas
  // mesh generator, but not necessarily for imported netcdf meshes. for netcdf meshes the current
  // best bet is to simply compare against a manual computation of the same FD laplacian
  //

  if(!readMeshFromDisk) {
    for(int i = 0; i < k_size; i++) {
      auto [Linf, L1, L2] = MeasureErrors(wrapper.innerEdges(mesh), nabla2_sol, nabla2, i);
      // printf("[lap] dx: %e L_inf: %e L_1: %e L_2: %e\n", 180. / w, Linf, L1, L2);
      printf("%e %e %e %e\n", 180. / w, Linf, L1, L2);
    }
  } else {
    std::vector<double> diamondManual(mesh.edges().size());
    for(int level = 0; level < k_size; level++) {
      for(int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
        auto diamondNbh = atlasInterface::getNeighbors(
            atlasInterface::atlasTag{}, mesh,
            {dawn::LocationType::Edges, dawn::LocationType::Cells, dawn::LocationType::Vertices},
            edgeIdx);

        if(diamondNbh.size() < 4) {
          diamondManual[edgeIdx] = 0.;
        }

        std::vector<double> diamondVals = {u(diamondNbh[0], level) + v(diamondNbh[0], level),
                                           u(diamondNbh[1], level) + v(diamondNbh[1], level),
                                           u(diamondNbh[2], level) + v(diamondNbh[2], level),
                                           u(diamondNbh[3], level) + v(diamondNbh[3], level)};
        double hx = 0.5 * wrapper.edgeLength(mesh, edgeIdx);
        double hy = 0.5 * wrapper.vertVertLength(mesh, edgeIdx);

        double fxx = (diamondVals[0] + diamondVals[1] - 2 * (vn(edgeIdx, level))) / (hx * hx);
        double fyy = (diamondVals[2] + diamondVals[3] - 2 * (vn(edgeIdx, level))) / (hy * hy);

        diamondManual[edgeIdx] = fxx + fyy;
      }
    }
    for(int i = 0; i < k_size; i++) {
      auto [Linf, L1, L2] = MeasureErrors(wrapper.innerEdges(mesh), diamondManual, nabla2, i);
      printf("%e %e %e %e\n", 180. / w, Linf, L1, L2);
    }
  }

  return 0;
}

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

std::tuple<double, double, double> MeasureErrors(std::vector<int> indices,
                                                 const std::vector<double>& ref,
                                                 const atlasInterface::Field<double>& sol,
                                                 int level) {
  double Linf = 0.;
  double L1 = 0.;
  double L2 = 0.;
  for(int idx : indices) {
    double dif = ref[idx] - sol(idx, level);
    Linf = fmax(fabs(dif), Linf);
    L1 += fabs(dif);
    L2 += dif * dif;
  }
  L1 /= indices.size();
  L2 = sqrt(L2) / sqrt(indices.size());
  return {Linf, L1, L2};
}