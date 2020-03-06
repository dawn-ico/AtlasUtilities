#include "AtlasToNetcdf.h"

#include <netcdf>
#include <numeric>
#include <vector>

#include <atlas/array.h>
#include <atlas/grid.h>
#include <atlas/meshgenerator.h>
#include <atlas/util/CoordinateEnums.h>

namespace {
template <std::size_t SizeT>
void writeField(const netCDF::NcFile& dataFile, const std::string& name,
                const atlas::array::ArrayView<double, SizeT> field, int fieldSize, int offset,
                std::optional<std::function<double(const double&)>> trafo = std::nullopt) {
  auto dim = dataFile.addDim("n" + name, fieldSize);
  netCDF::NcVar data = dataFile.addVar(name.c_str(), netCDF::ncDouble, dim);
  std::vector<double> dataOut;
  for(int i = 0; i < fieldSize; i++) {
    (trafo == std::nullopt) ? dataOut.push_back(field(i, offset))
                            : dataOut.push_back(trafo.value()(field(i, offset)));
  }
  data.putVar(dataOut.data());
}

template <typename ConnectivityT>
void writeConnectivityTable(const netCDF::NcFile& dataFile, const std::string& name,
                            const ConnectivityT& connectivity, int numEl, int numNbhPerEl) {
  auto dimX = dataFile.addDim("numEl" + name, numEl);
  auto dimYperX = dataFile.addDim("numNbh" + name, numNbhPerEl);
  netCDF::NcVar data = dataFile.addVar(name.c_str(), netCDF::ncInt, {dimYperX, dimX});
  std::vector<int> dataOut(numEl * numNbhPerEl);
  for(int elemIdx = 0; elemIdx < numEl; elemIdx++) {
    for(int innerIdx = 0; innerIdx < numNbhPerEl; innerIdx++) {
      // indices in netcdf are 1 based, data is column major
      int dbg = connectivity(elemIdx, innerIdx);
      dataOut[innerIdx * numEl + elemIdx] = connectivity(elemIdx, innerIdx) + 1;
    }
  }
  // printf("%d %d %d\n", dataOut[0 * numEl + 0], dataOut[1 * numEl + 0], dataOut[2 * numEl + 0]);
  data.putVar(dataOut.data());
}

bool isMinimalMesh(const atlas::Mesh& mesh) {
  return mesh.cells().edge_connectivity().rows() == 0 &&
         mesh.edges().node_connectivity().rows() == 0 &&
         mesh.edges().cell_connectivity().rows() == 0 &&
         mesh.nodes().cell_connectivity().rows() == 0 &&
         mesh.nodes().edge_connectivity().rows() == 0;
}

} // namespace

bool AtlasToNetCDF(const atlas::Mesh& mesh, const std::string& filename) {
  try {
    netCDF::NcFile dataFile(filename.c_str(), netCDF::NcFile::replace);

    auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
    auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());

    writeField<2>(dataFile, "x", xy, mesh.nodes().size(), atlas::LON);
    writeField<2>(dataFile, "y", xy, mesh.nodes().size(), atlas::LAT);
    auto lonToRad = [](double lon) { return lon * (M_PI) / 180; };
    auto latToRad = [](double lat) { return lat * (0.5 * M_PI) / 90; };
    writeField<2>(dataFile, "vlon", lonlat, mesh.nodes().size(), atlas::LON, lonToRad);
    writeField<2>(dataFile, "vlat", lonlat, mesh.nodes().size(), atlas::LAT, latToRad);

    // keeping DWD standard for naming
    const int nodesPerCell = 3;
    writeConnectivityTable(dataFile, "vertex_of_cell", mesh.cells().node_connectivity(),
                           mesh.cells().size(), nodesPerCell);

    // the minimal mesh is now written, either we are now done, or we are supposed to write the
    // complete mesh. that is, if any of the additional neighbor lists (which is both present in dwd
    // netcdf and atlas) is present, we expect _all_ of them to be present

    if(isMinimalMesh(mesh)) {
      return true;
    }

    const int edgesPerCell = 3;
    if(mesh.cells().edge_connectivity().rows() == 0) {
      std::cout << "WARNING: Partially complete mesh written to netcdf!\n";
    }
    writeConnectivityTable(dataFile, "edge_of_cell", mesh.cells().edge_connectivity(),
                           mesh.cells().size(), nodesPerCell);

    // edges
    const int cellPerEdge = 2;
    if(mesh.edges().cell_connectivity().rows() == 0) {
      std::cout << "WARNING: Partially complete mesh written to netcdf!\n";
    }
    writeConnectivityTable(dataFile, "adjacent_cell_of_edge", mesh.edges().cell_connectivity(),
                           mesh.edges().size(), cellPerEdge);

    const int nodesPerEdge = 2;
    if(mesh.edges().node_connectivity().rows() == 0) {
      std::cout << "WARNING: Partially complete mesh written to netcdf!\n";
    }
    writeConnectivityTable(dataFile, "edge_vertices", mesh.edges().node_connectivity(),
                           mesh.edges().size(), nodesPerEdge);

    // lets emulate the edge_idx field found in the DWD base grids
    {
      auto dim = dataFile.addDim("nEdgeIdx", mesh.edges().size());
      netCDF::NcVar data = dataFile.addVar("edge_index", netCDF::ncInt, dim);
      std::vector<int> dataOut(mesh.edges().size());
      std::iota(std::begin(dataOut), std::end(dataOut), 0);
      data.putVar(dataOut.data());
    }

    // nodes
    const int edgesPerNode = 6;
    if(mesh.nodes().edge_connectivity().rows() == 0) {
      std::cout << "WARNING: Partially complete mesh written to netcdf!\n";
    }
    writeConnectivityTable(dataFile, "edges_of_vertex", mesh.nodes().edge_connectivity(),
                           mesh.nodes().size(), edgesPerNode);

    const int cellsPerNode = 6;
    if(mesh.nodes().cell_connectivity().rows() == 0) {
      std::cout << "WARNING: Partially complete mesh written to netcdf!\n";
    }
    writeConnectivityTable(dataFile, "cells_of_vertex", mesh.nodes().cell_connectivity(),
                           mesh.nodes().size(), cellsPerNode);

    return true;
  } catch(netCDF::exceptions::NcException& e) {
    e.what();
    return false;
  }
}