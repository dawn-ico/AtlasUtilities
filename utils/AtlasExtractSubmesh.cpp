#include "AtlasExtractSubmesh.h"

#include <atlas/util/CoordinateEnums.h>

#include <atlas/mesh/ElementType.h>
#include <atlas/mesh/Elements.h>
#include <atlas/mesh/HybridElements.h>
#include <atlas/mesh/Nodes.h>

#include <numeric>
#include <optional>

namespace {
template <typename ConnectivityT>
void AllocNbhTable(ConnectivityT& connectivity, int numElements, int nbhPerElem) {
  std::vector<int> init(numElements * nbhPerElem, connectivity.missing_value());
  connectivity.add(numElements, nbhPerElem, init.data());
}

template <typename ConnectivityT>
void CopyNeighborTable(const ConnectivityT& connIn, const std::vector<int> keptIndices,
                       const std::unordered_map<int, int> oldToNewMap, ConnectivityT& connOut,
                       std::optional<std::set<int>> filter = std::nullopt) {
  int linearIdx = 0;
  for(int keptIdx : keptIndices) {
    int numNbh = connIn.cols(keptIdx);
    std::vector<int> newNbh;
    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int elemIdx = connIn(keptIdx, nbhIdx);
      if(filter != std::nullopt && filter.value().count(elemIdx) == 0) {
        newNbh.push_back(connOut.missing_value());
      } else {
        newNbh.push_back(oldToNewMap.at(elemIdx));
      }
    }
    connOut.set(linearIdx, newNbh.data());
    linearIdx++;
  }
}

void copyNodeData(const atlas::Mesh& meshIn, const std::vector<int> keptNodeIndices,
                  atlas::Mesh& meshOut) {
  const atlas::mesh::Nodes& nodesIn = meshIn.nodes();
  auto xyIn = atlas::array::make_view<double, 2>(nodesIn.xy());
  auto lonlatIn = atlas::array::make_view<double, 2>(nodesIn.lonlat());
  auto glbIdxNodeIn = atlas::array::make_view<atlas::gidx_t, 1>(nodesIn.global_index());
  auto remoteIdxIn = atlas::array::make_indexview<atlas::idx_t, 1>(nodesIn.remote_index());
  auto partIn = atlas::array::make_view<int, 1>(nodesIn.partition());
  auto ghostIn = atlas::array::make_view<int, 1>(nodesIn.ghost());
  auto flagsIn = atlas::array::make_view<int, 1>(nodesIn.flags());

  atlas::mesh::Nodes& nodes = meshOut.nodes();
  auto xy = atlas::array::make_view<double, 2>(nodes.xy());
  auto lonlat = atlas::array::make_view<double, 2>(nodes.lonlat());
  auto glbIdxNode = atlas::array::make_view<atlas::gidx_t, 1>(nodes.global_index());
  auto remoteIdx = atlas::array::make_indexview<atlas::idx_t, 1>(nodes.remote_index());
  auto part = atlas::array::make_view<int, 1>(nodes.partition());
  auto ghost = atlas::array::make_view<int, 1>(nodes.ghost());
  auto flags = atlas::array::make_view<int, 1>(nodes.flags());

  for(int nodeIdx = 0; nodeIdx < keptNodeIndices.size(); nodeIdx++) {
    xy(nodeIdx, atlas::LON) = xyIn(keptNodeIndices[nodeIdx], atlas::LON);
    xy(nodeIdx, atlas::LAT) = xyIn(keptNodeIndices[nodeIdx], atlas::LAT);
    lonlat(nodeIdx, atlas::LON) = lonlatIn(keptNodeIndices[nodeIdx], atlas::LON);
    lonlat(nodeIdx, atlas::LAT) = lonlatIn(keptNodeIndices[nodeIdx], atlas::LAT);
    glbIdxNode(nodeIdx) = glbIdxNodeIn(keptNodeIndices[nodeIdx]);
    remoteIdx(nodeIdx) = remoteIdxIn(keptNodeIndices[nodeIdx]);
    part(nodeIdx) = partIn(keptNodeIndices[nodeIdx]);
    ghost(nodeIdx) = ghostIn(keptNodeIndices[nodeIdx]);
    flags(nodeIdx) = flagsIn(keptNodeIndices[nodeIdx]);
  }
}

atlas::Mesh AtlasExtractSubMeshImpl(const atlas::Mesh& meshIn,
                                    const std::vector<int> keptCellIndices, bool complete = true) {

  // load old nbh tables
  const auto& cellToNodeIn = meshIn.cells().node_connectivity();
  const auto& cellToEdgeIn = meshIn.cells().edge_connectivity();
  const auto& edgeToNodeIn = meshIn.edges().node_connectivity();
  const auto& edgeToCellIn = meshIn.edges().cell_connectivity();
  const auto& nodeToEdgeIn = meshIn.nodes().edge_connectivity();
  const auto& nodeToCellIn = meshIn.nodes().cell_connectivity();

  // prepare maps
  // -------------

  // nodes
  std::set<int> keptNodeSet;
  for(int cellIdx : keptCellIndices) {
    int numNbh = cellToNodeIn.cols(cellIdx);
    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int nodeIdx = cellToNodeIn(cellIdx, nbhIdx);
      keptNodeSet.insert(nodeIdx);
    }
  }

  std::vector<int> keptNodeIndices(keptNodeSet.begin(), keptNodeSet.end());
  std::unordered_map<int, int> oldToNewNodeMap;
  for(int idx = 0; idx < keptNodeIndices.size(); idx++) {
    oldToNewNodeMap.emplace(keptNodeIndices[idx], idx);
  }

  // cells
  std::unordered_map<int, int> oldToNewCellMap;
  std::set<int> keptCellSet;
  for(int idx = 0; idx < keptCellIndices.size(); idx++) {
    oldToNewCellMap.emplace(keptCellIndices[idx], idx);
    keptCellSet.insert(keptCellIndices[idx]);
  }

  const int newSizeNodes = keptNodeIndices.size();
  const int newSizeCells = keptCellIndices.size();

  atlas::Mesh mesh;
  mesh.nodes().resize(newSizeNodes);
  mesh.cells().add(new atlas::mesh::temporary::Triangle(), newSizeCells);

  // copy node data
  copyNodeData(meshIn, keptNodeIndices, mesh);

  // cell data & cell nbh tables
  auto cellsPart = atlas::array::make_view<int, 1>(mesh.cells().partition());
  atlas::array::ArrayView<atlas::gidx_t, 1> glbIdxCell =
      atlas::array::make_view<atlas::gidx_t, 1>(mesh.cells().global_index());

  for(int cellIdx = 0; cellIdx < newSizeCells; cellIdx++) {
    glbIdxCell[cellIdx] = cellIdx;
    cellsPart[cellIdx] = cellIdx;
  }

  CopyNeighborTable(cellToNodeIn, keptCellIndices, oldToNewNodeMap,
                    mesh.cells().node_connectivity());

  // minimal mesh is now done
  if(!complete) {
    return mesh;
  }

  // edges (a minimal mesh may not contain edges, so we have to pull this down)
  std::set<int> keptEdgeSet;
  for(int cellIdx : keptCellIndices) {
    int numNbh = cellToEdgeIn.cols(cellIdx);
    for(int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = cellToEdgeIn(cellIdx, nbhIdx);
      keptEdgeSet.insert(edgeIdx);
    }
  }

  std::vector<int> keptEdgeIndices(keptEdgeSet.begin(), keptEdgeSet.end());
  std::unordered_map<int, int> oldToNewEdgeMap;
  for(int idx = 0; idx < keptEdgeIndices.size(); idx++) {
    oldToNewEdgeMap.emplace(keptEdgeIndices[idx], idx);
  }
  const int newSizeEdges = keptEdgeIndices.size();
  mesh.edges().add(new atlas::mesh::temporary::Line(), newSizeEdges);

  auto& cellToEdge = mesh.cells().edge_connectivity();

  const int nodesPerEdge = 2;
  const int cellsPerEdge = 2;
  const int cellsPerNode = 6; // maximum is 6, some with 5 exist
  const int edgesPerNode = 6; // maximum is 6, some with 5 exist
  const int edgesPerCell = 3;

  AllocNbhTable<atlas::mesh::HybridElements::Connectivity>(mesh.cells().edge_connectivity(),
                                                           mesh.cells().size(), edgesPerCell);
  CopyNeighborTable(cellToEdgeIn, keptCellIndices, oldToNewEdgeMap,
                    mesh.cells().edge_connectivity());

  // edge nbh tables
  AllocNbhTable<atlas::mesh::HybridElements::Connectivity>(mesh.edges().cell_connectivity(),
                                                           mesh.edges().size(), cellsPerEdge);
  AllocNbhTable<atlas::mesh::HybridElements::Connectivity>(mesh.edges().node_connectivity(),
                                                           mesh.edges().size(), nodesPerEdge);
  CopyNeighborTable(edgeToNodeIn, keptEdgeIndices, oldToNewNodeMap,
                    mesh.edges().node_connectivity());
  CopyNeighborTable(edgeToCellIn, keptEdgeIndices, oldToNewCellMap,
                    mesh.edges().cell_connectivity(), keptCellSet);

  // node nbh tables
  AllocNbhTable<atlas::mesh::Nodes::Connectivity>(mesh.nodes().cell_connectivity(),
                                                  mesh.nodes().size(), cellsPerNode);
  AllocNbhTable<atlas::mesh::Nodes::Connectivity>(mesh.nodes().edge_connectivity(),
                                                  mesh.nodes().size(), edgesPerNode);
  auto& nodeToEdge = mesh.nodes().edge_connectivity();
  auto& nodeToCell = mesh.nodes().cell_connectivity();

  CopyNeighborTable(nodeToEdgeIn, keptNodeIndices, oldToNewEdgeMap,
                    mesh.nodes().edge_connectivity(), keptEdgeSet);
  CopyNeighborTable(nodeToCellIn, keptNodeIndices, oldToNewCellMap,
                    mesh.nodes().cell_connectivity(), keptCellSet);

  return mesh;
};
} // namespace

atlas::Mesh AtlasExtractSubMeshComplete(const atlas::Mesh& meshIn, std::pair<int, int> range) {
  int lo = range.first;
  int hi = range.second;
  assert(lo >= 0 && lo < hi && hi <= meshIn.cells().size());
  const int newSizeCells = hi - lo;
  std::vector<int> keep(newSizeCells);
  std::iota(std::begin(keep), std::end(keep), lo);
  return AtlasExtractSubMeshImpl(meshIn, keep, true);
}

atlas::Mesh AtlasExtractSubMeshComplete(const atlas::Mesh& meshIn,
                                        const std::vector<int> keptCellIndices) {
  return AtlasExtractSubMeshImpl(meshIn, keptCellIndices, true);
}

atlas::Mesh AtlasExtractSubMeshMinimal(const atlas::Mesh& meshIn, std::pair<int, int> range) {
  int lo = range.first;
  int hi = range.second;
  assert(lo >= 0 && lo < hi && hi <= meshIn.cells().size());
  const int newSizeCells = hi - lo;
  std::vector<int> keep(newSizeCells);
  std::iota(std::begin(keep), std::end(keep), lo);
  return AtlasExtractSubMeshImpl(meshIn, keep, false);
}

atlas::Mesh AtlasExtractSubMeshMinimal(const atlas::Mesh& meshIn,
                                       const std::vector<int> keptCellIndices) {
  return AtlasExtractSubMeshImpl(meshIn, keptCellIndices, false);
}
