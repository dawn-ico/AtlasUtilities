#include "GenerateRectAtlasMesh.h"
#include <atlas/array.h>
#include <atlas/grid.h>
#include <atlas/meshgenerator.h>
#include <atlas/output/Gmsh.h>
#include <atlas/util/CoordinateEnums.h>

// determines if a cells is a downward triangle in the Atlas original indexing layout
bool cellIsDownwardTriangleOrig(const atlas::Mesh& mesh, int cell_idx) {
  auto ncols = mesh.cells().cell_connectivity().cols(cell_idx);
  for(int jc = 0; jc < ncols; ++jc) {
    auto neigh_idx = mesh.cells().cell_connectivity()(cell_idx, jc);
    if(neigh_idx < 0)
      continue;

    if(neigh_idx > cell_idx + 1) {
      return true;
    } else if(neigh_idx < cell_idx - 1) {
      return false;
    }
  }

  // If we are here is because we could not find a neighbour cell with offset > +-1
  // It can be that it is at the boundary. In this case we can determine its
  // triangle orientation based on the neighbour triangle orientation
  ncols = mesh.cells().cell_connectivity().cols(cell_idx);
  for(int jc = 0; jc < ncols; ++jc) {
    auto neigh_idx = mesh.cells().cell_connectivity()(cell_idx, jc);
    if(neigh_idx != (cell_idx - 1) && neigh_idx != (cell_idx + 1))
      continue;

    if(cellIsDownwardTriangleOrig(mesh, neigh_idx)) {
      return false;
    } else {
      return true;
    }
  }
  throw std::runtime_error("longest stride neighbour not found");
  return false;
}
// determines if a cells is a downward triangle in the structured indexing layout
bool cellIsDownwardTriangleStr(const atlas::Mesh& mesh, int cell_idx, int ncols) {
  return (cell_idx % ncols) >= ncols / 2;
}

bool cellIsDownwardTriangle(const atlas::Mesh& mesh, int cell_idx, int ncols,
                            bool structuredLayout) {
  if(structuredLayout) {
    return cellIsDownwardTriangleStr(mesh, cell_idx, ncols);
  } else {
    return cellIsDownwardTriangleOrig(mesh, cell_idx);
  }
}

// index replacement of cells in the cell connectivity
void replaceCellConnectivity(atlas::Mesh const& inmesh, atlas::Mesh& outmesh, int cidx, int ncidx) {
  auto ncols = inmesh.cells().node_connectivity().cols(cidx);
  for(int i = 0; i < ncols; ++i) {
    auto conn_value = inmesh.cells().node_connectivity()(cidx, i);

    outmesh.cells().node_connectivity().set(ncidx, i, conn_value);
  }
}

// compute the idx of a cell using a structured indexing layout
// col_idx is the column index considering that eah triangle
// (downward and upward) is a new column
int getStructuredIdxCell(int row_idx, int col_idx, int ncols) {

  // color 0 for upward triangles
  int color = -1;
  if(row_idx % 2 == 0) {
    color = ((col_idx + 1) % 2 == 0) ? 0 : 1;
  } else {
    color = ((col_idx) % 2 == 0) ? 0 : 1;
  }
  // new cell index offset wrt first cell in the row
  int ncidx_norm = col_idx / 2 + color * ncols / 2;

  return ncidx_norm + row_idx * ncols;
}

int computeEdgeNeighOfCell(atlas::Mesh& mesh, int row_idx, int col_idx, int color, int ncols) {

  auto cellIdx = getStructuredIdxCell(row_idx, col_idx, ncols);

  // leftmost edge of left most cell of rows that start with an upward triangle
  if(col_idx == 0 && (!cellIsDownwardTriangle(mesh, cellIdx, ncols, true)) && color == 0)
    return -1;
  // rightmost edge of right most cell of rows that end with an upward triangle
  if(col_idx == ncols - 1 && !cellIsDownwardTriangle(mesh, cellIdx, ncols, true) && color == 1)
    return -1;
  // bottom edges of domain
  if(row_idx == 0 && !cellIsDownwardTriangle(mesh, cellIdx, ncols, true) && color == 2)
    return -1;

  if(!cellIsDownwardTriangle(mesh, cellIdx, ncols, true)) {
    // upward triangle
    if(color == 0)
      return computeEdgeNeighOfCell(mesh, row_idx, col_idx - 1, 2, ncols);
    else if(color == 1)
      return computeEdgeNeighOfCell(mesh, row_idx, col_idx + 1, 0, ncols);
    else if(color == 2)
      return computeEdgeNeighOfCell(mesh, row_idx - 1, col_idx, 1, ncols);
  }
  auto j = col_idx / 2;
  auto njcols = ncols / 2;

  return color * njcols + j + row_idx * njcols * 3;
}

std::tuple<int, int> computeCartCoordStr(int cidx, int nx) {
  auto row_idx = cidx / nx;
  // cell index offset wrt first cell in the row
  int offset = 0;
  if(row_idx % 2 == 0) {
    if(cidx % nx < nx / 2) {
      offset = 1;
    }
  } else {
    if(cidx % nx >= nx / 2) {
      offset = 1;
    }
  }
  auto col = (cidx % (nx / 2)) * 2 + offset;

  return {row_idx, col};
}

void generateCell2CellTable(atlas::Mesh& mesh) {
  auto& cell2cell = mesh.cells().cell_connectivity();

  cell2cell.add(mesh.cells().size(), 3);

  auto& cell2edgei = mesh.cells().edge_connectivity();
  if(cell2edgei.blocks() != 1) {
    throw std::runtime_error("number of blocks for cell2edge should be one");
  }
  auto& cell2edge = cell2edgei.block(0);
  if(cell2edge.cols() != 3) {
    throw std::runtime_error("number of edge neighbours of a cell should be ==3");
  }

  auto& edge2celli = mesh.edges().cell_connectivity();
  if(edge2celli.blocks() != 1) {
    throw std::runtime_error("number of blocks for edge2cell should be one");
  }
  auto& edge2cell = edge2celli.block(0);
  if(edge2cell.cols() != 2) {
    throw std::runtime_error("number of cell neighbours of an edge should be ==2");
  }

  for(int cidx = 0; cidx < mesh.cells().size(); ++cidx) {
    // Construct here the cell to cell

    int c2c_idx = 0;
    for(int ec = 0; ec != cell2edge.cols(); ++ec) {
      int eidx = cell2edge(cidx, ec);
      for(int ce = 0; ce != edge2cell.cols(); ++ce) {
        int neigh_cell_idx = edge2cell(eidx, ce);
        if(neigh_cell_idx == cidx) {
          continue;
        }
        cell2cell.set(cidx, c2c_idx, neigh_cell_idx);
        c2c_idx++;
      }
    }
  }
}

atlas::Mesh AtlasStrIndxMesh(int nx, int ny) {

  auto mesh = AtlasMeshRect(nx, ny);
  generateCell2CellTable(mesh, true);

  auto meshstr = AtlasMeshRect(nx, ny);
  generateCell2CellTable(meshstr, true);

  atlas::mesh::Cells& cells = meshstr.cells();
  // we check some basic properties of the mesh generated
  if(nx % 2 != 0) {
    throw std::runtime_error("nx is required to be even integer");
  }
  if(meshstr.cells().edge_connectivity().blocks() != 1) {
    throw std::runtime_error("c->e connectivy contains more than one block");
  }

  if(meshstr.edges().cell_connectivity().blocks() != 1) {
    throw std::runtime_error("e->c connectivy contains more than one block");
  }

  if(meshstr.cells().cell_connectivity().blocks() != 1) {
    throw std::runtime_error("c->c connectivy contains more than one block");
  }

  if(meshstr.cells().node_connectivity().blocks() != 1) {
    throw std::runtime_error("c->n connectivy contains more than one block");
  }

  // nx is number of upward and downward triangles in a row
  if(cells.size() != nx * ny * 2) {
    throw std::runtime_error("number of cells do not match a structured layout");
  }

  // number of total edges has an excess (compared to the expected ny*ny*3) of
  // 1 edge per upward triangle at the bottom row and one edge per row (rightmost or leftmost)
  // in order to have complete cells
  if(meshstr.edges().size() != nx * ny * 3 + nx + ny) {
    throw std::runtime_error("number of edges do not match a structured layout");
  }
  if(meshstr.nodes().size() != (nx + 1) * (ny + 1)) {
    throw std::runtime_error("number of nodes do not match a structured layout");
  }

  // re-do the cell indexing
  // compute cell->node, cell->edge, cell->cell
  for(int cidx = 0; cidx != mesh.cells().size(); ++cidx) {
    auto row_idx = cidx / nx;
    // cell index offset wrt first cell in the row
    auto cidx_norm = cidx % nx;
    // new cell index offset wrt first cell in the row
    int ncidx = getStructuredIdxCell(row_idx, cidx_norm, nx);

    replaceCellConnectivity(mesh, meshstr, cidx, ncidx);

    auto& ocell2cell = meshstr.cells().cell_connectivity();
    bool isDownward = cellIsDownwardTriangle(mesh, cidx, nx, false);

    if(isDownward) {
      ocell2cell.set(ncidx, 0,
                     (row_idx == ny - 1) ? -1 : getStructuredIdxCell(row_idx + 1, cidx_norm, nx));
      ocell2cell.set(ncidx, 1,
                     (cidx_norm == 0) ? -1 : getStructuredIdxCell(row_idx, cidx_norm - 1, nx));
      ocell2cell.set(ncidx, 2,
                     (cidx_norm == nx - 1) ? -1 : getStructuredIdxCell(row_idx, cidx_norm + 1, nx));
    } else {
      ocell2cell.set(ncidx, 0,
                     (cidx_norm == 0) ? -1 : getStructuredIdxCell(row_idx, cidx_norm - 1, nx));
      ocell2cell.set(ncidx, 1,
                     (cidx_norm == nx - 1) ? -1 : getStructuredIdxCell(row_idx, cidx_norm + 1, nx));
      ocell2cell.set(ncidx, 2,
                     (row_idx == 0) ? -1 : getStructuredIdxCell(row_idx - 1, cidx_norm, nx));
    }
  }

  auto& cell2edges = meshstr.cells().edge_connectivity();
  auto& edge2cells = meshstr.edges().cell_connectivity();

  // we initialize the connectivity since the algorithm
  // below sets the e->c when there is a neighbour cell
  // but does not reset to -1 when there is no neighbour
  for(int eidx = 0; eidx != meshstr.edges().size(); ++eidx) {
    for(int h = 0; h < 2; ++h) {
      meshstr.edges().cell_connectivity().set(eidx, h, -1);
    }
  }
  // re-do the edge indexing
  // re-compute the cell -> edge and edge -> cell connectivity
  int ghost = ny * (nx / 2) * 3;
  for(int cidx = 0; cidx != meshstr.cells().size(); ++cidx) {
    bool isDownward = cellIsDownwardTriangle(meshstr, cidx, nx, true);

    auto row_idx = cidx / nx;
    // cell index offset wrt first cell in the row
    int offset = 0;
    if(row_idx % 2 == 0) {
      if(cidx % nx < nx / 2) {
        offset = 1;
      }
    } else {
      if(cidx % nx >= nx / 2) {
        offset = 1;
      }
    }
    auto col = (cidx % (nx / 2)) * 2 + offset;

    for(int h = 0; h < 3; ++h) {
      auto edgeIdx = computeEdgeNeighOfCell(meshstr, row_idx, col, h, nx);
      if(edgeIdx == -1)
        edgeIdx = ghost++;

      cell2edges.set(cidx, h, edgeIdx);
      edge2cells.set(edgeIdx, isDownward ? 0 : 1, cidx);
    }
  }
  return meshstr;
}