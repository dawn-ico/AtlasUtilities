#pragma once

#include <atlas/mesh.h>
#include <tuple>

void generateCell2CellTable(atlas::Mesh& mesh, bool allocate);
// compute a structured index layout mesh
// Only the following connectivities are supported:
// c->n, c->e, c->c, e->c
atlas::Mesh AtlasStrIndxMesh(int ny);
