#ifndef __BKJSLDFLJHLERW_MESH_GEN_H
#define __BKJSLDFLJHLERW_MESH_GEN_H
#include "mesh.h"

Mesh cube_mesh();
Mesh sphere_mesh(std::array<double,3> center, double radius, bool interior = true);

#endif
