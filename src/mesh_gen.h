#ifndef __BKJSLDFLJHLERW_MESH_GEN_H
#define __BKJSLDFLJHLERW_MESH_GEN_H
#include "mesh.h"

Mesh cube_mesh();
Mesh sphere_mesh(const Vec3<double>& center, double radius, bool interior = true);
Mesh rect_mesh(const Vec3<double>& lower_left,
               const Vec3<double>& upper_left, 
               const Vec3<double>& upper_right, 
               const Vec3<double>& lower_right);

#endif
