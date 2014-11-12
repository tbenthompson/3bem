#ifndef __BKJSLDFLJHLERW_MESH_GEN_H
#define __BKJSLDFLJHLERW_MESH_GEN_H
#include "mesh.h"

Mesh<3> cube_mesh();
Mesh<3> sphere_mesh(const Vec3<double>& center,
                           double radius, bool interior = true);
Mesh<3> rect_mesh(const Vec3<double>& lower_left,
               const Vec3<double>& upper_left, 
               const Vec3<double>& upper_right, 
               const Vec3<double>& lower_right);

Mesh<2> square_mesh();
Mesh<2> circle_mesh(std::array<double,2> center, double r);

#endif
