#ifndef __BKJSLDFLJHLERW_MESH_GEN_H
#define __BKJSLDFLJHLERW_MESH_GEN_H
#include "mesh.h"

Mesh3D cube_mesh();
Mesh3D sphere_mesh(const Vec3<double>& center, double radius, bool interior = true);
Mesh3D rect_mesh(const Vec3<double>& lower_left,
               const Vec3<double>& upper_left, 
               const Vec3<double>& upper_right, 
               const Vec3<double>& lower_right);

NewMesh<3> cube_mesh_new();
NewMesh<3> sphere_mesh_new(const Vec3<double>& center,
                           double radius, bool interior = true);
NewMesh<3> rect_mesh_new(const Vec3<double>& lower_left,
               const Vec3<double>& upper_left, 
               const Vec3<double>& upper_right, 
               const Vec3<double>& lower_right);

NewMesh<2> square_mesh();
NewMesh<2> circle_mesh(std::array<double,2> center, double r);

#endif
