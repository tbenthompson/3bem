#ifndef __BKJSLDFLJHLERW_MESH_GEN_H
#define __BKJSLDFLJHLERW_MESH_GEN_H
#include "mesh.h"

namespace tbem {

Mesh<3> cube_mesh();
Mesh<3> sphere_mesh(const Vec3<double>& center,
                           double radius, bool interior = true);

Mesh<2> line_mesh(const Vec2<double>& a, const Vec2<double>& b);
Mesh<2> square_mesh();
Mesh<2> circle_mesh(std::array<double,2> center, double r);

} //END NAMESPACE tbem

#endif
