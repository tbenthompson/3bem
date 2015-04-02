#ifndef __BKJSLDFLJHLERW_MESH_GEN_H
#define __BKJSLDFLJHLERW_MESH_GEN_H
#include "mesh.h"
#include "vec_ops.h"

namespace tbem {

Mesh<3> sphere_mesh(const Vec3<double>& center, double radius, int refinements);
Mesh<3> rect_mesh(const Vec3<double>& lower_left, const Vec3<double>& upper_left,
    const Vec3<double>& upper_right, const Vec3<double>& lower_right);

Mesh<2> line_mesh(const Vec2<double>& a, const Vec2<double>& b);
//TODO: use Vec2<double> for circle_mesh center
Mesh<2> circle_mesh(std::array<double,2> center, double r, int refinements);

} //END NAMESPACE tbem

#endif
