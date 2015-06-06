#ifndef TBEMBKJSLDFLJHLERW_MESH_GEN_H
#define TBEMBKJSLDFLJHLERW_MESH_GEN_H
#include "vec.h"

namespace tbem {

template <size_t dim> struct Mesh;

Mesh<3> sphere_mesh(const Vec3<double>& center, double radius, int refinements);
Mesh<3> rect_mesh(const Vec3<double>& lower_left, const Vec3<double>& upper_left,
    const Vec3<double>& upper_right, const Vec3<double>& lower_right);

Mesh<2> line_mesh(const Vec2<double>& a, const Vec2<double>& b);
Mesh<2> circle_mesh(const Vec2<double>& center, double r, int refinements);

} //END NAMESPACE tbem

#endif
