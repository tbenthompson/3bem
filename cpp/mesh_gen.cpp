#include "mesh_gen.h"
#include "mesh.h"
#include "vec_ops.h"
#include "geometry.h"

namespace tbem {

Vec3<double> spherify(const Vec3<double>& center, double r, const Vec3<double>& x) {
    return (r / dist(x, center)) * (x - center) + center;
}

Mesh<3> sphere_mesh(const Vec3<double>& center, double r, int refinements) {
    std::vector<Vec3<double>> vertices =
    {
        {0.0, -r, 0.0}, {r, 0.0, 0.0}, {0.0, 0.0, r},
        {-r, 0.0, 0.0}, {0.0, 0.0, -r}, {0.0, r, 0.0}
    };

    for (unsigned int i = 0; i < vertices.size(); i++) {
        vertices[i] += center;
    }

    std::vector<std::array<int,3>> faces = 
    {
        {1, 0, 2}, {2, 0, 3}, {3, 0, 4}, {4, 0, 1},
        {5, 1, 2}, {5, 2, 3}, {5, 3, 4}, {5, 4, 1}
    };

    auto prototype = Mesh<3>::from_vertices_faces(vertices, faces)
        .refine_repeatedly(refinements);

    std::vector<Facet<3>> spherified_facets;
    for (const auto& f: prototype.facets) {
        spherified_facets.push_back({
            spherify(center, r, f[0]),
            spherify(center, r, f[1]),
            spherify(center, r, f[2])
        });
    }
    return {spherified_facets};
}
        
Mesh<3> rect_mesh(const Vec3<double>& lower_left, const Vec3<double>& upper_left,
    const Vec3<double>& upper_right, const Vec3<double>& lower_right) {

    std::vector<Vec3<double>> vertices = {
        lower_left, upper_left, upper_right, lower_right
    };
    std::vector<std::array<int,3>> faces = {
        {0, 3, 2}, {0, 2, 1}
    };
    return Mesh<3>::from_vertices_faces(vertices, faces);
}

Mesh<2> line_mesh(const Vec2<double>& a, const Vec2<double>& b) {
    std::vector<Vec<double,2>> vertices = {a, b};
    std::vector<std::array<int,2>> segs = {{0, 1}};
    return Mesh<2>::from_vertices_faces(vertices, segs);
}


Vec2<double> circlify(const Vec2<double>& center, double r, const Vec2<double>& x) {
    return (r / dist(x, center)) * (x - center) + center;
}

Mesh<2> circle_mesh(std::array<double,2> c, double r, int refinements) {
    std::vector<std::array<double,2>> vertices = {
        {c[0] + r, c[1]}, {c[0], c[1] + r}, {c[0] - r, c[1]}, {c[0], c[1] - r},
    };
    std::vector<std::array<int, 2>> segs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };

    auto prototype = Mesh<2>::from_vertices_faces(vertices, segs)
        .refine_repeatedly(refinements);

    std::vector<Facet<2>> circlified_facets;
    for (const auto& f: prototype.facets) {
        circlified_facets.push_back({
            circlify(c, r, f[0]),
            circlify(c, r, f[1])
        });
    }

    return {circlified_facets};
}

} //END namespace tbem
