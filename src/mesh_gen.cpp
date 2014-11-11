#include "mesh_gen.h"
#include "vec.h"

Mesh3D cube_mesh() {
    std::vector<std::array<double, 3>> vertices = {
        {0.0,0.0,0.0}, {0.0,0.0,1.0}, {0.0,1.0,1.0},
        {1.0,1.0,0.0}, {0.0,0.0,0.0}, {0.0,1.0,0.0},
        {1.0,0.0,1.0}, {0.0,0.0,0.0}, {1.0,0.0,0.0},
        {1.0,1.0,0.0}, {1.0,0.0,0.0}, {0.0,0.0,0.0},
        {0.0,0.0,0.0}, {0.0,1.0,1.0}, {0.0,1.0,0.0},
        {1.0,0.0,1.0}, {0.0,0.0,1.0}, {0.0,0.0,0.0},
        {0.0,1.0,1.0}, {0.0,0.0,1.0}, {1.0,0.0,1.0},
        {1.0,1.0,1.0}, {1.0,0.0,0.0}, {1.0,1.0,0.0},
        {1.0,0.0,0.0}, {1.0,1.0,1.0}, {1.0,0.0,1.0},
        {1.0,1.0,1.0}, {1.0,1.0,0.0}, {0.0,1.0,0.0},
        {1.0,1.0,1.0}, {0.0,1.0,0.0}, {0.0,1.0,1.0},
        {1.0,1.0,1.0}, {0.0,1.0,1.0}, {1.0,0.0,1.0}
    };

    std::vector<std::array<int, 3>> faces;
    for (int i = 0; i < 12; i++) {
        faces.push_back({3 * i, 3 * i + 1, 3 * i + 2});
    }

    Mesh3D cube{vertices, faces, false, nullptr};
    cube = clean_mesh(cube);

    return cube;
}


Mesh3D sphere_mesh(const Vec3<double>& center, double r, bool interior) {
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
    if (!interior) {
        for (unsigned int f = 0; f < faces.size(); f++) {
            std::swap(faces[f][0], faces[f][1]);
        }
    }

    Mesh3D octahedron{vertices, faces, true, 
        [=](std::array<double,3> x) {
            double dist = std::sqrt(dist2(x, center));
            x[0] = (r / dist) * (x[0] - center[0]) + center[0];
            x[1] = (r / dist) * (x[1] - center[1]) + center[1];
            x[2] = (r / dist) * (x[2] - center[2]) + center[2];
            return x;
        }
    };

    return octahedron;
}

Mesh3D rect_mesh(const Vec3<double>& lower_left,
               const Vec3<double>& upper_left, 
               const Vec3<double>& upper_right, 
               const Vec3<double>& lower_right) {
    std::vector<Vec3<double>> vertices = {
        lower_left, upper_left, upper_right, lower_right
    };

    std::vector<std::array<int,3>> faces = {
        {0, 3, 2}, {0, 2, 1}
    };

    Mesh3D rect{vertices, faces, false, nullptr};
    return rect;
}

NewMesh<3> cube_mesh_new() {
    std::vector<std::array<double, 3>> vertices = {
        {0.0,0.0,0.0}, {0.0,0.0,1.0}, {0.0,1.0,1.0},
        {1.0,1.0,0.0}, {0.0,0.0,0.0}, {0.0,1.0,0.0},
        {1.0,0.0,1.0}, {0.0,0.0,0.0}, {1.0,0.0,0.0},
        {1.0,1.0,0.0}, {1.0,0.0,0.0}, {0.0,0.0,0.0},
        {0.0,0.0,0.0}, {0.0,1.0,1.0}, {0.0,1.0,0.0},
        {1.0,0.0,1.0}, {0.0,0.0,1.0}, {0.0,0.0,0.0},
        {0.0,1.0,1.0}, {0.0,0.0,1.0}, {1.0,0.0,1.0},
        {1.0,1.0,1.0}, {1.0,0.0,0.0}, {1.0,1.0,0.0},
        {1.0,0.0,0.0}, {1.0,1.0,1.0}, {1.0,0.0,1.0},
        {1.0,1.0,1.0}, {1.0,1.0,0.0}, {0.0,1.0,0.0},
        {1.0,1.0,1.0}, {0.0,1.0,0.0}, {0.0,1.0,1.0},
        {1.0,1.0,1.0}, {0.0,1.0,1.0}, {1.0,0.0,1.0}
    };

    std::vector<std::array<int, 3>> faces;
    for (int i = 0; i < 12; i++) {
        faces.push_back({3 * i, 3 * i + 1, 3 * i + 2});
    }

    return NewMesh<3>::from_vertices_faces(vertices, faces, false, nullptr);
}

NewMesh<3> rect_mesh_new(const Vec3<double>& lower_left,
               const Vec3<double>& upper_left, 
               const Vec3<double>& upper_right, 
               const Vec3<double>& lower_right) {
    std::vector<Vec3<double>> vertices = {
        lower_left, upper_left, upper_right, lower_right
    };

    std::vector<std::array<int,3>> faces = {
        {0, 3, 2}, {0, 2, 1}
    };

    return NewMesh<3>::from_vertices_faces(vertices, faces, false, nullptr);
}

NewMesh<3> sphere_mesh_new(const Vec3<double>& center, double r, bool interior) {
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
    if (!interior) {
        for (unsigned int f = 0; f < faces.size(); f++) {
            std::swap(faces[f][0], faces[f][1]);
        }
    }

    return NewMesh<3>::from_vertices_faces(vertices, faces, true,
        [=](std::array<double,3> x) {
            double dist = std::sqrt(dist2(x, center));
            x[0] = (r / dist) * (x[0] - center[0]) + center[0];
            x[1] = (r / dist) * (x[1] - center[1]) + center[1];
            x[2] = (r / dist) * (x[2] - center[2]) + center[2];
            return x;
        });
}

NewMesh<2> square_mesh() {
    std::vector<Vec<double,2>> vertices = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0},
    };
    
    std::vector<std::array<int, 2>> segs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };

    return NewMesh<2>::from_vertices_faces(vertices, segs, false, nullptr);
}

NewMesh<2> circle_mesh(std::array<double,2> c, double r) {
    std::vector<std::array<double,2>> vertices = {
        {c[0] + r, c[1]}, {c[0], c[1] + r}, {c[0] - r, c[1]}, {c[0], c[1] - r},
    };
    std::vector<std::array<int, 2>> segs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };
    return NewMesh<2>::from_vertices_faces(vertices, segs, true, 
        [=](const Vec2<double>& x) {
            double dist = std::sqrt(dist2(x, c));
            return Vec2<double>{
                (r / dist) * (x[0] - c[0]) + c[0],
                (r / dist) * (x[1] - c[1]) + c[1]
            };
        });
}
