#include "mesh_gen.h"
#include "vec.h"

Mesh cube_mesh() {
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

    Mesh cube{vertices, faces, false, nullptr};
    cube = clean_mesh(cube);

    return cube;
}

Mesh sphere_mesh(Vec3<double> center, double r, bool interior) {
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

    Mesh octahedron{vertices, faces, true, 
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