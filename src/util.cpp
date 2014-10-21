#include "util.h"
#include "bem.h"

Mesh square_mesh() {
    std::vector<std::array<double, 2>> vertices = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0},
    };
    
    std::vector<std::array<int, 2>> segs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };

    Mesh m = {vertices, segs};
    return m;
}

Mesh refined_square_mesh(int levels) {

    Mesh m = square_mesh();
    for (int i = 0; i < levels; i++) {
        m = refine_mesh(m, naturals(m.segments.size()));
    }
    return m;
}

Mesh circle_mesh(std::array<double, 2> center, double r, int n_segments) {
    std::vector<std::array<double, 2>> vertices(n_segments);
    std::vector<std::array<int, 2>> segs(n_segments);

    for (int i = 0; i < n_segments; i++) {
        double theta = (2 * M_PI) * (i / (double)n_segments);
        vertices[i] = {center[0] + r * cos(theta), center[1] + r * sin(theta)};
        segs[i] = {i, i + 1};
    }
    segs[n_segments - 1][1] = 0;

    Mesh m = {vertices, segs};
    return m;
}
