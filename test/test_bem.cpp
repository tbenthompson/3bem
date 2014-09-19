#include "UnitTest++.h"
#include "bem.h"
#include "numerics.h"
#include <math.h>
#include <iostream>

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

TEST(RefineMesh) {
    Mesh m = square_mesh();
    Mesh refined = refine_mesh(m, {0, 2});
    double correct[2][2] = {{0.5, 0.0}, {0.5, 1.0}};
    CHECK_ARRAY_CLOSE(refined.vertices[4], correct[0], 2, 1e-14);
    CHECK_ARRAY_CLOSE(refined.vertices[5], correct[1], 2, 1e-14);

    int new_segs[4][2] = {{0, 4}, {4, 1}, {2, 5}, {5, 3}};
    CHECK_ARRAY_EQUAL(refined.segments[0], new_segs[0], 2);
    CHECK_ARRAY_EQUAL(refined.segments[1], new_segs[1], 2);
    CHECK_ARRAY_EQUAL(refined.segments[3], new_segs[2], 2);
    CHECK_ARRAY_EQUAL(refined.segments[4], new_segs[3], 2);
}

void check_mesh_to_points(QuadratureRule quad) {
    Mesh m = square_mesh();
    // m = refine_mesh(m, {0, 1, 2, 3});

    auto subsegs = get_src_obs(m, quad);
    CHECK(subsegs.center.size() == quad.size() * 4);
    for (int i = 0; i < 5; i++) {
        CHECK_CLOSE(subsegs.center[i][0], ((quad[i].first + 1.0) / 2.0), 1e-14);
        CHECK_CLOSE(subsegs.center[i][1], 0.0, 1e-14);
    }
    // TODO: Good tests for successful subsegmentation?
    // for (auto s: subsegs.center) {
    //     std::cout << "X: " << s[0] << " Y:" << s[1] << std::endl;
    // }
}
TEST(MeshToPointsGauss) {
    check_mesh_to_points(gauss(5));
}

TEST(MeshToPointsTanhSinh) {
    check_mesh_to_points(double_exp(8, 0.3));
}
