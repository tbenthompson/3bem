#include "UnitTest++.h"
#include "bem.h"
#include "numerics.h"
#include "test_shared.h"
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

Mesh refined_square_mesh(int levels) {
    Mesh m = square_mesh();
    for (int i = 0; i < levels; i++) {
        m = refine_mesh(m, naturals(m.segments.size()));
    }
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

    auto subsegs = get_src_obs(m, quad);
    /* CHECK(subsegs.center.size() == quad.size() * 4); */
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

TEST(HeavyLoadMeshRefine) {
    auto m = refined_square_mesh(17);
    
    TIC
    m = refine_mesh(m, naturals(m.segments.size()));
    TOC("refine_mesh to " + std::to_string(m.segments.size()) + " segments");
}

TEST(HeavyLoadMeshToPoints) {
    auto m = refined_square_mesh(1);

    auto quad = gauss(5);
    TIC
    auto subsegs = get_src_obs(m, quad);
    TOC("get_src_obs on " + std::to_string(m.segments.size()) + " segments");
}

// TEST(InteractOneKernel) {
//     for (int refinement = 0; refinement < 5; refinement++) {
//         for (int q_order = 0; q_order < 3; q_order++) {
//             for (int n_near_steps = 1; n_near_steps < 4; n_near_steps++) {
//                 auto m = refined_square_mesh(refinement); 
//                 int n_verts = m.vertices.size();
//                 std::vector<double> src_str(n_verts);
//                 for (int i = 0; i < n_verts; i++) {
//                     src_str[i] = 1.0;
//                 }
//                 auto obs = get_src_obs(m, double_exp(q_order, 0.3));
//                 auto srcs = get_src_obs(m, gauss(q_order + 1));
//                 auto obs_values = direct_interact(m, srcs, obs, src_str, n_near_steps);
//                 for (auto o: obs_values) {
//                     CHECK_CLOSE(o, 4.0, 1e-8);
//                 }
//             }
//         }
//     }
// }

TEST(Interact) {
    auto m = refined_square_mesh(0); 
    auto obs = get_src_obs(m, double_exp(0, 0.3));
    auto srcs = get_src_obs(m, gauss(3));

    int n_verts = m.vertices.size();
    std::vector<double> src_str(n_verts);
    for (int i = 0; i < n_verts; i++) {
        src_str[i] = 1.0;
    }

    TIC
    auto obs_values = direct_interact(m, srcs, obs, src_str, 6);
    TOC("direct_interact on " + std::to_string(m.segments.size()) + " segments");
    for(auto o: obs_values) 
    {
        std::cout << o << std::endl;
    }
    //TODO: Codify the perimeter tests for the boundary element summation using
    // the one kernel.
}
